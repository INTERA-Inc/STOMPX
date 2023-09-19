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
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     This engineering program numerically simulates the transport
!     of H2O, NaCl, CH4 and CO2 through multifluid subsurface 
!     environments under non-isothermal conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      ISUB_LOG = 1
      SUB_LOG(1) = 'STOMP-HYD'
      ICODE = 37
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
      CALL RDINPT37
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
      CALL CHK37
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
      IF( ISLC(50).NE.0 ) CALL JCBP_GM
!
!---  Compute primary variable increments  ---
!
      CALL INCRM37
!
!---  Initial hydrologic and thermodynamic properties on interior nodes
!     and boundary surfaces  ---
!
      CALL PROP37
      CALL BCP37
!
!---  Compute initial solute concentrations  ---
!
      CALL CISC37
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
!---  Compute initial fluxes on interior and boundary surfaces  ---
!
      ISVF = 1
      CALL DRCVL
      CALL DRCVG
!
!---  Iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DRCVN
!
!---  Iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFGA
      CALL DFFGO
      CALL DFFGW
!
!---  Iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFLA
      CALL DFFLO
!
!---  Isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS
      CALL DFFLW
!
!---  Iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFNO
      IF( ISLC(45).EQ.0 ) CALL DFFNW
      IF( ISLC(45).EQ.0 ) CALL DFFNA
!
!---  Isothermal option  ---
!
      IF( ISLC(30).EQ.0 ) THEN
        CALL THD37
        CALL THL37
        CALL THG37
!
!---    Iso-CO2 option  ---
!
        IF( ISLC(45).EQ.0 )  CALL THN37
!        CALL THDL37
!        CALL THDG37
      ENDIF
      CALL BCF37
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
      CALL LDO37
!
!---  Load old time step arrays for the volumetric stress
!     and pore pressure  ---
!
      IF( ISLC(50).NE.0 ) THEN
        INDX = 1
        CALL LD_GM( INDX )
      ENDIF
!
!---  Reference node(s) output  ---
!
      IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 .OR.
     &  MOD( (NSTEP-NRST),IFQO ).EQ.0 ) CALL REFNOD
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
        CALL BCP37
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
!---  Recalculate zero-iterate phases and properties for kinetic
!     guest molecule exchange  ---
!
      IF( NITER.EQ.1 ) THEN
        CALL INCRM37
        CALL PROP37
      ENDIF
!
!---  Compute boundary saturation, relative permeability, and
!     thermodynamic properties  ---
!
      CALL BCP37
!
!---  Compute source contributions  ---
!
      CALL SORC37
!
!---  Compute aqueous-phase volumetric flux (interior surfaces)  ---
!
      CALL DRCVL
!
!---  Compute gas-phase volumetric flux (interior surfaces)  ---
!
      CALL DRCVG
!
!---  Compute liquid-CO2 phase volumetric flux (interior surfaces),
!     iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DRCVN
!
!---  Compute water vapor diffusion flux through the gas phase
!     (interior surfaces), iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFGA
!
!---  Compute CH4 diffusion flux through the gas phase
!     (interior surfaces), iso-CO2 option  ---
!
      CALL DFFGO
      CALL DFFGW
!
!---  Compute dissolved-CO2 diffusion flux through the aqueous phase
!     (interior surfaces), iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFLA
!
!---  Compute dissolved-CH4 diffusion flux through the aqueous phase
!     (interior surfaces), iso-CO2 option  ---
!
      CALL DFFLO
!
!---  Compute dissolved-NaCl diffusion flux through the aqueous phase
!     (interior surfaces), isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS
      CALL DFFLW
!
!---  Compute CH4 diffusion flux through the liquid-CO2 phase
!     (interior surfaces), iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFNO
!
!---  Compute H2O diffusion flux through the liquid-CO2 phase
!     (interior surfaces), iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFNW
!
!---  Compute CO2 diffusion flux through the liquid-CO2 phase
!     (interior surfaces), iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFNA
!
!---  Isothermal option  ---
!
      IF( ISLC(30).EQ.0 ) THEN
!
!---    Compute diffusive thermal flux (interior surfaces)  ---
!
        CALL THD37
!
!---    Compute aqueous-phase advective thermal flux
!       (interior surfaces)  ---
!
        CALL THL37
!
!---    Compute gas-phase advective thermal flux
!       (interior surfaces)  ---
!
        CALL THG37
!
!---    Compute liquid-CO2-phase advective thermal flux
!       (interior surfaces), iso-CO2 option  ---
!
        IF( ISLC(45).EQ.0 ) CALL THN37
!
!---    Compute aqueous diffusive thermal flux
!       (interior surfaces)  ---
!
!        CALL THDL37
!
!---    Compute gas diffusive thermal flux
!       (interior surfaces)  ---
!
!        CALL THDG37
      ENDIF
!
!---  Compute aqueous-phase volumetric flux, gas-phase volumetric flux,
!     water vapor mass flux, aqueous-phase salt flux
!     (boundary surfaces)  ---
!
      CALL BCF37
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
!     (zero flux boundary), isothermal option  ---
!
      IF( ISLC(30).EQ.0 ) CALL JCBT_HYD
!
!---  Load Jacobian matrix for the water equation
!     (zero flux boundary)  ---
!
      CALL JCBW_HYD
!
!---  Load Jacobian matrix for the CO2 equation
!     (zero flux boundary), iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL JCBA_HYD
!
!---  Load Jacobian matrix for the CH4 equation
!     (zero flux boundary)  ---
!
      CALL JCBO_HYD
!
!---  Load Jacobian matrix for the salt equation
!     (zero flux boundary), isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL JCBS_HYD
!
!---  Modify the Jacobian matrix for boundary conditions  ---
!
      CALL BCJ37
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
      CALL UPDT37
!
!---  Compute convergence from maximum relative residuals  ---
!
      CALL RSDL37
!
!---  Compute primary variable increments, saturation,
!     relative permeability, porosity, tortuosity,
!     thermodynamic properties for interior nodes,
!     except immediately after a new time step  ---
!
      CALL INCRM37
      CALL PROP37
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
!---  Compute current fluxes for transport solutions or flux
!     integrations  ---
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
!---  Compute liquid-CO2 phase volumetric flux (interior surfaces),
!     iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DRCVN
!
!---  Compute water vapor diffusion flux through the gas phase
!     (interior surfaces), iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFGA
!
!---  Compute CH4 diffusion flux through the gas phase
!     (interior surfaces), iso-CO2 option  ---
!
      CALL DFFGO
      CALL DFFGW
!
!---  Compute dissolved-CO2 diffusion flux through the aqueous phase
!     (interior surfaces), iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFLA
!
!---  Compute dissolved-CH4 diffusion flux through the aqueous phase
!     (interior surfaces), iso-CO2 option  ---
!
      CALL DFFLO
!
!---  Compute dissolved-NaCl diffusion flux through the aqueous phase
!     (interior surfaces), isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS
      CALL DFFLW
!
!---  Compute CH4 diffusion flux through the liquid-CO2 phase
!     (interior surfaces), iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFNO
!
!---  Compute H2O diffusion flux through the liquid-CO2 phase
!     (interior surfaces), iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFNW
!
!---  Compute CO2 diffusion flux through the liquid-CO2 phase
!     (interior surfaces), iso-CO2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFNA
!
!---  Isothermal option  ---
!
      IF( ISLC(30).EQ.0 ) THEN
!
!---    Compute diffusive thermal flux (interior surfaces)  ---
!
        CALL THD37
!
!---    Compute aqueous-phase advective thermal flux
!       (interior surfaces)  ---
!
        CALL THL37
!
!---    Compute gas-phase advective thermal flux
!       (interior surfaces)  ---
!
        CALL THG37
!
!---    Compute liquid-CO2-phase advective thermal flux
!       (interior surfaces), iso-CO2 option  ---
!
        IF( ISLC(45).EQ.0 ) CALL THN37
!
!---    Compute aqueous diffusive thermal flux
!       (interior surfaces)  ---
!
!        CALL THDL37
!
!---    Compute gas diffusive thermal flux
!       (interior surfaces)  ---
!
!        CALL THDG37
      ENDIF
!
!---  Compute aqueous-phase volumetric flux, gas-phase volumetric flux,
!     water vapor mass flux, aqueous-phase salt flux
!     (boundary surfaces)  ---
!
      CALL BCF37
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
          CALL SPRP37( NSL )
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
          CALL SORT37( NSL )
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
          CALL SBND37( NSL )
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
          CALL SORIT37( NSL )
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
      SUBROUTINE BCF37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Compute boundary surface fluxes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      USE FLUXN
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/BCF37'
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
            WN(M,NPZ) = 0.D+0
            WLA(M,NPZ) = 0.D+0
            WLO(M,NPZ) = 0.D+0
            WLW(M,NPZ) = 0.D+0
            WGA(M,NPZ) = 0.D+0
            WGO(M,NPZ) = 0.D+0
            WGW(M,NPZ) = 0.D+0
            WNC(1,M,NPZ) = 0.D+0
            WNC(2,M,NPZ) = 0.D+0
            WNC(3,M,NPZ) = 0.D+0
            WDLA(M,NPZ) = 0.D+0
            WDLO(M,NPZ) = 0.D+0
            WDGA(M,NPZ) = 0.D+0
            WDGO(M,NPZ) = 0.D+0
            WDGW(M,NPZ) = 0.D+0
            WDNC(1,M,NPZ) = 0.D+0
            WDNC(2,M,NPZ) = 0.D+0
            WDNC(3,M,NPZ) = 0.D+0
            WQ(M,NPZ) = 0.D+0
            WS(M,NPZ) = 0.D+0
            WDS(M,NPZ) = 0.D+0
   10     CONTINUE
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          DO 20 M = 1,ISVF
            VL(M,NPY) = 0.D+0
            VG(M,NPY) = 0.D+0
            VN(M,NPY) = 0.D+0
            VLA(M,NPY) = 0.D+0
            VLO(M,NPY) = 0.D+0
            VLW(M,NPY) = 0.D+0
            VGA(M,NPY) = 0.D+0
            VGO(M,NPY) = 0.D+0
            VGW(M,NPY) = 0.D+0
            VNC(1,M,NPY) = 0.D+0
            VNC(2,M,NPY) = 0.D+0
            VNC(3,M,NPY) = 0.D+0
            VDLA(M,NPY) = 0.D+0
            VDLO(M,NPY) = 0.D+0
            VDGA(M,NPY) = 0.D+0
            VDGO(M,NPY) = 0.D+0
            VDGW(M,NPY) = 0.D+0
            VDNC(1,M,NPY) = 0.D+0
            VDNC(2,M,NPY) = 0.D+0
            VDNC(3,M,NPY) = 0.D+0
            VQ(M,NPY) = 0.D+0
            VS(M,NPY) = 0.D+0
            VDS(M,NPY) = 0.D+0
   20     CONTINUE
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          DO 30 M = 1,ISVF
            UL(M,NPX) = 0.D+0
            UG(M,NPX) = 0.D+0
            UN(M,NPX) = 0.D+0
            ULA(M,NPX) = 0.D+0
            ULO(M,NPX) = 0.D+0
            ULW(M,NPX) = 0.D+0
            UGA(M,NPX) = 0.D+0
            UGO(M,NPX) = 0.D+0
            UGW(M,NPX) = 0.D+0
            UNC(1,M,NPX) = 0.D+0
            UNC(2,M,NPX) = 0.D+0
            UNC(3,M,NPX) = 0.D+0
            UDLA(M,NPX) = 0.D+0
            UDLO(M,NPX) = 0.D+0
            UDGA(M,NPX) = 0.D+0
            UDGO(M,NPX) = 0.D+0
            UDGW(M,NPX) = 0.D+0
            UDNC(1,M,NPX) = 0.D+0
            UDNC(2,M,NPX) = 0.D+0
            UDNC(3,M,NPX) = 0.D+0
            UQ(M,NPX) = 0.D+0
            US(M,NPX) = 0.D+0
            UDS(M,NPX) = 0.D+0
   30     CONTINUE
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          DO 40 M = 1,ISVF
            UL(M,NQX) = 0.D+0
            UG(M,NQX) = 0.D+0
            UN(M,NQX) = 0.D+0
            ULA(M,NQX) = 0.D+0
            ULO(M,NQX) = 0.D+0
            ULW(M,NQX) = 0.D+0
            UGA(M,NQX) = 0.D+0
            UGO(M,NQX) = 0.D+0
            UGW(M,NQX) = 0.D+0
            UNC(1,M,NQX) = 0.D+0
            UNC(2,M,NQX) = 0.D+0
            UNC(3,M,NQX) = 0.D+0
            UDLA(M,NQX) = 0.D+0
            UDLO(M,NQX) = 0.D+0
            UDGA(M,NQX) = 0.D+0
            UDGO(M,NQX) = 0.D+0
            UDGW(M,NQX) = 0.D+0
            UDNC(1,M,NQX) = 0.D+0
            UDNC(2,M,NQX) = 0.D+0
            UDNC(3,M,NQX) = 0.D+0
            UQ(M,NQX) = 0.D+0
            US(M,NQX) = 0.D+0
            UDS(M,NQX) = 0.D+0
   40     CONTINUE
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          DO 50 M = 1,ISVF
            VL(M,NQY) = 0.D+0
            VG(M,NQY) = 0.D+0
            VN(M,NQY) = 0.D+0
            VLA(M,NQY) = 0.D+0
            VLO(M,NQY) = 0.D+0
            VLW(M,NQY) = 0.D+0
            VGA(M,NQY) = 0.D+0
            VGO(M,NQY) = 0.D+0
            VGW(M,NQY) = 0.D+0
            VNC(1,M,NQY) = 0.D+0
            VNC(2,M,NQY) = 0.D+0
            VNC(3,M,NQY) = 0.D+0
            VDLA(M,NQY) = 0.D+0
            VDLO(M,NQY) = 0.D+0
            VDGA(M,NQY) = 0.D+0
            VDGO(M,NQY) = 0.D+0
            VDGW(M,NQY) = 0.D+0
            VDNC(1,M,NQY) = 0.D+0
            VDNC(2,M,NQY) = 0.D+0
            VDNC(3,M,NQY) = 0.D+0
            VQ(M,NQY) = 0.D+0
            VS(M,NQY) = 0.D+0
            VDS(M,NQY) = 0.D+0
   50     CONTINUE
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          DO 60 M = 1,ISVF
            WL(M,NQZ) = 0.D+0
            WG(M,NQZ) = 0.D+0
            WN(M,NQZ) = 0.D+0
            WLA(M,NQZ) = 0.D+0
            WLO(M,NQZ) = 0.D+0
            WLW(M,NQZ) = 0.D+0
            WGA(M,NQZ) = 0.D+0
            WGO(M,NQZ) = 0.D+0
            WGW(M,NQZ) = 0.D+0
            WNC(1,M,NQZ) = 0.D+0
            WNC(2,M,NQZ) = 0.D+0
            WNC(3,M,NQZ) = 0.D+0
            WDLA(M,NQZ) = 0.D+0
            WDLO(M,NQZ) = 0.D+0
            WDGA(M,NQZ) = 0.D+0
            WDGO(M,NQZ) = 0.D+0
            WDGW(M,NQZ) = 0.D+0
            WDNC(1,M,NQZ) = 0.D+0
            WDNC(2,M,NQZ) = 0.D+0
            WDNC(3,M,NQZ) = 0.D+0
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
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSB( N,NB )
            ENDIF
          ENDIF
!
!---      Aqueous Neumann boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 110 M = 1,ISVF
              WL(M,NSZ(N)) = BCX(3)
  110       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAB( N,NB )
            CALL DFFLOB( N,NB )
            CALL DFFLWB( N,NB )
!
!---      Aqueous not zero flux boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLB( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAB( N,NB )
            CALL DFFLOB( N,NB )
            CALL DFFLWB( N,NB )
          ENDIF
!
!---      Gas Neumann boundary  ---
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 111 M = 1,ISVF
              WG(M,NSZ(N)) = BCX(4)
  111       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAB( N,NB )
            CALL DFFGOB( N,NB )
            CALL DFFGWB( N,NB )
!
!---      Gas not zero flux boundary  ---
!
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVGB( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAB( N,NB )
            CALL DFFGOB( N,NB )
            CALL DFFGWB( N,NB )
          ENDIF
!
!---      Iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
!
!---        Liquid-CO2 Neumann else Dirichlet, Saturated,
!           Unit Gradient  ---
!
            IF( IBCT(4,NB).EQ.2 ) THEN
              DO 112 M = 1,ISVF
                WN(M,NSZ(N)) = BCX(11)
  112         CONTINUE
              CALL DFFNOB( N,NB )
              CALL DFFNWB( N,NB )
              CALL DFFNAB( N,NB )
            ELSEIF( IBCT(4,NB).NE.3 ) THEN
              CALL DRCVNB( N,NB )
              CALL DFFNOB( N,NB )
              CALL DFFNWB( N,NB )
              CALL DFFNAB( N,NB )
            ENDIF
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO 113 M = 1,ISVF
                WQ(M,NSZ(N)) = BCX(2)
  113         CONTINUE
!
!---        Energy outflow or inflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13 ) THEN
              CALL THLB37( N,NB )
              CALL THGB37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNB37( N,NB )
!              CALL THDLB37( N,NB )
!              CALL THDGB37( N,NB )
!
!---        Energy zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDB37( N,NB )
              CALL THLB37( N,NB )
              CALL THGB37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNB37( N,NB )
!              CALL THDLB37( N,NB )
!              CALL THDGB37( N,NB )
            ENDIF
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSS( N,NB )
            ENDIF
          ENDIF
!
!---      Aqueous Neumann boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 120 M = 1,ISVF
              VL(M,NSY(N)) = BCX(3)
  120       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAS( N,NB )
            CALL DFFLOS( N,NB )
            CALL DFFLWS( N,NB )
!
!---      Aqueous not zero flux boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLS( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAS( N,NB )
            CALL DFFLOS( N,NB )
            CALL DFFLWS( N,NB )
          ENDIF
!
!---      Gas Neumann boundary  ---
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 121 M = 1,ISVF
              VG(M,NSY(N)) = BCX(4)
  121       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAS( N,NB )
            CALL DFFGOS( N,NB )
            CALL DFFGWS( N,NB )
!
!---      Gas not zero flux boundary  ---
!
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVGS( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAS( N,NB )
            CALL DFFGOS( N,NB )
            CALL DFFGWS( N,NB )
          ENDIF
!
!---      Iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
!
!---        Liquid-CO2 Neumann boundary  ---
!
            IF( IBCT(4,NB).EQ.2 ) THEN
              DO 122 M = 1,ISVF
                VN(M,NSY(N)) = BCX(11)
  122         CONTINUE
              CALL DFFNOS( N,NB )
              CALL DFFNWS( N,NB )
              CALL DFFNAS( N,NB )
!
!---        Liquid-CO2 not zero flux boundary  ---
!
            ELSEIF( IBCT(4,NB).NE.3 ) THEN
              CALL DRCVNS( N,NB )
              CALL DFFNOS( N,NB )
              CALL DFFNWS( N,NB )
              CALL DFFNAS( N,NB )
            ENDIF
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO 123 M = 1,ISVF
                VQ(M,NSY(N)) = BCX(2)
  123         CONTINUE
!
!---        Energy outflow or inflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13 ) THEN
              CALL THLS37( N,NB )
              CALL THGS37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNS37( N,NB )
!              CALL THDLS37( N,NB )
!              CALL THDGS37( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDS37( N,NB )
              CALL THLS37( N,NB )
              CALL THGS37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNS37( N,NB )
!              CALL THDLS37( N,NB )
!              CALL THDGS37( N,NB )
            ENDIF
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSW( N,NB )
            ENDIF
          ENDIF
!
!---     Aqueous Neumann boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 130 M = 1,ISVF
              UL(M,NSX(N)) = BCX(3)
  130       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAW( N,NB )
            CALL DFFLOW( N,NB )
            CALL DFFLWW( N,NB )
!
!---      Aqueous not zero flux boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLW( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAW( N,NB )
            CALL DFFLOW( N,NB )
            CALL DFFLWW( N,NB )
          ENDIF
!
!---      Gas Neumann boundary  ---
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 131 M = 1,ISVF
              UG(M,NSX(N)) = BCX(4)
  131       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAW( N,NB )
            CALL DFFGOW( N,NB )
            CALL DFFGWW( N,NB)
!
!---      Gas not zero flux boundary  ---
!
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVGW( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAW( N,NB )
            CALL DFFGOW( N,NB )
            CALL DFFGWW( N,NB)
          ENDIF
!
!---      Iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
!
!---        Liquid-CO2 Neumann boundary  ---
!
            IF( IBCT(4,NB).EQ.2 ) THEN
              DO 132 M = 1,ISVF
                UN(M,NSX(N)) = BCX(11)
  132         CONTINUE
              CALL DFFNOW( N,NB )
              CALL DFFNWW( N,NB )
              CALL DFFNAW( N,NB )
!
!---        Liquid-CO2 not zero flux boundary  ---
!
            ELSEIF( IBCT(4,NB).NE.3 ) THEN
              CALL DRCVNW( N,NB )
              CALL DFFNOW( N,NB )
              CALL DFFNWW( N,NB )
              CALL DFFNAW( N,NB )
            ENDIF
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO 133 M = 1,ISVF
                UQ(M,NSX(N)) = BCX(2)
  133         CONTINUE
!
!---        Energy outflow or inflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13 ) THEN
              CALL THLW37( N,NB )
              CALL THGW37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNW37( N,NB )
!              CALL THDLW37( N,NB )
!              CALL THDGW37( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDW37( N,NB )
              CALL THLW37( N,NB )
              CALL THGW37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNW37( N,NB )
!              CALL THDLW37( N,NB )
!              CALL THDGW37( N,NB )
            ENDIF
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSE( N,NB )
            ENDIF
          ENDIF
!
!---      Aqueous Neumann boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 140 M = 1,ISVF
              UL(M,NSX(N)+1) = BCX(3)
  140       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAE( N,NB )
            CALL DFFLOE( N,NB )
            CALL DFFLWE( N,NB )
!
!---      Aqueous not zero flux boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLE( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAE( N,NB )
            CALL DFFLOE( N,NB )
            CALL DFFLWE( N,NB )
          ENDIF
!
!---      Gas Neumann boundary  ---
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 141 M = 1,ISVF
              UG(M,NSX(N)+1) = BCX(4)
  141       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAE( N,NB )
            CALL DFFGOE( N,NB )
            CALL DFFGWE( N,NB )
!
!---      Gas not zero flux boundary  ---
!
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVGE( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAE( N,NB )
            CALL DFFGOE( N,NB )
            CALL DFFGWE( N,NB )
          ENDIF
!
!---      Iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
!
!---        Liquid-CO2 Neumann boundary  ---
!
            IF( IBCT(4,NB).EQ.2 ) THEN
              DO 142 M = 1,ISVF
                UN(M,NSX(N)+1) = BCX(11)
  142         CONTINUE
              CALL DFFNOE( N,NB )
              CALL DFFNWE( N,NB )
              CALL DFFNAE( N,NB )
            ELSEIF( IBCT(4,NB).NE.3 ) THEN
              CALL DRCVNE( N,NB )
              CALL DFFNOE( N,NB )
              CALL DFFNWE( N,NB )
              CALL DFFNAE( N,NB )
            ENDIF
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO 143 M = 1,ISVF
                UQ(M,NSX(N)+1) = BCX(2)
  143         CONTINUE
!
!---        Energy outflow or inflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13 ) THEN
              CALL THLE37( N,NB )
              CALL THGE37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNE37( N,NB )
!              CALL THDLE37( N,NB )
!              CALL THDGE37( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDE37( N,NB )
              CALL THLE37( N,NB )
              CALL THGE37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNE37( N,NB )
!              CALL THDLE37( N,NB )
!              CALL THDGE37( N,NB )
            ENDIF
          ENDIF
!
!---    North boundary
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSN( N,NB )
            ENDIF
          ENDIF
!
!---      Aqueous Neumann boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 150 M = 1,ISVF
              VL(M,NSY(N)+IFLD) = BCX(3)
  150       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAN( N,NB )
            CALL DFFLON( N,NB )
            CALL DFFLWN( N,NB )
!
!---      Aqueous not zero flux boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLN( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAN( N,NB )
            CALL DFFLON( N,NB )
            CALL DFFLWN( N,NB )
          ENDIF
!
!---      Gas Neumann boundary  ---
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 151 M = 1,ISVF
              VG(M,NSY(N)+IFLD) = BCX(4)
  151       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAN( N,NB )
            CALL DFFGON( N,NB )
            CALL DFFGWN( N,NB )
!
!---      Gas not zero flux boundary  ---
!
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVGN( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAN( N,NB )
            CALL DFFGON( N,NB )
            CALL DFFGWN( N,NB )
          ENDIF
!
!---      Iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
!
!---        Liquid-CO2 Neumann boundary  ---
!
            IF( IBCT(4,NB).EQ.2 ) THEN
              DO 152 M = 1,ISVF
                VN(M,NSY(N)+IFLD) = BCX(11)
  152         CONTINUE
              CALL DFFNAN( N,NB )
              CALL DFFNAN( N,NB )
              CALL DFFNAN( N,NB )
!
!---        Liquid-CO2 not zero flux boundary  ---
!
            ELSEIF( IBCT(4,NB).NE.3 ) THEN
              CALL DRCVNN( N,NB )
              CALL DFFNAN( N,NB )
              CALL DFFNAN( N,NB )
              CALL DFFNAN( N,NB )
            ENDIF
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO 153 M = 1,ISVF
                VQ(M,NSY(N)+IFLD) = BCX(2)
  153         CONTINUE
!
!---        Energy outflow or inflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13 ) THEN
              CALL THLN37( N,NB )
              CALL THGN37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNN37( N,NB )
!              CALL THDLN37( N,NB )
!              CALL THDGN37( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDN37( N,NB )
              CALL THLN37( N,NB )
              CALL THGN37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNN37( N,NB )
!              CALL THDLN37( N,NB )
!              CALL THDGN37( N,NB )
            ENDIF
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLST( N,NB )
            ENDIF
          ENDIF
!
!---      Aqueous Neumann boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 160 M = 1,ISVF
              WL(M,NSZ(N)+IJFLD) = BCX(3)
  160       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAT( N,NB )
            CALL DFFLOT( N,NB )
            CALL DFFLWT( N,NB )
!
!---      Aqueous not zero flux boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLT( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFLAT( N,NB )
            CALL DFFLOT( N,NB )
            CALL DFFLWT( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 161 M = 1,ISVF
              WG(M,NSZ(N)+IJFLD) = BCX(4)
  161       CONTINUE
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAT( N,NB )
            CALL DFFGOT( N,NB )
            CALL DFFGWT( N,NB )
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVGT( N,NB )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL DFFGAT( N,NB )
            CALL DFFGOT( N,NB )
            CALL DFFGWT( N,NB )
          ENDIF
!
!---      Iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
!
!---        Liquid-CO2 Neumann boundary  ---
!
            IF( IBCT(4,NB).EQ.2 ) THEN
              DO 162 M = 1,ISVF
                WN(M,NSZ(N)+IJFLD) = BCX(11)
  162         CONTINUE
              CALL DFFNOT( N,NB )
              CALL DFFNWT( N,NB )
              CALL DFFNAT( N,NB )
!
!---        Liquid-CO2 not zero flux boundary  ---
!
            ELSEIF( IBCT(4,NB).NE.3 ) THEN
              CALL DRCVNT( N,NB )
              CALL DFFNOT( N,NB )
              CALL DFFNWT( N,NB )
              CALL DFFNAT( N,NB )
            ENDIF
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO 163 M = 1,ISVF
                WQ(M,NSZ(N)+IJFLD) = BCX(2)
  163         CONTINUE
!
!---        Energy outflow or inflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13 ) THEN
              CALL THLT37( N,NB )
              CALL THGT37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNT37( N,NB )
!              CALL THDLT37( N,NB )
!              CALL THDGT37( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDT37( N,NB )
              CALL THLT37( N,NB )
              CALL THGT37( N,NB )
!
!---          Iso-CO2 option  ---
!
              IF( ISLC(45).EQ.0 ) CALL THNT37( N,NB )
!              CALL THDLT37( N,NB )
!              CALL THDGT37( N,NB )
            ENDIF
          ENDIF
        ENDIF
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCF37 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCJ37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Modify the Jacobian matrix for boundary conditions
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/BCJ37'
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
!---      Aqueous  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLW_B( N,NB,NPZ )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBLA_B( N,NB,NPZ )
            CALL JCBLO_B( N,NB,NPZ )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBGW_B( N,NB,NPZ )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBGA_B( N,NB,NPZ )
            CALL JCBGO_B( N,NB,NPZ )
          ENDIF
!
!---      Liquid CO2  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBNW_B( N,NB,NPZ )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBNA_B( N,NB,NPZ )
            CALL JCBNO_B( N,NB,NPZ )
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBT_B( N,NB,NPZ )
            ENDIF
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBS_B( N,NB,NPZ )
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLW_S( N,NB,NPY )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBLA_S( N,NB,NPY )
            CALL JCBLO_S( N,NB,NPY )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBGW_S( N,NB,NPY )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBGA_S( N,NB,NPY )
            CALL JCBGO_S( N,NB,NPY )
          ENDIF
!
!---      Liquid CO2  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBNW_S( N,NB,NPY )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBNA_S( N,NB,NPY )
            CALL JCBNO_S( N,NB,NPY )
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBT_S( N,NB,NPY )
            ENDIF
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBS_S( N,NB,NPY )
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLW_W( N,NB,NPX )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBLA_W( N,NB,NPX )
            CALL JCBLO_W( N,NB,NPX )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBGW_W( N,NB,NPX )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBGA_W( N,NB,NPX )
            CALL JCBGO_W( N,NB,NPX )
          ENDIF
!
!---      Liquid CO2  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBNW_W( N,NB,NPX )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBNA_W( N,NB,NPX )
            CALL JCBNO_W( N,NB,NPX )
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBT_W( N,NB,NPX )
            ENDIF
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBS_W( N,NB,NPX )
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLW_E( N,NB,NQX )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBLA_E( N,NB,NQX )
            CALL JCBLO_E( N,NB,NQX )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBGW_E( N,NB,NQX )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBGA_E( N,NB,NQX )
            CALL JCBGO_E( N,NB,NQX )
          ENDIF
!
!---      Liquid CO2  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBNW_E( N,NB,NQX )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBNA_E( N,NB,NQX )
            CALL JCBNO_E( N,NB,NQX )
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBT_E( N,NB,NQX )
            ENDIF
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBS_E( N,NB,NQX )
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLW_N( N,NB,NQY )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBLA_N( N,NB,NQY )
            CALL JCBLO_N( N,NB,NQY )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBGW_N( N,NB,NQY )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBGA_N( N,NB,NQY )
            CALL JCBGO_N( N,NB,NQY )
          ENDIF
!
!---      Liquid CO2  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBNW_N( N,NB,NQY )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBNA_N( N,NB,NQY )
            CALL JCBGW_N( N,NB,NQY )
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBT_N( N,NB,NQY )
            ENDIF
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBS_N( N,NB,NQY )
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLW_T( N,NB,NQZ )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBLA_T( N,NB,NQZ )
            CALL JCBLO_T( N,NB,NQZ )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBGW_T( N,NB,NQZ )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBGA_T( N,NB,NQZ )
            CALL JCBGO_T( N,NB,NQZ )
          ENDIF
!
!---      Liquid CO2  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBNW_T( N,NB,NQZ )
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) CALL JCBNA_T( N,NB,NQZ )
            CALL JCBNO_T( N,NB,NQZ )
          ENDIF
!
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBT_T( N,NB,NQZ )
            ENDIF
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBS_T( N,NB,NQZ )
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCJ37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCP37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Compute saturation, relative permeability and thermodynamic
!     properties for boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FDVT
      USE FDVP
      USE FDVN
      USE FDVI
      USE FDVH
      USE FDVG
      USE CONST
      USE BCVT
      USE BCVS
      USE BCVP
      USE BCVN
      USE BCVI
      USE BCVH
      USE BCVGC
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
      REAL*8 BCX(LBCV),YLSX(2)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCP37'
      EPSLX = 1.D-4
!
!---  Assign values for initial condition type boundary conditions  ---
!
      IF( NSTEP-NRST.LE.1 .AND. NITER.LE.1 ) THEN
       DO 50 NB = 1,NBC
        IF( IBCT(1,NB).EQ.12 .OR. IBCT(2,NB).EQ.12 .OR.
     &    IBCT(3,NB).EQ.12 .OR. IBCT(4,NB).EQ.12 ) THEN
          N = IBCN(NB)
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
        ENDIF
        IF( IBCT(1,NB).EQ.12 ) TB(1,NB) = T(1,N)
        IF( IBCT(2,NB).EQ.12 ) PLB(1,NB) = PL(2,N) + RHOL(2,N)*GB
        IF( IBCT(3,NB).EQ.12 ) PGB(1,NB) = PG(2,N) + RHOG(2,N)*GB
        IF( IBCT(4,NB).EQ.12 ) PNB(1,NB) = PN(2,N) + RHON(2,N)*GB
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
!
!---         Energy boundary  ---
!
             IF( IBCT(1,NB).EQ.2 ) THEN
               BCX(2) = BCX(2)-5.D-1*DTBC*(BC(2,M,MB)-BC(2,M-1,MB))/TDBC
             ENDIF
!
!---         Aqueous boundary  ---
!
             IF( IBCT(2,NB).EQ.2 ) THEN
               BCX(3) = BCX(3)-5.D-1*DTBC*(BC(3,M,MB)-BC(3,M-1,MB))/TDBC
             ENDIF
!
!---         Gas boundary  ---
!
             IF( IBCT(3,NB).EQ.2 ) THEN
               BCX(4) = BCX(4)-5.D-1*DTBC*(BC(4,M,MB)-BC(4,M-1,MB))/TDBC
             ENDIF
!
!---         Liquid-CO2 boundary  ---
!
             IF( IBCT(4,NB).EQ.2 ) THEN
               BCX(11) = BCX(11)-5.D-1*DTBC*(BC(11,M,MB)-BC(11,M-1,MB))
     &           /TDBC
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
        IF( IBCT(2,NB).EQ.12 ) BCX(3) = PLB(1,NB)
        IF( IBCT(3,NB).EQ.12 ) BCX(4) = PGB(1,NB)
        IF( IBCT(4,NB).EQ.12 ) BCX(11) = PNB(1,NB)
        IF( IBCT(1,NB).EQ.12 ) BCX(2) = TB(1,NB)
        N = IBCN(NB)
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
          TX = T(2,N)
          PLX = PL(M,N)
          PGX = PG(M,N)
          PNX = PN(M,N)
!
!---      Gas Dirichlet  ---
!
          IF( IBCT(3,NB).EQ.1 .OR. IBCT(3,NB).EQ.12
     &      .OR. IBCT(3,NB).EQ.26 .OR. IBCT(3,NB).EQ.42 ) THEN
            PGX = BCX(4)
!
!---      Gas Neumann  ---
!
          ELSEIF( IBCT(3,NB).EQ.2 ) THEN
            PGX = PGX + BCX(4)*DB*VISG(M,N)/PERM(IBD,IZN)
     &        + RHOG(M,N)*GB
!
!---      Gas Zero Flux  ---
!
          ELSEIF( IBCT(3,NB).EQ.3 ) THEN
            IF( (BCX(4)+PATM).GT.-1.D+16 ) THEN
              PGX = BCX(4)
            ELSE
              PGX = PGX + RHOG(M,N)*GB
            ENDIF
!
!---      Gas Unit Gradient  ---
!
          ELSEIF( IBCT(3,NB).EQ.5 ) THEN
            PGX = PGX
!
!---      Gas Hydraulic Gradient Base Node  ---
!
          ELSEIF( IBCT(3,NB).EQ.-11 ) THEN
            PGX = BCX(4)
!
!---      Gas Hydraulic Gradient  ---
!
          ELSEIF( IBCT(3,NB).EQ.11 ) THEN
            IF( M.EQ.2 ) THEN
              PGX = HGBG37( NBHG(3,NB),NB,M )
            ELSE
              PGX = PGB(2,NB)
            ENDIF
          ENDIF
!
!---      Aqueous Dirichlet  ---
!
          IF( IBCT(2,NB).EQ.1 .OR. IBCT(2,NB).EQ.12 .OR.
     &      IBCT(2,NB).EQ.26 .OR. IBCT(2,NB).EQ.42 ) THEN
            PLX = BCX(3)
!
!---      Aqueous Neumann  ---
!
          ELSEIF( IBCT(2,NB).EQ.2 ) THEN
            PLX = PLX + BCX(3)*DB*VISL(M,N)/PERM(IBD,IZN)
     &        + RHOL(M,N)*GB
!
!---      Aqueous Zero Flux  ---
!
          ELSEIF( IBCT(2,NB).EQ.3 ) THEN
            IF( (BCX(3)+PATM).GT.-1.D+16 ) THEN
              PLX = BCX(3)
            ELSE
              PLX = PLX + RHOL(M,N)*GB
            ENDIF
!
!---      Aqueous Unit Gradient  ---
!
          ELSEIF( IBCT(2,NB).EQ.5 ) THEN
            PLX = PLX
!
!---      Aqueous Hydraulic Gradient Base Node  ---
!
          ELSEIF( IBCT(2,NB).EQ.-11 ) THEN
            PLX = BCX(3)
!
!---      Aqueous Hydraulic Gradient  ---
!
          ELSEIF( IBCT(2,NB).EQ.11 ) THEN
            IF( M.EQ.2 ) THEN
              PLX = HGBL37( NBHG(2,NB),NB,M )
            ELSE
              PLX = PLB(2,NB)
            ENDIF
          ENDIF
!
!---      Liquid-CO2 Dirichlet  ---
!
          IF( IBCT(4,NB).EQ.1 .OR. IBCT(4,NB).EQ.12 .OR.
     &      IBCT(4,NB).EQ.26 .OR. IBCT(4,NB).EQ.42 ) THEN
            PNX = BCX(11)
!
!---      Liquid-CO2 Neumann; Logic prevents liquid-CO2 
!         flow from two-phase nodes  ---
!
          ELSEIF( IBCT(4,NB).EQ.2 ) THEN
            PNX = PNX + BCX(11)*DB*VISN(M,N)/PERM(IBD,IZN)
     &        + RHON(M,N)*GB
!
!---      Liquid-CO2 Zero Flux  ---
!
          ELSEIF( IBCT(4,NB).EQ.3 ) THEN
            IF( (BCX(11)+PATM).GT.-1.D+16 ) THEN
              PNX = BCX(11)
            ELSE
              PNX = PNX + RHON(M,N)*GB
            ENDIF
!
!---      Liquid-CO2 Unit Gradient  ---
!
          ELSEIF( IBCT(4,NB).EQ.5 ) THEN
            PNX = PNX
!
!---      Liquid-CO2 Hydraulic Gradient Base Node  ---
!
          ELSEIF( IBCT(4,NB).EQ.-11 ) THEN
            PNX = BCX(11)
!
!---      Liquid-CO2 Hydraulic Gradient  ---
!
          ELSEIF( IBCT(4,NB).EQ.11 ) THEN
            IF( M.EQ.2 ) THEN
              PNX = HGBN37( NBHG(4,NB),NB,M )
            ELSE
              PNX = PNB(2,NB)
            ENDIF
          ENDIF
!
!---      Energy Dirichlet, outflow, initial conditions
!         or inflow  ---
!
          IF( IBCT(1,NB).EQ.1 .OR. IBCT(1,NB).EQ.7 
     &      .OR. IBCT(1,NB).EQ.12 .OR. IBCT(1,NB).EQ.13 ) THEN
            TX = BCX(2)
!
!---      Energy Neumann  ---
!
          ELSEIF( IBCT(1,NB).EQ.2 ) THEN
            INDX = ABS( IBCD(NB) )
            TKP = MAX(1.D+0-PORD(M,N),0.D+0)*THKS(INDX,IZ(N)) + 
     &        PORD(M,N)*(THKL(M,N)*SL(M,N) + 
     &        THKH(M,N)*SH(M,N) + THKI(M,N)*SI(M,N) + 
     &        THKN(M,N)*SN(M,N) + THKG(M,N)*SG(M,N))
            TX = TX + BCX(2)*DB/TKP
!
!---      Zero Flux  ---
!
          ELSEIF( IBCT(1,NB).EQ.3 ) THEN
            TX = TX
          ENDIF
!
!---      Salt concentration  ---
!
          CALL SOL_LS( T(2,N),XLSMX )
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).EQ.34 .OR. IBCT(IEQS,NB).EQ.35 ) THEN
              YLSB(M,NB) = BCX(6)*XLSMX
              XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
            ELSEIF( IBCT(IEQS,NB).EQ.36 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
              YLSB(M,NB) = BCX(6)
              XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
            ENDIF
          ELSE
            YLSB(M,NB) = 0.D+0
            XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
          ENDIF
!
!---      Interfacial surface tensions  ---
!
          SIG_GL = SCHR(5,IZN)
          SIG_NL = SCHR(6,IZN)
          SIG_GN = SCHR(7,IZN)
          SIG_HL = SCHR(9,IZN)
          SIG_IL = SCHR(9,IZN)
          CALL SFT_L( TX,XLSB(M,NB),SIG_GL )
          CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &      /(2.D+0*SIG_GL*SIG_NL)
          CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &      /(2.D+0*SIG_GL*SIG_GN)
!
!---      Interphase scaling factors  ---
!
          BGL = SCHR(8,IZN)/SIG_GL
          BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
          BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
          BHL = SCHR(8,IZN)/SCHR(9,IZN)
          BIL = SCHR(8,IZN)/SCHR(9,IZN)
!
!---      Hydrate saturation  ---
!
          SHB(M,NB) = BCX(12)
!
!---      Entry pressure with hydrate barrier adjustment  ---
!
          HCPBFX = 1.D+0
          IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
            HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &        SIN(SCHR(12,IZN)*SHB(M,NB) + SCHR(13,IZN)) )
          ENDIF
          ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---      System pressure  ---
!
          PX = MAX( PGX,PLX,PNX ) + PATM
!
!---      Brine vapor pressure  ---
!
          CALL SP_B( TX,XLSB(M,NB),PSBX )
          PPLX = MAX( PLX+PATM,PSBX )
          CALL DENS_B( TX,PPLX,XLSB(M,NB),RHOBX )
          PCX = MAX( PGX-PLX,0.D+0 )
          IF( ISLC(44).EQ.1 ) THEN
            PVWX = PSBX
          ELSE
            CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSB(M,NB) )
          ENDIF
          PVWB(M,NB) = PVWX*BCX(5)
!
!---      Saturated CO2 pressure  ---
!
          CALL SP_A( TX,PSAX )
!
!---      Temperature above the CO2 critical point,
!         no liquid CO2 conditions  ---
!
          TKX = TX+TABS
          IF( TKX.GT.TCRA ) THEN
!
!---        Liquid-CO2 check  ---
!
            PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),PLX+ENPR/BNL )
            IF( PNX.GT.PNC ) THEN
              INDX = 17
              N_DB = N
              RLMSG = TX
              CHMSG = 'Boundary Liquid CO2 Conditions'
     &          // ' w/ Temperature > Critical CO2 Temperature: '
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Micro-emulsion  ---
!
            IF( BCX(10).GT.EPSL ) THEN
              PGX = PLX + ENPR/BGL
              PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),PLX+ENPR/BNL )
              PNX = PNC
              YGAX = 1.D+0
            ELSE
              PGX = MAX( (PLX+ENPR/BGL),(PNX+ENPR/BGN),PGX )
              PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),PLX+ENPR/BNL )
              PNX = MAX( PNC,PNX )
              YGAX = MAX( 0.D+0,MIN( BCX(9),1.D+0 ) )
            ENDIF
!
!---        Gas mole fraction of hydrate formers  ---
!
            YGOX = MAX( 1.D+0-YGAX,0.D+0 )
            YMGAX = (YGAX/WTMA)/((YGAX/WTMA)+(YGOX/WTMO))
            YMGOX = MAX( 1.D+0-YMGAX,0.D+0 )
            PVX = PGX+PATM-PVWB(M,NB)
            PVAB(M,NB) = YMGAX*PVX
            PVOB(M,NB) = YMGOX*PVX
!
!---        Micro-emulsion  ---
!
            IF( BCX(10).GT.EPSL ) ZLAB(M,NB) = -BCX(10)
!
!---        Hydrate saturation   ---
!
            IF( SHB(M,NB).GT.EPSL ) THEN
              NPHAZX = 2
            ELSE
              NPHAZX = 1
            ENDIF
!
!---      Temperature below the CO2 critical point  ---
!
          ELSE
            PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),PLX+ENPR/BNL )
!
!---        Liquid-CO2 pressure greater than critical pressure,
!           liquid-CO2 conditions  ---
!
            IF( (PNX-PNC).GT.EPSL ) THEN
              PVX = MAX( PGX+PATM-PVWB(M,NB),PSAX+PVWB(M,NB) )
              PVAB(M,NB) = PSAX
              PVOB(M,NB) = MAX( PVX-PSAX-PVAB(M,NB),0.D+0 )
              YMGOX = PVOB(M,NB)/(PVAB(M,NB)+PVOB(M,NB))
              PGX = MAX( (PLX+ENPR/BGL),(PNX+ENPR/BGN),PGX )
!
!---          Micro-emulsion  ---
!
              IF( BCX(10).GT.EPSL ) ZLAB(M,NB) = BCX(10)
!
!---          Hydrate saturation   ---
!
              IF( SHB(M,NB).GT.EPSL ) THEN
                NPHAZX = 6
              ELSE
                NPHAZX = 5
              ENDIF
!
!---        Liquid-CO2 pressure equal to critical pressure
!           no liquid-CO2 conditions  ---
!
            ELSE
              YGAX = MAX( 0.D+0,MIN( BCX(9),1.D+0 ) )
              YGOX = MAX( 1.D+0-YGAX,0.D+0 )
              YMGAX = (YGAX/WTMA)/((YGAX/WTMA)+(YGOX/WTMO))
              YMGOX = MAX( 1.D+0-YMGAX,0.D+0 )
              PVX = PGX+PATM-PVWB(M,NB)
              PVAB(M,NB) = YMGAX*PVX
              PVOB(M,NB) = YMGOX*PVX
!
!---          Micro-emulsion  ---
!
              IF( BCX(10).GT.EPSL ) ZLAB(M,NB) = -BCX(10)
!
!---          Hydrate saturation   ---
!
              IF( SHB(M,NB).GT.EPSL ) THEN
                NPHAZX = 2
              ELSE
                NPHAZX = 1
              ENDIF
            ENDIF
          ENDIF
!
!---      Set boundary temperature and pressures  ---
!
          TB(M,NB) = TX
          PLB(M,NB) = PLX
          PGB(M,NB) = PGX
          PNB(M,NB) = PNX
!
!---      Aqueous or aqueous-gas, no hydrate, no liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - CO2 vapor pressure
!         CH4 mass - CH4 vapor pressure
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          IF( NPHAZX.EQ.1 ) THEN
            PSBX = -PSBX
            CALL FLH_PC1( PEQHX,PGB(M,NB),PHB(M,NB),PIB(M,NB),
     &        PLB(M,NB),PNB(M,NB),PSBX,PVAB(M,NB),PVOB(M,NB),
     &        PVWB(M,NB),RHOHB(M,NB),RKGB(M,NB),RKLB(1,M,NB),
     &        RKNB(M,NB),SHB(M,NB),SIB(M,NB),SGB(M,NB),
     &        SLB(M,NB),SNB(M,NB),TB(M,NB),TEQHX,XHAB(M,NB),
     &        XHOB(M,NB),XHWB(M,NB),XLSB(M,NB),YLSB(M,NB),
     &        YMGOX,IZN,N )
!
!---      Aqueous or aqueous-gas, hydrate, no liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - Gas CH4 mole fraction of formers
!         CH4 mass - hydrate saturation
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZX.EQ.2 ) THEN
            PSBX = -PSBX
#ifdef eq_hyd
            CALL FLH_PC2( PEQHX,PGB(M,NB),PHB(M,NB),PIB(M,NB),
     &        PLB(M,NB),PNB(M,NB),PSBX,PVAB(M,NB),PVOB(M,NB),
     &        PVWB(M,NB),RHOHB(M,NB),RKGB(M,NB),RKLB(1,M,NB),
     &        RKNB(M,NB),SHB(M,NB),SIB(M,NB),SGB(M,NB),
     &        SLB(M,NB),SNB(M,NB),TB(M,NB),TEQHX,XHAB(M,NB),
     &        XHOB(M,NB),XHWB(M,NB),XLSB(M,NB),YLSB(M,NB),
     &        YMGOX,IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = YMGOX
            CALL FLH_PC2( PEQHX,PGB(M,NB),PHB(M,NB),PIB(M,NB),
     &        PLB(M,NB),PNB(M,NB),PSBX,PVAB(M,NB),PVOB(M,NB),
     &        PVWB(M,NB),RHOHB(M,NB),RKGB(M,NB),RKLB(1,M,NB),
     &        RKNB(M,NB),SHB(M,NB),SIB(M,NB),SGB(M,NB),
     &        SLB(M,NB),SNB(M,NB),TB(M,NB),TEQHX,XHAB(M,NB),
     &        XHOB(M,NB),XHWB(M,NB),XLSB(M,NB),YLSB(M,NB),
     &        YMGOX,YMGOY,YMGOOX,IZN,N )
#endif
#ifdef k2_hyd
            PVAHOX = PVAB(M,NB)
            PVOHOX = PVOB(M,NB)
            CALL FLH_PC2( PEQHX,PGB(M,NB),PHB(M,NB),PIB(M,NB),
     &        PLB(M,NB),PNB(M,NB),PSBX,PVAB(M,NB),PVAHX,PVAHOX,
     &        PVOB(M,NB),PVOHX,PVOHOX,PVWB(M,NB),
     &        RHOHB(M,NB),RKGB(M,NB),RKLB(1,M,NB),
     &        RKNB(M,NB),SHB(M,NB),SIB(M,NB),SGB(M,NB),
     &        SLB(M,NB),SNB(M,NB),TB(M,NB),TEQHX,XHAB(M,NB),
     &        XHOB(M,NB),XHWB(M,NB),XLSB(M,NB),YLSB(M,NB),
     &        YMGOX,IZN,N )
#endif
!
!---      Aqueous or aqueous-gas, no hydrate, liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - liquid CO2 pressure
!         CH4 mass - CH4 vapor pressure
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZX.EQ.5 ) THEN
            PSBX = -PSBX
            CALL FLH_PC5( PEQHX,PGB(M,NB),PHB(M,NB),PIB(M,NB),
     &        PLB(M,NB),PNB(M,NB),PSBX,PVAB(M,NB),PVOB(M,NB),
     &        PVWB(M,NB),RHOHB(M,NB),RKGB(M,NB),RKLB(1,M,NB),
     &        RKNB(M,NB),SHBX,SIB(M,NB),SGB(M,NB),
     &        SLB(M,NB),SNB(M,NB),TB(M,NB),TEQHX,XHAB(M,NB),
     &        XHOB(M,NB),XHWB(M,NB),XLSB(M,NB),YLSB(M,NB),
     &        YMGOX,IZN,N )
!
!---      Aqueous or aqueous-gas, hydrate, liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - liquid CO2 pressure
!         CH4 mass - hydrate saturation
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZX.EQ.6 ) THEN
            PSBX = -PSBX
#ifdef eq_hyd
            CALL FLH_PC6( PEQHX,PGB(M,NB),PHB(M,NB),PIB(M,NB),
     &        PLB(M,NB),PNB(M,NB),PSBX,PVAB(M,NB),PVOB(M,NB),
     &        PVWB(M,NB),RHOHB(M,NB),RKGB(M,NB),RKLB(1,M,NB),
     &        RKNB(M,NB),SHBX,SIB(M,NB),SGB(M,NB),
     &        SLB(M,NB),SNB(M,NB),TB(M,NB),TEQHX,XHAB(M,NB),
     &        XHOB(M,NB),XHWB(M,NB),XLSB(M,NB),YLSB(M,NB),
     &        YMGOX,IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = YMGOX
            CALL FLH_PC6( PEQHX,PGB(M,NB),PHB(M,NB),PIB(M,NB),
     &        PLB(M,NB),PNB(M,NB),PSBX,PVAB(M,NB),PVOB(M,NB),
     &        PVWB(M,NB),RHOHB(M,NB),RKGB(M,NB),RKLB(1,M,NB),
     &        RKNB(M,NB),SHBX,SIB(M,NB),SGB(M,NB),
     &        SLB(M,NB),SNB(M,NB),TB(M,NB),TEQHX,XHAB(M,NB),
     &        XHOB(M,NB),XHWB(M,NB),XLSB(M,NB),YLSB(M,NB),
     &        YMGOX,YMGOY,YMGOOX,IZN,N )
#endif
#ifdef k2_hyd
            PVAHOX = PVAB(M,NB)
            PVOHOX = PVOB(M,NB)
            CALL FLH_PC6( PEQHX,PGB(M,NB),PHB(M,NB),PIB(M,NB),
     &        PLB(M,NB),PNB(M,NB),PSBX,PVAB(M,NB),PVAHX,PVAHOX,
     &        PVOB(M,NB),PVOHX,PVOHOX,PVWB(M,NB),
     &        RHOHB(M,NB),RKGB(M,NB),RKLB(1,M,NB),
     &        RKNB(M,NB),SHBX,SIB(M,NB),SGB(M,NB),
     &        SLB(M,NB),SNB(M,NB),TB(M,NB),TEQHX,XHAB(M,NB),
     &        XHOB(M,NB),XHWB(M,NB),XLSB(M,NB),YLSB(M,NB),
     &        YMGOX,IZN,N )
#endif
          ENDIF
!
!---      Hydrate temperature check  ---
!
          IF( SHB(M,NB).GT.EPSL .AND. PVX.LT.PEQHX ) THEN
            INDX = 17
            N_DB = N
            RLMSG = PGX+PATM
            CHMSG = 'Boundary Gas Vapor Pressure < Hydrate Equilibrium'
     &        // ' Pressure w/ Hydrate Saturation: '
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Convert pressure to absolute prior to computing physical
!         properties  ---
!
          PLX = PLB(M,NB) + PATM
          PGX = PGB(M,NB) + PATM
          PNX = PNB(M,NB) + PATM
!
!---      Dissolved-CO2 and Dissolved-CH4 aqueous mole fractions   ---
!
          INCG = 1
          PVAX = BCX(7)*PVAB(M,NB)
          CALL SOL_LA( TX,PVAX,XLSB(M,NB),XMLAB(M,NB) )
          INCG = 2
          PVOX = BCX(8)*PVOB(M,NB)
          CALL SOL_LO( TX,PVOX,XLSB(M,NB),XMLOB(M,NB) )
!
!---      Dissolved-CO2 and Dissolved-CH4 aqueous mass fractions  ---
!
          WTMLX = (XMLOB(M,NB)*WTMO + XMLAB(M,NB)*WTMA + 
     &      WTMW*(1.D+0-XMLOB(M,NB)-XMLAB(M,NB)))/
     &      (1.D+0-XLSB(M,NB)+WTMW*XLSB(M,NB)/WTMS)
          XLAB(M,NB) = XMLAB(M,NB)*WTMA/WTMLX
          XLOB(M,NB) = XMLOB(M,NB)*WTMO/WTMLX
          CALL DENS_L( TX,RHOBX,XLAB(M,NB),XLOB(M,NB),RHOLB(M,NB) )
          XLWB(M,NB) = MAX( 1.D+0-XLAB(M,NB)-XLOB(M,NB)-XLSB(M,NB),
     &      0.D+0 )
          WTMLX = XLSB(M,NB)/WTMS + XLAB(M,NB)/WTMA + 
     &      XLOB(M,NB)/WTMO + XLWB(M,NB)/WTMW
          XMLSB(M,NB) = XLSB(M,NB)/WTMS/WTMLX
          XMLAB(M,NB) = XLAB(M,NB)/WTMA/WTMLX
          XMLOB(M,NB) = XLOB(M,NB)/WTMO/WTMLX
          XMLWB(M,NB) = XLWB(M,NB)/WTMW/WTMLX
!
!---      Aqueous molar density  ---
!
          RHOMLB(M,NB) = RHOLB(M,NB)*WTMLX
!
!---      Porous-media porosity  ---
!
          CALL PORSTY( N,PX,PCMP(N),PORDB(M,NB),PORTB(M,NB) )
          PORDB(M,NB) = MAX( PORDB(M,NB),EPSL )
          PORTB(M,NB) = MAX( PORTB(M,NB),PORDB(M,NB) )
!
!---      Aqueous, gas, and liquid-CO2 tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SLB(M,NB),SGB(M,NB),
     &      SNB(M,NB),PORDB(M,NB),TORLB(M,NB),TORGB(M,NB),TORNB(M,NB) )
!
!---      Gas density and component fractions  ---
!
          ISRX = 2
          CALL DENS_W( TX,PVWB(M,NB),RHOX,RHOGWX,ISRX )
          CALL DENS_A( TX,PVAB(M,NB),RHOGAX,ISRX )
          CALL DENS_O( TX,PVOB(M,NB),RHOGOX )
          RHOGB(M,NB) = RHOGWX+RHOGAX+RHOGOX
          XGAB(M,NB) = RHOGAX/RHOGB(M,NB)
          XGOB(M,NB) = RHOGOX/RHOGB(M,NB)
          XGWB(M,NB) = RHOGWX/RHOGB(M,NB)
          WTMGX = XGAB(M,NB)/WTMA + XGOB(M,NB)/WTMO + XGWB(M,NB)/WTMW
          XMGAB(M,NB) = XGAB(M,NB)/WTMA/WTMGX
          XMGOB(M,NB) = XGOB(M,NB)/WTMO/WTMGX
          XMGWB(M,NB) = XGWB(M,NB)/WTMW/WTMGX
          RHOMGB(M,NB) = RHOGB(M,NB)*WTMGX
!
!---      Liquid-CO2 density and viscosity  ---
!
          ISRX = 1
          CALL DENS_A( TX,PX,RHONB(M,NB),ISRX )
          CALL VISC_A( TX,RHONB(M,NB),VISNB(M,NB) )
          CALL SOL_N( PX,RHOGB(M,NB),RHOMGB(M,NB),RHONB(M,NB),
     &      RHOMNB(M,NB),TX,WTMNX,XMGAB(M,NB),XMGOB(M,NB),
     &      XMGWB(M,NB),XNCB(1,M,NB),XNCB(2,M,NB),XNCB(3,M,NB) )
          XMNCB(1,M,NB) = XNCB(1,M,NB)/WTMA/WTMNX
          XMNCB(2,M,NB) = XNCB(2,M,NB)/WTMO/WTMNX
          XMNCB(3,M,NB) = XNCB(3,M,NB)/WTMW/WTMNX
          XNAB(M,NB) = XNCB(1,M,NB)
          XNOB(M,NB) = XNCB(2,M,NB)
          XNWB(M,NB) = XNCB(3,M,NB)
          CALL DIFC_N( DFNCB(1,M,NB),TX,VISNB(M,NB),WTMNX )
!
!---      Ice density  ---
!
          CALL DENS_I( TX,RHOIB(M,NB) )
!
!---      Gas viscosity  ---
!
          CALL VISC_A( TX,RHOGAX,VISGAX )
          CALL VISC_O( TX,RHOGOX,VISGOOX,VISGOEX,VISGOX )
          CALL VISC_W( TX,PVWB(M,NB),RHOGWX,VISGWX )
          CALL VISC_G( VISGAX,VISGOX,VISGWX,XMGAB(M,NB),
     &      XMGOB(M,NB),XMGWB(M,NB),VISGB(M,NB) )
!
!---      Gas diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            DFGOB(M,NB) = DFGOC
            DFGAB(M,NB) = DFGAC
          ELSEIF( ISLC(2).EQ.2 ) THEN
            CALL BNDFAW( TX,PGX,DFGAW )
            CALL BNDFAO( TX,PGX,DFGAO )
            CALL BNDFOW( TX,PGX,DFGOW )
            IF( (1.D+0-XMGOB(M,NB)).LT.EPSL ) THEN
              DFGOB(M,NB) = 2.D+0*DFGAO*DFGOW/(DFGAO+DFGOW)
            ELSE
              DFGOB(M,NB) = (1.D+0-XMGOB(M,NB))/
     &          ((XMGAB(M,NB)/DFGAO) + (XMGWB(M,NB)/DFGOW))
            ENDIF
            IF( (1.D+0-XMGAB(M,NB)).LT.EPSL ) THEN
              DFGAB(M,NB) = 2.D+0*DFGAO*DFGAW/(DFGAO+DFGAW)
            ELSE
              DFGAB(M,NB) = (1.D+0-XMGAB(M,NB))/
     &          ((XMGOB(M,NB)/DFGAO) + (XMGWB(M,NB)/DFGAW))
            ENDIF
          ENDIF
!
!---      Aqueous viscosity  ---
!
          ISRX = 1
          CALL DENS_W( TX,PX,RHOLWX,RHOX,ISRX )
          CALL VISC_W( TX,PX,RHOLWX,VISLWX )
          CALL VISC_B( TX,XLSB(M,NB),VISLWX,VISBX )
          CALL VISC_L( XMLAB(M,NB),XMLOB(M,NB),VISBX,VISGAX,VISGOX,
     &      VISLB(M,NB) )
!
!---      Dissolved-CO2 and -NaCl diffusion coefficient  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLAB(M,NB) = DFLAC
            DFLOB(M,NB) = DFLOC
            DFLSB(M,NB) = DFLSC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            INCG = 1
            CALL DIFC_LA( TX,VISLB(M,NB),VISGAX,DFLAB(M,NB) )
            INCG = 2
            CALL DIFC_LA( TX,VISLB(M,NB),VISGOX,DFLOB(M,NB) )
            CALL DIFC_LS( TX,XLSB(M,NB),VISLB(M,NB),DFLSB(M,NB) )
          ENDIF
!
!---      Nonisothermal simulation  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!  
!---        Aqueous thermal conductivity  ---
!
            CALL THK_W( TX,PX,RHOLWX,THKLWX )
            CALL THK_B( TX,XLSB(M,NB),THKLWX,THKLB(M,NB) )
!
!---        Gas thermal conductivity  ---
!
            CALL THK_A( TX,PGX,RHOGAX,THKGAX )
            CALL THK_O( TX,RHOGOX,THKGOX,VISGOOX )
            CALL THK_W( TX,PGX,RHOGWX,THKGWX )
            CALL THK_G( TX,THKGAX,THKGOX,THKGWX,
     &        XMGAB(M,NB),XMGOB(M,NB),XMGWB(M,NB),THKGB(M,NB) )
!
!---        Liquid-CO2 thermal conductivity  ---
!
            CALL THK_A( TX,PX,RHONB(M,NB),THKNB(M,NB) )
!
!---        Ice thermal conductivity  ---
!
            CALL THK_I( TX,THKIB(M,NB) )
!
!---        Hydrate thermal conductivity  ---
!
            CALL THK_H( TX,THKHB(M,NB) )
!
!---        Gas enthalpy and internal energy  ---
!
            ISRX = 2
            CALL ENTH_A( TX,PVAB(M,NB),HGAB(M,NB),UEGAB(M,NB),ISRX )
            CALL ENTH_O( TX,PVOB(M,NB),HGOB(M,NB),UGOX )
            CALL ENTH_W( TX,PVWB(M,NB),HLWX,HGWB(M,NB),ISRX )
            UGWX = HGWB(M,NB) - PVWB(M,NB)/RHOGWX
            UEGB(M,NB) = XGAB(M,NB)*UEGAB(M,NB) + 
     &        XGOB(M,NB)*UGOX + XGWB(M,NB)*UGWX
            HGB(M,NB) = XGAB(M,NB)*HGAB(M,NB) +  
     &        XGOB(M,NB)*HGOB(M,NB) + XGWB(M,NB)*HGWB(M,NB)
!
!---        Liquid water enthalpy  ---
!
            ISRX = 1
            CALL ENTH_W( TX,PX,HLWB(M,NB),HGWX,ISRX )
!
!---        Liquid-CO2 enthalpy  ---
!
            ISRX = 1
            CALL ENTH_A( TX,PX,HNAX,UNX,ISRX )
            FGWX = SGB(M,NB)*RHOGB(M,NB)*XGWB(M,NB) + SMALL
            FLWX = SLB(M,NB)*RHOLB(M,NB)*XLWB(M,NB) + SMALL
            HNWX = (HGWB(M,NB)*FGWX + HLWB(M,NB)*FLWX)/(FGWX + FLWX)
            HNB(M,NB) = XNAB(M,NB)*HNAX + XNWB(M,NB)*HNWX + 
     &        XNOB(M,NB)*HGOB(M,NB)
!
!---        Aqueous enthalpy  ---
!
            CALL ENTH_B( TX,XLSB(M,NB),HLWB(M,NB),HBX )
            FGAX = SGB(M,NB)*RHOGB(M,NB)*XGAB(M,NB) + SMALL
            FNAX = SNB(M,NB)*RHONB(M,NB)*XNAB(M,NB) + SMALL
            HGAX = (HGAB(M,NB)*FGAX + HNB(M,NB)*FNAX)/(FGAX + FNAX)
            CALL ENTH_L( TX,XLSB(M,NB),XLAB(M,NB),XLOB(M,NB),HBX,
     &        HGAX,HGOB(M,NB),HLB(M,NB) )
!
!---        Ice enthalpy  ---
!
            CALL ENTH_I( TX,HIB(M,NB) )
!
!---        Hydrate enthalpy  ---
!
            CALL ENTH_H( TX,XHAB(M,NB),XHOB(M,NB),XHWB(M,NB),
     &        HGAB(M,NB),HGOB(M,NB),HHB(M,NB),HLWB(M,NB) )
            
          ENDIF
!
!---      NaCl volumetric concentration  ---
!
          TMSB(M,NB) = YLSB(M,NB)*RHOLB(M,NB)*SLB(M,NB)*PORDB(M,NB)
  300   CONTINUE
  400 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCP37 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CAP37( CPGLX,CPGNX,CPHLX,CPILX,CPNLX,SHX,SIX,SLX,SNX,
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
!     H2O-CO2-CH4-NaCl-E-HYD1 Mode
!
!     Compute the gas-aqueous and hydrate aqueous capillary pressures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 27 September 2004
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
      REAL*8 GX(8)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CAP37'
!
!---  Residual aqueous saturation  ---
!
      SLRX = SCHR(4,IZN)
!
!---  van Genuchten saturation function
!
      IF( ISCHR(IZN).EQ.1 ) THEN
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
        ASTX = (((SLX+SNX)/(1.D+0-SHX-SIX))-SLRX)/(1.D+0-SLRX)
        ASLX = ((SLX/(1.D+0-SHX-SIX))-SLRX)/(1.D+0-SLRX)
        ASNX = ASTX - ASLX
        ASHX = SHX/(1.D+0-SLRX)
        ASIX = SIX/(1.D+0-SLRX)
        ESLHX = MAX( ASLX-ASHX,0.D+0 )
        ESLIX = MAX( ASLX-ASHX-ASIX,0.D+0 )
        IF( ASTX.GT.EPSL ) THEN
          HDGN = (((1.D+0/ASTX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))
     &      /SCHR(1,IZN)
        ELSE
          HDGN = 0.D+0
        ENDIF
        IF( ASLX.GT.EPSL ) THEN
          HDNL = (((1.D+0/ASLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))
     &      /SCHR(1,IZN)
        ELSE
          HDNL = 0.D+0
        ENDIF
        IF( ESLHX.GT.EPSL .AND. ASHX.GT.EPSL ) THEN
          HDHL = (((1.D+0/ESLHX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))
     &      /SCHR(1,IZN)
        ELSE
          HDHL = 0.D+0
        ENDIF
        IF( ESLIX.GT.EPSL .AND. ASIX.GT.EPSL ) THEN
          HDIL = (((1.D+0/ESLIX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))
     &      /SCHR(1,IZN)
        ELSE
          HDIL = 0.D+0
        ENDIF
        CPGNX = HDGN*RHORL*GRAV/BGN
        CPHLX = HDHL*RHORL*GRAV/BHL
        CPILX = HDIL*RHORL*GRAV/BIL
        CPNLX = HDNL*RHORL*GRAV/BNL
        CPGLX = CPGNX+CPNLX
!
!---  Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 ) THEN
        CL = MAX( SCHR(3,IZN),SMALL )
        ASTX = (((SLX+SNX)/(1.D+0-SHX-SIX))-SLRX)/(1.D+0-SLRX)
        ASLX = ((SLX/(1.D+0-SHX-SIX))-SLRX)/(1.D+0-SLRX)
        ASNX = ASTX - ASLX
        ASHX = SHX/(1.D+0-SLRX)
        ASIX = SIX/(1.D+0-SLRX)
        ESLHX = MAX( ASLX-ASHX,0.D+0 )
        ESLIX = MAX( ASLX-ASHX-ASIX,0.D+0 )
        IF( ASTX.GT.EPSL ) THEN
          HDGN = SCHR(2,IZN)*(1.D+0/ASTX)**(1.D+0/CL)
        ELSE
          HDGN = 0.D+0
        ENDIF
        IF( ASLX.GT.EPSL ) THEN
          HDNL = SCHR(2,IZN)*(1.D+0/ASLX)**(1.D+0/CL)
        ELSE
          HDNL = 0.D+0
        ENDIF
        IF( ESLHX.GT.EPSL .AND. ASHX.GT.EPSL ) THEN
          HDHL = SCHR(2,IZN)*(1.D+0/ESLHX)**(1.D+0/CL)
        ELSE
          HDHL = 0.D+0
        ENDIF
        IF( ESLIX.GT.EPSL .AND. ASIX.GT.EPSL ) THEN
          HDIL = SCHR(2,IZN)*(1.D+0/ESLIX)**(1.D+0/CL)
        ELSE
          HDIL = 0.D+0
        ENDIF
        HDGL = HDGN - HDNL
        CPGNX = HDGN*RHORL*GRAV/BGN
        CPHLX = HDHL*RHORL*GRAV/BHL
        CPILX = HDIL*RHORL*GRAV/BIL
        CPNLX = HDNL*RHORL*GRAV/BNL
        CPGLX = CPGNX+CPNLX
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Check the thermodynamic and hydrologic states declared through
!     user inputs.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FDVS
      USE FDVP
      USE FDVH
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
      REAL*8 GX(4,2),RX(2,2),RPX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHK37'
      EPSLX = 1.D-4
!
!---  Read CO2 and CH4 property files  ---
!
      IF( LPTA.LE.0 ) THEN
        INDX = 22
        CHMSG = 'Parameter LPTA = 0 for CO2-CH4 Operational Mode'
        CALL WRMSGS( INDX )
      ENDIF
      CALL RDPF_37
      WTMA = 44.010D+0
      TCRA = 304.1D+0
      TTPA = 216.592D+0
      PCRA = 73.8D+5
      PTPA = 5.1795D+5
      VCRA = 93.9D+0
      TBA = 194.68D+0
      WTMO = 16.043D+0
      RCO = RCU/WTMO
      TFPO = 90.7D+0
      TBO = 111.6D+0
      TCRO = 190.564D+0
      PCRO = 45.992D+5
      VCRO = 99.2D+0
      ZCRO = 0.288D+0
      PAFO = 0.011D+0
      DPMO = 0.D+0
      CPOC(1) = 1.925D+1
      CPOC(2) = 5.213D-2
      CPOC(3) = 1.197D-5
      CPOC(4) = -1.132D-8
      IVAPO = 1
      SATOC(1) = -6.00435D+0
      SATOC(2) = 1.1885D+0
      SATOC(3) = -0.83408D+0
      SATOC(4) = -1.22833D+0
!
!---  Number of components  ---
!
      NGC = 3
      IF( NGC.GT.LNGC ) THEN
        INDX = 5
        CHMSG = 'Number of Components > Parameter LNGC'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Component name  ---
!
      GCNM(1) = 'co2'
      GCNM(2) = 'ch4'
      GCNM(3) = 'h2o'
!
!---  Component molecular weight, kg/kmol  ---
!
      GCPP(1,1) = 44.010D+0
      GCPP(1,2) = 16.043D+0
      GCPP(1,3) = 18.015D+0
!
!---  Component normal boiling point, K  ---
!
      GCPP(2,1) = 0.D+0
      GCPP(2,2) = 111.6D+0
      GCPP(2,3) = 373.2D+0
!
!---  Component critical temperature, K  ---
!
      GCPP(3,1) = 304.1D+0
      GCPP(3,2) = 190.4D+0
      GCPP(3,3) = 647.3D+0
!
!---  Component critical pressure, Pa  ---
!
      GCPP(4,1) = 73.8D+5
      GCPP(4,2) = 46.0D+5
      GCPP(4,3) = 221.2D+5
!
!---  Component critical molar volume, cm^3/mole ---
!
      GCPP(5,1) = 93.9D+0
      GCPP(5,2) = 99.2D+0
      GCPP(5,3) = 57.1D+0
!
!---  Component critical compressibility  ---
!
      GCPP(6,1) = 0.274D+0
      GCPP(6,2) = 0.288D+0
      GCPP(6,3) = 0.235D+0
!
!---  Component Pitzer acentric factor  ---
!
      GCPP(7,1) = 0.239D+0
      GCPP(7,2) = 0.011D+0
      GCPP(7,3) = 0.344D+0
!
!---  Component dipole moment, debyes  ---
!
      GCPP(8,1) = 0.D+0
      GCPP(8,2) = 0.D+0
      GCPP(8,3) = 1.8D+0
!
!---  Peng-Robinson pure component parameters (a, b, fw)  ---
!
      DO 80 IGC = 1,NGC
        GCPP(13,IGC) = 0.45724D+0*(RCU**2)*(GCPP(3,IGC)**2)/
     &    GCPP(4,IGC)
        GCPP(14,IGC) = 0.07780D+0*RCU*GCPP(3,IGC)/GCPP(4,IGC)
        GCPP(15,IGC) = 0.37464D+0 + 1.54226D+0*GCPP(7,IGC) -
     &    0.26992D+0*(GCPP(7,IGC)**2)
   80 CONTINUE
!
!---  Check initial temperature, aqueous pressure, gas pressure,
!     and aqueous saturation  ---
!
      INDX = 0
      DO 100 N = 1,NFLD
        DO 90 M = 1,ISVC+2
          T(M,N) = T(2,N)
          PG(M,N) = PG(2,N)
          PL(M,N) = PL(2,N)
          PN(M,N) = PN(2,N)
          SG(M,N) = SG(2,N)
          SL(M,N) = SL(2,N)
          SN(M,N) = SN(2,N)
   90   CONTINUE
        N_DB = N
        IF( IXP(N).EQ.0 ) GOTO 100
        IF( T(2,N).GT.374.14D+0 .OR. T(2,N).LT.-20.D+0 ) THEN
          INDX = 16
          IMSG = N
          RLMSG = T(2,N)
          CHMSG = 'Out of Range Initial Temperature(C) @ Node'
          CALL WRMSGS( INDX )
        ENDIF
        IF( PL(2,N).GT.8.D+8-PATM ) THEN
          INDX = 16
          IMSG = N
          CHMSG = 'Out of Range Initial Aqueous Pressure(Pa) @ Node'
          RLMSG = PL(2,N)+PATM
          CALL WRMSGS( INDX )
        ENDIF
        IF( PG(2,N).GT.8.D+8-PATM ) THEN
          INDX = 16
          IMSG = N
          CHMSG = 'Out of Range Initial Gas Pressure(Pa) @ Node'
          RLMSG = PG(2,N)+PATM
          CALL WRMSGS( INDX )
        ENDIF
  100 CONTINUE
      IF( INDX.GT.0 ) STOP
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
          PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
        ENDIF
  200 CONTINUE
!
!---  Convert initial conditions through flash calculations  ---
!
      DO 3000 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 3000
        IZN = IZ(N)
        N_DB = N
!
!---    Initial condition options
!       IC_NCG(1,N) CO2 concentration option  
!       IC_NCG(2,N) CH4 concentration option
!       IC_NCG(1,N) or IC_NCG(2,N) > 10 gas pressure specified
!       IC_NCG(1,N) or IC_NCG(2,N) > 100 temperature specified
!       IC_OPT(2,N) hydrate saturation option
!       ICBRN(N) salt concentration option  ---
!
        IC_PGAX = (MOD(IC_NCG(1,N),100)-MOD(IC_NCG(1,N),10))/10
        IC_PGOX = (MOD(IC_NCG(2,N),100)-MOD(IC_NCG(2,N),10))/10
        IC_PGX = MAX( IC_PGAX,IC_PGOX )
        IC_TAX = (IC_NCG(1,N)-MOD(IC_NCG(1,N),100))/100
        IC_TOX = (IC_NCG(2,N)-MOD(IC_NCG(2,N),100))/100
        IC_TX = MAX( IC_TAX,IC_TOX )
        IC_AX = MOD(IC_NCG(1,N),10)
        IC_OX = MOD(IC_NCG(2,N),10)
        IC_SX = ICBRN(N)
        IC_SHX = IC_OPT(2,N)
!
!---    CO2 initial conditions specified with no CO2 solution  ---
!
        IF( IC_AX.NE.0 .AND. ISLC(45).NE.0 ) THEN
          INDX = 4
          CHMSG = 'CO2 Initial Conditions Declared for ' //
     &      'No CO2 Solution.'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Hydrate saturation declared  ---
!
        IF( IC_SHX.EQ.1 ) THEN
!
!---      Hydrate former CO2 fraction declared  ---
!
          IF( IC_AX.EQ.5 .OR. IC_AX.EQ.6 ) THEN
!
!---        CH4 concentration declared  ---
!
            IF( IC_OX.NE.0 ) THEN
              INDX = 4
              CHMSG = 'Specification Error: ' // 
     &          'Hydrate Saturation + ' //
     &          'Hydrate CO2 Fraction + ' //
     &          'CH4 Concentration Declared.'
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Gas pressure declared  ---
!
            IF( IC_PGX.NE.0 ) THEN
              INDX = 4
              CHMSG = 'Specification Error: ' // 
     &          'Hydrate Saturation + ' //
     &          'Hydrate CO2 Fraction + ' //
     &          'Gas Pressure Declared.'
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
!
!---      Hydrate former CH4 fraction declared  ---
!
          IF( IC_OX.EQ.5 .OR. IC_OX.EQ.6 ) THEN
!
!---        Liquid-CO2 pressure declared  ---
!
            IF( IC_AX.EQ.9 ) THEN
!
!---          Temperature declared  ---
!
              IF( IC_TX.EQ.1 ) THEN
                INDX = 4
                CHMSG = 'Specification Error: ' // 
     &            'Hydrate Saturation + ' //
     &            'Hydrate CH4 Fraction + ' //
     &            'Liquid-CO2 Pressure + ' //
     &            'Temperature Declared.'
                CALL WRMSGS( INDX )
              ENDIF
!
!---        CO2 concentration declared  ---
!
            ELSEIF( IC_AX.NE.0 ) THEN
              INDX = 4
              CHMSG = 'Specification Error: ' // 
     &          'Hydrate Saturation + ' //
     &          'Hydrate CH4 Fraction + ' //
     &          'CO2 Concentration Declared.'
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Gas pressure declared  ---
!
            IF( IC_PGX.NE.0 ) THEN
              INDX = 4
              CHMSG = 'Specification Error: ' // 
     &          'Hydrate Saturation + ' //
     &          'Hydrate CH4 Fraction + ' //
     &          'Gas Pressure Declared.'
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Liquid-CO2 pressure and temperatur declared  ---
!
            IF( IC_AX.EQ.9 .AND. IC_TX.NE.0 ) THEN
              INDX = 4
              CHMSG = 'Specification Error: ' // 
     &          'Hydrate Saturation + ' //
     &          'Hydrate CH4 Fraction + ' //
     &          'Liquid-CO2 Pressure + ' //
     &          'Temperature Declared.'
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
!
!---    No hydrate saturation declared  ---
!
        ELSE
!
!---      Hydrate former CO2 fraction declared  ---
!
          IF( IC_AX.EQ.5 .OR. IC_AX.EQ.6 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Hydrate Saturation Not + ' //
     &        'Hydrate CO2 Fraction Declared.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate former CH4 fraction declared  ---
!
          IF( IC_OX.EQ.5 .OR. IC_OX.EQ.6 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Hydrate Saturation Not + ' //
     &        'Hydrate CH4 Fraction Declared.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Gas pressure declared  ---
!
          IF( IC_PGX.EQ.1 ) THEN
!
!---        Gas pressure declared with both CO2 and CH4
!           concentrations declared  ---
!
            IF( IC_AX.NE.0 .AND. IC_OX.NE.0 ) THEN
              INDX = 4
              CHMSG = 'Specification Error: ' // 
     &          'Gas Pressured + ' //
     &          'Both CO2 and CH4 Concentrations Declared.'
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
        ENDIF
!
!---    Salt concentration undeclared  ---
!
        IF( IC_SX.EQ.0 ) IC_SX = 3
!
!---    Aqueous salt relative saturation or
!       aqueous salt mass fraction  ---
!
        IF( IC_SX.EQ.2 .OR. IC_SX.EQ.3 ) THEN
!
!---      Aqueous salt relative saturation  ---
!
          IF( IC_SX.EQ.2 ) THEN
            PHILSX = TMS(2,N)
            CALL SOL_LS( T(2,N),XLSMX )
            YLS(2,N) = PHILSX*XLSMX
            XLS(2,N) = MIN( YLS(2,N),XLSMX )
!
!---      Aqueous salt mass fraction  ---
!
          ELSEIF( IC_SX.EQ.3 ) THEN
            YLS(2,N) = TMS(2,N)
            CALL SOL_LS( T(2,N),XLSMX )
            XLS(2,N) = MIN( YLS(2,N),XLSMX )
          ENDIF
!
!---      Hydrate saturation declared   ---
!
          IF( IC_SHX.EQ.1 ) THEN
            IF( INHYD(IZN).EQ.1 ) THEN
              INDX = 21
              RLMSG = SH(2,N)
              CHMSG = 'Hydrate Saturation in No-Hydrate Rock: ' //
     &          ROCK(IZN) // 'SH = '
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Hydrate CH4 mass/mole fraction of formers  ---
!
            IF( IC_OX.EQ.5 .OR. IC_OX.EQ.6 ) THEN
!
!---          Hydrate CH4 mass fraction of formers  ---
!
              IF( IC_OX.EQ.5 ) THEN
                YHOX = PVO(2,N)
!
!---          Hydrate CH4 mole fraction of formers  ---
!
              ELSEIF( IC_OX.EQ.6 ) THEN
                YMHOX = PVO(2,N)
                YHOX = YMHOX*WTMO/(YMHOX*WTMO + (1.D+0-YMHOX)*WTMA)
              ENDIF
!
!---          Liquid-CO2 pressure declared  ---
!
              IF( IC_AX.EQ.9 .AND. IC_TX.EQ.0 ) THEN
!
!---            Initial condition flash calculation, assuming 
!               liquid-CO2 and hydrate equilibrium, solving
!               for temperature  ---
!
                PN(2,N) = PVA(2,N)
                CALL FLH_395( PG(2,N),PH(2,N),PL(2,N),PN(2,N),
     &            PVA(2,N),PVO(2,N),T(2,N),RHOH(2,N),
     &            RHOL(2,N),SH(2,N),SI(2,N),SL(2,N),SN(2,N),
     &            YHOX,YLS(2,N),YMGO(2,N),XLA(2,N),XLO(2,N),XHA(2,N),
     &            XHO(2,N),XHW(2,N),XMLA(2,N),XMLO(2,N),INDX,IZN,N )
                IF( INDX.NE.0 ) CALL WRMSGS( INDX )
              ELSE
!
!---            Initial condition flash calculation, assuming
!               hydrate equilibrium  ---
!
                CALL FLH_305( PG(2,N),PH(2,N),PL(2,N),PN(2,N),
     &            PVA(2,N),PVO(2,N),T(2,N),RHOH(2,N),
     &            RHOL(2,N),SH(2,N),SI(2,N),SL(2,N),
     &            YHOX,YLS(2,N),YMGO(2,N),XLA(2,N),XLO(2,N),XHA(2,N),
     &            XHO(2,N),XHW(2,N),XMLA(2,N),XMLO(2,N),INDX,IZN,N )
                IF( INDX.NE.0 ) CALL WRMSGS( INDX )
              ENDIF
!
!---        Hydrate CO2 mass/mole fraction of formers  ---
!
            ELSEIF( IC_AX.EQ.5 .OR. IC_AX.EQ.6 ) THEN
!
!---          Hydrate CO2 mass fraction of formers  ---
!
              IF( IC_AX.EQ.5 ) THEN
                YHAX = PVA(2,N)
!
!---          Hydrate CH4 mass fraction of formers  ---
!
              ELSEIF( IC_AX.EQ.6 ) THEN
                YMHAX = PVA(2,N)
                YHAX = YMHAX*WTMA/(YMHAX*WTMA + (1.D+0-YMHAX)*WTMO)
              ENDIF
!
!---          Initial condition flash calculation, assuming
!             hydrate equilibrium  ---
!
              CALL FLH_350( PG(2,N),PH(2,N),PL(2,N),PN(2,N),
     &          PVA(2,N),PVO(2,N),T(2,N),RHOH(2,N),RHOL(2,N),SH(2,N),
     &          SI(2,N),SL(2,N),YHAX,YLS(2,N),YMGO(2,N),
     &          XLA(2,N),XLO(2,N),XHA(2,N),XHO(2,N),
     &          XHW(2,N),XMLA(2,N),XMLO(2,N),INDX,IZN,N )
              IF( INDX.NE.0 ) CALL WRMSGS( INDX )
!
!---        Liquid-CO2 pressure declared  ---
!
            ELSEIF( IC_AX.EQ.9 ) THEN
!
!---          Initial condition flash calculation, assuming 
!             liquid-CO2 and hydrate equilibrium  ---
!
              PN(2,N) = PVA(2,N)
              CALL FLH_390( PG(2,N),PH(2,N),PL(2,N),PN(2,N),
     &          PVA(2,N),PVO(2,N),T(2,N),RHOH(2,N),
     &          RHOL(2,N),SH(2,N),SI(2,N),SL(2,N),SN(2,N),
     &          YHOX,YLS(2,N),YMGO(2,N),XLA(2,N),XLO(2,N),XHA(2,N),
     &          XHO(2,N),XHW(2,N),XMLA(2,N),XMLO(2,N),INDX,IZN,N )
              IF( INDX.NE.0 ) CALL WRMSGS( INDX )
            ENDIF
!
!---      Gas pressure or no gas pressure, no hydrate saturation;
!         for no gas pressure specification, gas pressure set
!         equal to aqueous or liquid CO2 pressure plus entry
!         pressure  ---
!
          ELSEIF( IC_PGX.EQ.0 .OR. IC_PGX.EQ.1 ) THEN
!
!---        Gas CH4 mass/mole fraction of formers  ---
!
            IF( IC_OX.EQ.7 .OR. IC_OX.EQ.8 ) THEN
              IF( IC_OX.EQ.7 ) THEN
                YGOX = PVO(2,N)
                YMGO(2,N) = (YGOX/WTMO)/(YGOX/WTMO + (1.D+0-YGOX)/WTMA)
              ELSEIF( IC_OX.EQ.8 ) THEN
                YMGO(2,N) = PVO(2,N)
              ENDIF
              CALL FLH_387( PG(2,N),PH(2,N),PL(2,N),PN(2,N),PVA(2,N),
     &          PVO(2,N),PVWX,T(2,N),RHOH(2,N),RHOL(2,N),SH(2,N),
     &          SI(2,N),SL(2,N),YLS(2,N),YMGO(2,N),XLA(2,N),XLO(2,N),
     &          XHA(2,N),XHO(2,N),XHW(2,N),XMLA(2,N),XMLO(2,N),
     &          INDX,IZN,N )
              IF( INDX.NE.0 ) CALL WRMSGS( INDX )
!
!---        Gas CO2 mass/mole fraction of formers  ---
!
            ELSEIF( IC_AX.EQ.7 .OR. IC_AX.EQ.8 ) THEN
              IF( IC_AX.EQ.7 ) THEN
                YGAX = PVA(2,N)
                YMGAX = (YGAX/WTMA)/(YGAX/WTMA + (1.D+0-YGAX)/WTMO)
              ELSEIF( IC_AX.EQ.8 ) THEN
                YMGAX = PVA(2,N)
              ENDIF
              CALL FLH_378( PG(2,N),PH(2,N),PL(2,N),PN(2,N),PVA(2,N),
     &          PVO(2,N),PVWX,T(2,N),RHOH(2,N),RHOL(2,N),SH(2,N),
     &          SI(2,N),SL(2,N),YLS(2,N),YMGAX,YMGO(2,N),XLA(2,N),
     &          XLO(2,N),XHA(2,N),XHO(2,N),XHW(2,N),XMLA(2,N),
     &          XMLO(2,N),INDX,IZN,N )
              IF( INDX.NE.0 ) CALL WRMSGS( INDX )
            ENDIF
          ENDIF
        ELSE
          INDX = 4
          CHMSG = 'Uncoded Initial Condition Option:  ' //
     &      'Dissolved-Salt Volumetric Concentration.'
          CALL WRMSGS( INDX )
        ENDIF
        SG(2,N) = MAX( 1.D+0-SH(2,N)-SI(2,N)-SL(2,N)-SN(2,N),0.D+0 )
        IF( SG(2,N).LT.1.D-10 ) SG(2,N)=0.D+0
        SL(1,N) = SL(2,N)
        SI(1,N) = SI(2,N)
        SH(1,N) = SH(2,N)
        SN(1,N) = SN(2,N)
        SG(1,N) = SG(2,N)
        T(1,N) = T(2,N)
        PL(1,N) = PL(2,N)
        PG(1,N) = PG(2,N)
        PN(1,N) = PN(2,N)
        XLS(1,N) = XLS(2,N)
        YLS(1,N) = YLS(2,N)
        PVA(1,N) = PVA(2,N)
        PVO(1,N) = PVO(2,N)
        YMGO(1,N) = YMGO(2,N)
#ifdef k1_hyd
        VISDG(1,N) = YMGO(2,N)
        VISDG(2,N) = YMGO(2,N)
#endif
#ifdef k2_hyd
        VISDG(1,N) = PVO(2,N)
        VISDG(2,N) = PVO(2,N)
        VISDG(3,N) = PVA(2,N)
        VISDG(4,N) = PVA(2,N)
#endif
!
!---    Aqueous or aqueous-gas, no hydrate, no liquid CO2  ---
!
        IF( SH(2,N).LE.EPSL .AND. SN(2,N).LE.EPSL ) THEN
          NPHAZ(2,N) = 1
!
!---    Aqueous or aqueous-gas, hydrate, no liquid CO2
!
        ELSEIF( SH(2,N).GT.EPSL .AND. SN(2,N).LE.EPSL ) THEN
          NPHAZ(2,N) = 2
!
!---    Aqueous or aqueous-gas, no hydrate, low liquid CO2
!
        ELSEIF( SH(2,N).LE.EPSL .AND. SN(2,N).LE.1.D-2 ) THEN
          NPHAZ(2,N) = 3
!
!---    Aqueous or aqueous-gas, hydrate, low liquid CO2
!
        ELSEIF( SH(2,N).GT.EPSL .AND. SN(2,N).LE.1.D-2 ) THEN
          NPHAZ(2,N) = 4
!
!---    Aqueous or aqueous-gas, no hydrate, liquid CO2
!
        ELSEIF( SH(2,N).LE.EPSL .AND. SN(2,N).GT.1.D-2 ) THEN
          NPHAZ(2,N) = 5
!
!---    Aqueous or aqueous-gas, hydrate, liquid CO2
!
        ELSEIF( SH(2,N).GT.EPSL .AND. SN(2,N).GT.1.D-2 ) THEN
          NPHAZ(2,N) = 6
        ENDIF
!
!---    Assign dissolved-salt mass fractions,
!       isobrine option  ---
!
        IF( ISLC(32).EQ.1 ) THEN
          DO 300 M = 3,ISVC+2
            XLS(M,N) = XLS(2,N)
            YLS(M,N) = YLS(2,N)
  300     CONTINUE
        ENDIF
!
!---    Assign temperature,
!       isothermal option  ---
!
        IF( ISLC(30).EQ.1 ) THEN
          DO 310 M = 3,ISVC+2
            T(M,N) = T(2,N)
  310     CONTINUE
        ENDIF
 3000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CISC37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Compute initial solute concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/CISC37'
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
!---  End of CISC37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFINA37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Reset temporarily inactive nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 21 November 2007
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
      SUB_LOG(ISUB_LOG) = '/DFINA37'
      IF( ICNV.EQ.3 ) THEN
        DO 100 N = 1,NFLD
          IXP(N) = ABS(IXP(N))
  100   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFINA37 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Compute primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      USE FDVS
      USE FDVP
      USE FDVN
      USE FDVH
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
      REAL*8 GX(3,2),AJX(6,6),BJX(6),RKLX(3)
      INTEGER IJX(6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM37'
      PETA = 1.25D+0
      EPSLX = 1.D-4
      CAX = 0.114543D+0
      CBX = 17.09028D+0
      CCX = 0.881712D+0
!
!---  Phase options, compute phase condition   ---
!
      DO 100 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 100
        IZN = IZ(N)
        N_DB = N
!
!---    Interfacial surface tensions  ---
!
        CALL SOL_LS( T(2,N),XLSMX )
        XLSX = MIN( YLS(2,N),XLSMX )
        SIG_GL = SCHR(5,IZN)
        SIG_NL = SCHR(6,IZN)
        SIG_GN = SCHR(7,IZN)
        SIG_HL = SCHR(9,IZN)
        SIG_IL = SCHR(9,IZN)
        CALL SFT_L( T(2,N),XLSX,SIG_GL )
        CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &    /(2.D+0*SIG_GL*SIG_NL)
        CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &    /(2.D+0*SIG_GL*SIG_GN)
!
!---    Interphase scaling factors  ---
!
        BGL = SCHR(8,IZN)/SIG_GL
        BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
        BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
        BHL = SCHR(8,IZN)/SCHR(9,IZN)
        BIL = SCHR(8,IZN)/SCHR(9,IZN)
        HGL = MAX( 0.D+0,BGL*(PG(2,N)-PL(2,N)) )
        HNL = MAX( 0.D+0,BNL*(PN(2,N)-PL(2,N)) )
        HGN = MAX( 0.D+0,BGN*(PG(2,N)-PN(2,N)) )
!
!---    Saturated CO2 vapor pressure  ---
!
        CALL SP_A( T(2,N),PSAX )
        TKX = T(2,N)+TABS
!
!---    Aqueous or aqueous-gas, no hydrate, no liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - CO2 vapor pressure
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          PSBX = 1.D+0
          SGX = SG(2,N)
          CALL FLH_PC1( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &      RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &      SI(2,N),SGX,SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &      XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      Transition to hydrate, low liquid-CO2 conditions (1 -> 4)  ---
!
          IF( T(2,N).LE.TEQHX
     &      .AND. PVA(2,N).GT.PSAX .AND. PVA(2,N).LT.PCRA
     &      .AND. TKX.GT.TTPA .AND. TKX.LT.TCRA ) THEN
!            SN(2,N) = 1.D-9
!            NPHAZ(2,N) = 4
#ifdef eq_hyd
            NPHAZ(2,N) = 6
            PSBX = 1.D+0
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &        PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
            PN(2,N) = PN(2,N) + 1.D-1
!
!---      Transition to low liquid-CO2 conditions (1 -> 3)  ---
!
          ELSEIF( PVA(2,N).GT.PSAX .AND. PVA(2,N).LT.PCRA
     &      .AND. TKX.GT.TTPA .AND. TKX.LT.TCRA ) THEN
!            SN(2,N) = 1.D-9
!            NPHAZ(2,N) = 3
            NPHAZ(2,N) = 5
            PVA(2,N) = PSAX
            PSBX = 1.D+0
            CALL FLH_PC5( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            PN(2,N) = PN(2,N) + 1.D-1
!
!---      Transition to hydrate conditions (1 -> 2)  ---
!
          ELSEIF( T(2,N).LE.TEQHX ) THEN
            NPHAZ(2,N) = 2
            PSBX = 1.D+0
#ifdef eq_hyd
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &        PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
!
!---      No transition (1 -> 1)  ---
!
          ELSE
!
!---        Gas appears  ---
!
            IF( SGX.GT.EPSL .AND. SG(2,N).LT.EPSL ) THEN
              SGX = 1.D-2
              SLX = 1.D+0-SH(2,N)-SI(2,N)-SN(2,N)-SGX
              CALL CAP37( CPGLX,CPGNX,CPHLX,CPILX,CPNLX,
     &          SH(2,N),SI(2,N),SLX,SN(2,N),IZN,N )            
              PL(2,N) = PG(2,N)-CPGLX
              PN(2,N) = PG(2,N)-CPGNX
            ENDIF
            NPHAZ(2,N) = 1
          ENDIF
!
!---    Aqueous or aqueous-gas, hydrate, no liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - gas CH4 mole fraction of formers
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          PSBX = 1.D+0
#ifdef eq_hyd
          SGX = SG(2,N)
          CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &      RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SGX,
     &      SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &      XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
          YMGOOX = VISDG(2,N)
          SGX = SG(2,N)
          CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &    PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &    RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SGX,
     &    SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &    XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
          VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
          PVOHOX = VISDG(2,N)
          PVAHOX = VISDG(4,N)
          SGX = SG(2,N)
          CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &      PVW(2,N),RHOH(2,N),RKG(2,N),
     &      RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SGX,
     &      SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &      XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
          VISDG(1,N) = PVOHX
          VISDG(3,N) = PVAHX
#endif
!
!---      Transition to no hydrate, 
!         low liquid-CO2 conditions (2 -> 3)  ---
!
          IF( SH(2,N).LE.EPSL
     &      .AND. PVA(2,N).GT.PSAX .AND. PVA(2,N).LT.PCRA
     &      .AND. TKX.GT.TTPA .AND. TKX.LT.TCRA ) THEN
!            SN(2,N) = 1.D-9
!            NPHAZ(2,N) = 3
            NPHAZ(2,N) = 5
            PVA(2,N) = PSAX
            PSBX = 1.D+0
            CALL FLH_PC5( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            PN(2,N) = PN(2,N) + 1.D-1
!
!---      Transition to hydrate, low liquid-CO2 conditions (2 -> 4)  ---
!
          ELSEIF( PVA(2,N).GT.PSAX .AND. PVA(2,N).LT.PCRA
     &      .AND. TKX.GT.TTPA .AND. TKX.LT.TCRA ) THEN
!            SN(2,N) = 1.D-9
!            NPHAZ(2,N) = 4
            NPHAZ(2,N) = 6
            SNX = 1.D-2
            SGLX = 1.D+0-SH(2,N)-SI(2,N)-SNX
            SLX = SGLX*SL(2,N)/(SGX+SL(2,N))
            CALL CAP37( CPGLX,CPGNX,CPHLX,CPILX,CPNLX,SH(2,N),SI(2,N),
     &        SLX,SNX,IZN,N )            
            PG(2,N) = MAX(PSAX+PVO(2,N)+PVW(2,N)-PATM,PL(2,N)+ENPR/BGL)
            PN(2,N) = PL(2,N) + CPNLX
            PSBX = 1.D+0
#ifdef eq_hyd
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &        PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
!
!---      Transition to no hydrate, 
!         no liquid-CO2 conditions (2 -> 1)  ---
!
          ELSEIF( SH(2,N).LE.EPSL ) THEN
            NPHAZ(2,N) = 1
            PSBX = 1.D+0
            CALL FLH_PC1( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      No transition from hydrate, 
!         no liquid-CO2 conditions (2 -> 2)  ---
!
          ELSE
!
!---        Gas appears  ---
!
            IF( SGX.GT.EPSL .AND. SG(2,N).LT.EPSL ) THEN
              SGX = 1.D-2
              SLX = 1.D+0-SH(2,N)-SI(2,N)-SN(2,N)-SGX
              CALL CAP37( CPGLX,CPGNX,CPHLX,CPILX,CPNLX,
     &          SH(2,N),SI(2,N),SLX,SN(2,N),IZN,N )            
              PL(2,N) = PG(2,N)-CPGLX
              PN(2,N) = PG(2,N)-CPGNX
            ENDIF
            NPHAZ(2,N) = 2
          ENDIF
!
!---    Aqueous or aqueous-gas, no hydrate, low liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 saturation
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
          PVA(2,N) = PSAX
          PSBX = 1.D+0
          CALL FLH_PC3( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &      RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &      SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &      XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      Transition to hydrate, liquid-CO2  conditions (3 -> 6)  ---
!
          IF( T(2,N).LE.TEQHX .AND. TKX.GT.TTPA
     &       .AND. TKX.LT.TCRA .AND. SN(2,N).GT.1.D-2 ) THEN
            NPHAZ(2,N) = 6
            PSBX = 1.D+0
#ifdef eq_hyd
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &        PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
            PN(2,N) = PN(2,N) + 1.D+1
!
!---      Transition to no hydrate, liquid-CO2 conditions (3 -> 5)  ---
!
          ELSEIF( TKX.GT.TTPA .AND. TKX.LT.TCRA 
     &      .AND. SN(2,N).GT.1.D-2 ) THEN
            NPHAZ(2,N) = 5
            PVA(2,N) = PSAX
            PSBX = 1.D+0
            CALL FLH_PC5( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            PN(2,N) = PN(2,N) + 1.D+1
!
!---      Transition to hydrate, no liquid-CO2 conditions (3 -> 2)  ---
!
          ELSEIF( T(2,N).LE.TEQHX .AND. SN(2,N).LT.EPSL ) THEN
            NPHAZ(2,N) = 2
            PSBX = 1.D+0
#ifdef eq_hyd
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &        PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
!
!---      Transition to no hydrate, no liquid-CO2 (3 -> 1)  ---
!
          ELSEIF( SN(2,N).LT.EPSL ) THEN
            NPHAZ(2,N) = 1
            PSBX = 1.D+0
            CALL FLH_PC1( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      Transition to hydrate, 
!         low liquid-CO2  conditions (3 -> 4)  ---
!
          ELSEIF( T(2,N).LE.TEQHX .AND. TKX.GT.TTPA 
     &       .AND. TKX.LT.TCRA ) THEN
            NPHAZ(2,N) = 4
            PSBX = 1.D+0
            CALL FLH_PC4( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      No transition from no hydrate, 
!         low liquid-CO2 conditions (3 -> 3)  ---
!
          ELSEIF( TKX.GT.TTPA .AND. TKX.LT.TCRA ) THEN
            NPHAZ(2,N) = 3
          ENDIF
!
!---    Aqueous or aqueous-gas, hydrate, low liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 saturation
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
          PSBX = 1.D+0
          CALL FLH_PC4( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &      RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &      SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &      XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      Transition to no hydrate, liquid-CO2 conditions (4 -> 5)  ---
!
          IF( TKX.GT.TTPA .AND. TKX.LT.TCRA 
     &      .AND. SN(2,N).GT.1.D-2 .AND. SH(2,N).LE.EPSL ) THEN
            NPHAZ(2,N) = 5
            PSBX = 1.D+0
            CALL FLH_PC5( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            PN(2,N) = PN(2,N) + 1.D+1
!
!---      Transition to hydrate, liquid-CO2  conditions (4 -> 6)  ---
!
          ELSEIF( TKX.GT.TTPA .AND. TKX.LT.TCRA 
     &      .AND. SN(2,N).GT.1.D-2 ) THEN
            NPHAZ(2,N) = 6
            PSBX = 1.D+0
#ifdef eq_hyd
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &        PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
            PN(2,N) = PN(2,N) + 1.D+1
!
!---      Transition to no hydrate, no liquid-CO2 (4 -> 1)  ---
!
          ELSEIF( SN(2,N).LT.EPSL .AND. SH(2,N).LT.EPSL ) THEN
            NPHAZ(2,N) = 1
            PSBX = 1.D+0
            CALL FLH_PC1( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      Transition to hydrate, no liquid-CO2 conditions (4 -> 2)  ---
!
          ELSEIF( SN(2,N).LT.EPSL ) THEN
            NPHAZ(2,N) = 2
            PSBX = 1.D+0
#ifdef eq_hyd
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),
     &        XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),
     &        XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &        PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),
     &        XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
!
!---      Transition to no hydrate, 
!         low liquid-CO2 conditions (4 -> 3)  ---
!
          ELSEIF( TKX.GT.TTPA .AND. TKX.LT.TCRA 
     &      .AND. SH(2,N).LT.EPSL ) THEN
            NPHAZ(2,N) = 3
            PSBX = 1.D+0
            CALL FLH_PC3( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      No transition from hydrate, 
!         low liquid-CO2 conditions (4 -> 4)  ---
!
          ELSEIF( TKX.GT.TTPA .AND. TKX.LT.TCRA ) THEN
            NPHAZ(2,N) = 4
          ENDIF
!
!---    Aqueous or aqueous-gas, no hydrate, liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 pressure
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
          PVA(2,N) = PSAX
          PSBX = 1.D+0
          SGX = SG(2,N)
          CALL FLH_PC5( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &      RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &      SI(2,N),SGX,SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &      XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      Transition to hydrate, 
!         low liquid-CO2  conditions (5 -> 4)  ---
!
!          IF( T(2,N).LE.TEQHX .AND. SN(2,N).LE.1.D-2 ) THEN
!            NPHAZ(2,N) = 4
          IF( T(2,N).LE.TEQHX .AND. SN(2,N).LE.EPSL ) THEN
            NPHAZ(2,N) = 2
            PSBX = 1.D+0
#ifdef eq_hyd
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),
     &        XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),
     &        XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &        PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),
     &        XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
!
!---      Transition from no hydrate, 
!         low liquid-CO2 conditions (5 -> 3)  ---
!
!          ELSEIF( SN(2,N).LE.1.D-2 ) THEN
!            NPHAZ(2,N) = 3
          ELSEIF( SN(2,N).LE.EPSL ) THEN
            NPHAZ(2,N) = 1
            PVA(2,N) = PSAX
            PSBX = 1.D+0
            CALL FLH_PC1( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      Transition to hydrate, liquid-CO2  conditions (5 -> 6)  ---
!
          ELSEIF( T(2,N).LE.TEQHX ) THEN
            NPHAZ(2,N) = 6
            PSBX = 1.D+0
#ifdef eq_hyd
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &        PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
!
!---      No transition to no hydrate, 
!         liquid-CO2 conditions (5 -> 5)  ---
!
          ELSE
!
!---        Gas appears  ---
!
            IF( SGX.GT.EPSL .AND. SG(2,N).LT.EPSL ) THEN
              SGX = 1.D-2
              SLX = 1.D+0-SH(2,N)-SI(2,N)-SN(2,N)-SGX
              CALL CAP37( CPGLX,CPGNX,CPHLX,CPILX,CPNLX,
     &          SH(2,N),SI(2,N),SLX,SN(2,N),IZN,N )            
              PN(2,N) = PG(2,N)-CPGNX
              PL(2,N) = PN(2,N)-CPNLX
            ENDIF
            NPHAZ(2,N) = 5
          ENDIF
!
!---    Aqueous or aqueous-gas, hydrate, liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 pressure
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
          PSBX = 1.D+0
#ifdef eq_hyd
          SGX = SG(2,N)
          SLX = SL(2,N)
          SNX = SN(2,N)
          CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &      RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SGX,
     &      SLX,SNX,T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &      XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
          SGX = SG(2,N)
          SLX = SL(2,N)
          SNX = SN(2,N)
          YMGOOX = VISDG(2,N)
          CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &      RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SGX,
     &      SLX,SNX,T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &      XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
          VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
          SGX = SG(2,N)
          SLX = SL(2,N)
          SNX = SN(2,N)
          PVOHOX = VISDG(2,N)
          PVAHOX = VISDG(4,N)
          CALL FLH_PC6( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &      PVW(2,N),RHOH(2,N),RKG(2,N),
     &      RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SGX,
     &      SLX,SNX,T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &      XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
          VISDG(1,N) = PVOHX
          VISDG(3,N) = PVAHX
#endif
!
!---      Limit changes in gas and liquid-CO2 saturation  ---
!
          DSGMX = MAX( MIN( 6.D-1*SG(2,N),5.D-2 ),1.D-3 )
          DSNMX = MAX( MIN( 6.D-1*SN(2,N),5.D-2 ),1.D-3 )
          IF( ABS(SGX-SG(2,N)).GT.DSGMX ) THEN
            SGX = SG(2,N) + SIGN(DSGMX,SGX-SG(2,N))
          ENDIF
          IF( ABS(SNX-SN(2,N)).GT.DSNMX ) THEN
            SNX = SN(2,N) + SIGN(DSNMX,SNX-SN(2,N))
          ENDIF
          SLX = 1.D+0-SH(2,N)-SI(2,N)-SGX-SNX
          CALL CAP37( CPGLX,CPGNX,CPHLX,CPILX,CPNLX,
     &      SH(2,N),SI(2,N),SLX,SNX,IZN,N )            
          PN(2,N) = PG(2,N)-CPGNX
          PL(2,N) = PN(2,N)-CPNLX
!!
!!---      Gas appears and liquid-CO2 disappears, 
!!         no transition from hydrate,
!!         liquid-CO2 conditions (6 -> 6)  ---
!!
!          IF( SGX.GT.EPSL .AND. SG(2,N).LT.EPSL .AND.
!     &      SNX.LT.EPSL .AND. SN(2,N).GT.1.D-3 ) THEN
!              SGLNX = 1.D+0-SH(2,N)-SI(2,N)
!              SGX = 1.D-2*SGLNX
!              SNX = 6.D-1*SN(2,N)
!              SLX = SGLNX-SGX-SNX
!              CALL CAP37( CPGLX,CPGNX,CPHLX,CPILX,CPNLX,
!     &          SH(2,N),SI(2,N),SLX,SNX,IZN,N )            
!              PN(2,N) = PG(2,N)-CPGNX
!              PL(2,N) = PN(2,N)-CPNLX
!              NPHAZ(2,N) = 6
!!
!!---      Gas appears, no transition from hydrate,
!!         liquid-CO2 conditions (6 -> 6)  ---
!!
!          ELSEIF( SGX.GT.EPSL .AND. SG(2,N).LT.EPSL ) THEN
!            SGLNX = 1.D+0-SH(2,N)-SI(2,N)
!            SGX = 1.D-2*SGLNX
!            SLNX = 1.D+0-SH(2,N)-SI(2,N)-SGX
!            SLX = SLNX*SLX/(SLX+SNX)
!            SNX = SLNX*SNX/(SLX+SNX)
!            CALL CAP37( CPGLX,CPGNX,CPHLX,CPILX,CPNLX,
!     &        SH(2,N),SI(2,N),SLX,SNX,IZN,N )            
!            PN(2,N) = PG(2,N)-CPGNX
!            PL(2,N) = PN(2,N)-CPNLX
!            NPHAZ(2,N) = 6
!!
!!---      Liquid-CO2 disappears, 
!!         no transition from hydrate,
!!         liquid-CO2 conditions (6 -> 6)  ---
!!
!          ELSEIF( SNX.LT.EPSL .AND. SN(2,N).GT.1.D-3 ) THEN
!            SNX = 6.D-1*SN(2,N)
!            SGLX = 1.D+0-SH(2,N)-SI(2,N)-SNX
!            SLX = SGLX*SL(2,N)/(SG(2,N)+SL(2,N))
!            SGX = SGLX*SG(2,N)/(SL(2,N)+SL(2,N))
!            CALL CAP37( CPGLX,CPGNX,CPHLX,CPILX,CPNLX,
!     &        SH(2,N),SI(2,N),SLX,SNX,IZN,N )            
!            PN(2,N) = PG(2,N)-CPGNX
!            PL(2,N) = PN(2,N)-CPNLX
!            NPHAZ(2,N) = 6
!
!---      Transition to no hydrate, 
!         low liquid-CO2 conditions (6 -> 3)  ---
!
!          IF( SH(2,N).LT.EPSL .AND. SN(2,N).LE.1.D-2 ) THEN
!            NPHAZ(2,N) = 3
          IF( SH(2,N).LT.EPSL .AND. SNX.LE.EPSL ) THEN
            NPHAZ(2,N) = 1
            PVA(2,N) = PSAX
            PSBX = 1.D+0
            CALL FLH_PC1( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      Transition to hydrate, low liquid-CO2 conditions (6 -> 4)  ---
!
!          ELSEIF( SN(2,N).LE.1.D-2 ) THEN
!            NPHAZ(2,N) = 4
          ELSEIF( SNX.LE.EPSL ) THEN
            NPHAZ(2,N) = 2
            PSBX = 1.D+0
#ifdef eq_hyd
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC2( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVAHX,PVAHOX,PVO(2,N),PVOHX,PVOHOX,
     &        PVW(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),
     &        SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),XHO(2,N),XHW(2,N),
     &        XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
!
!---      Transition to no hydrate, liquid-CO2 conditions (6 -> 5)  ---
!
          ELSEIF( SH(2,N).LE.EPSL ) THEN
            NPHAZ(2,N) = 5
            PSBX = 1.D+0
            CALL FLH_PC5( PEQHX,PG(2,N),PH(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PSBX,PVA(2,N),PVO(2,N),PVW(2,N),
     &        RHOH(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &        SI(2,N),SG(2,N),SL(2,N),SN(2,N),T(2,N),TEQHX,XHA(2,N),
     &        XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),YMGO(2,N),IZN,N )
!
!---      No transition from hydrate, 
!         liquid-CO2 conditions (6 -> 6)  ---
!
          ELSE
            NPHAZ(2,N) = 6
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
!---    Interfacial surface tensions  ---
!
        CALL SOL_LS( T(2,N),XLSMX )
        XLSX = MIN( YLS(2,N),XLSMX )
        SIG_GL = SCHR(5,IZN)
        SIG_NL = SCHR(6,IZN)
        SIG_GN = SCHR(7,IZN)
        SIG_HL = SCHR(9,IZN)
        CALL SFT_L( T(2,N),XLSX,SIG_GL )
        CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &    /(2.D+0*SIG_GL*SIG_NL)
        CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &    /(2.D+0*SIG_GL*SIG_GN)
!
!---    Interphase scaling factors  ---
!
        BGL = SCHR(8,IZN)/SIG_GL
        BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
        BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
        BHL = SCHR(8,IZN)/SCHR(9,IZN)
        BIL = SCHR(8,IZN)/SCHR(9,IZN)
!
!---    Energy (temperature), isothermal option  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          IF( SH(2,N).GT.EPSL .OR. SI(2,N).GT.EPSL ) THEN
            IF( SG(2,N).GT.EPSL ) THEN
              DNR(IEQT,N) = 1.D-7
            ELSE
              DNR(IEQT,N) = 1.D-6
            ENDIF
          ELSE
            DNR(IEQT,N) = 1.D-5
          ENDIF
        ENDIF
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
!---    Aqueous or aqueous-gas, no hydrate, no liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - CO2 vapor pressure
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
!
!---      Water mass (aqueous pressure)  ---
!
          DNR(IEQW,N) = -MAX( 1.D-2,1.D-6*(PG(2,N)-PL(2,N)) )
!
!---      CO2 mass (CO2 vapor pressure), iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) DNR(IEQA,N) = -DNR(IEQW,N)
!
!---      CH4 mass (CH4 vapor pressure)  ---
!
          DNR(IEQDO,N) = -DNR(IEQW,N)
!
!---    Aqueous or aqueous-gas, hydrate, no liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - gas CH4 mole fraction of formers
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
!
!---      Water mass (aqueous pressure)  ---
!
          DNR(IEQW,N) = -MAX( 1.D-2,1.D-6*(PG(2,N)-PL(2,N)) )
!
!---      CO2 mass (gas CH4 mole fraction of formers),
!         iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) DNR(IEQA,N) = SIGN(1.D-6,5.D-1-YMGO(2,N))
!
!---      CH4 mass (hydrate saturation)  ---
!
          SHX = MAX( 1.D-9,1.D-5*ABS(SH(2,N)) )
          DNR(IEQDO,N) = SIGN(SHX,5.D-1-SH(2,N))
!
!---    Aqueous or aqueous-gas, no hydrate, low liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid-CO2 saturation
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
!
!---      Water mass (aqueous pressure)  ---
!
          DNR(IEQW,N) = -1.D-2
!
!---      CO2 mass (liquid-CO2 saturation), iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            SNX = MAX( 1.D-9,1.D-5*ABS(SN(2,N)) )
            DNR(IEQA,N) = SIGN(SNX,5.D-1-SN(2,N))
          ENDIF
!
!---      CH4 mass (CH4 vapor pressure)  ---
!
          DNR(IEQDO,N) = 1.D-2
!
!---    Aqueous or aqueous-gas, hydrate, low liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 saturation
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
!
!---      Water mass (aqueous pressure)  ---
!
          DNR(IEQW,N) = -1.D-4
!
!---      CO2 mass (liquid-CO2 saturation), iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            SNX = MAX( 1.D-9,1.D-5*ABS(SN(2,N)) )
            DNR(IEQA,N) = SIGN(SNX,5.D-1-SN(2,N))
          ENDIF
!
!---      CH4 mass (hydrate saturation)  ---
!
          SHX = MAX( 1.D-9,1.D-5*ABS(SH(2,N)) )
          DNR(IEQDO,N) = SIGN(SHX,5.D-1-SH(2,N))
!
!---    Aqueous or aqueous-gas, no hydrate, liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 pressure
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
!
!---      Water mass (aqueous pressure)  ---
!
          DNR(IEQW,N) = -MAX( 1.D-2,1.D-6*ABS(PN(2,N)-PL(2,N)) )
!
!---      CO2 mass (liquid CO2 pressure), iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) DNR(IEQA,N) = -DNR(IEQW,N)
!
!---      CH4 mass (CH4 vapor pressure)  ---
!
          DNR(IEQDO,N) = 1.D-2
!
!---    Aqueous or aqueous-gas, hydrate, liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 pressure
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
!
!---      Water mass (aqueous pressure)  ---
!
          DNR(IEQW,N) = -MAX( 1.D-2,1.D-6*ABS(PN(2,N)-PL(2,N)) )
!
!---      CO2 mass (liquid CO2 pressure), iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) DNR(IEQA,N) = -DNR(IEQW,N)
!
!---      CH4 mass (hydrate saturation)  ---
!
          SHX = MAX( 1.D-9,1.D-5*ABS(SH(2,N)) )
          DNR(IEQDO,N) = SIGN(SHX,5.D-1-SH(2,N))
        ENDIF
!
!---    Increment the primary variables  ---
!
        DO 200 M = 3,ISVC+2
          T(M,N) = T(2,N)
          PG(M,N) = PG(2,N)
          PL(M,N) = PL(2,N)
          PN(M,N) = PN(2,N)
          SN(M,N) = SN(2,N)
          SL(M,N) = SL(2,N)
          SI(M,N) = SI(2,N)
          SH(M,N) = SH(2,N)
          PVA(M,N) = PVA(2,N)
          PVO(M,N) = PVO(2,N)
          YLS(M,N) = YLS(2,N)
          YMGO(M,N) = YMGO(2,N)
!
!---      Aqueous or aqueous-gas, no hydrate, no liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - CO2 vapor pressure
!         CH4 mass - CH4 vapor pressure
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          IF( NPHAZ(2,N).EQ.1 ) THEN
!
!---        Isothermal option  ---
!
            IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
!
!---        Iso-CO2 option  ---
!
            ELSEIF( M.EQ.IEQA+2 .AND. ISLC(45).EQ.0 ) THEN
              PVA(M,N) = PVA(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQDO+2 ) THEN
              PVO(M,N) = PVO(M,N) + DNR(IEQDO,N)
!
!---        Isobrine option  ---
!
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Aqueous or aqueous-gas, hydrate, no liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - CH4 hydrate mole fraction of formers
!         CH4 mass - hydrate saturation
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
!
!---        Isothermal option  ---
!
            IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
!
!---        Iso-CO2 option  ---
!
            ELSEIF( M.EQ.IEQA+2 .AND. ISLC(45).EQ.0 ) THEN
              YMGO(M,N) = YMGO(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQDO+2 ) THEN
              SH(M,N) = SH(M,N) + DNR(IEQDO,N)
!
!---        Isobrine option  ---
!
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Aqueous or aqueous-gas, no hydrate, low liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - liquid CO2 saturation
!         CH4 mass - CH4 vapor pressure
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
!
!---        Isothermal option  ---
!
            IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
!
!---        Iso-CO2 option  ---
!
            ELSEIF( M.EQ.IEQA+2 .AND. ISLC(45).EQ.0 ) THEN
              SN(M,N) = SN(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQDO+2 ) THEN
              PVO(M,N) = PVO(M,N) + DNR(IEQDO,N)
!
!---        Isobrine option  ---
!
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Aqueous or aqueous-gas, hydrate, low liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - liquid CO2 saturation
!         CH4 mass - hydrate saturation
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
!
!---        Isothermal option  ---
!
            IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
!
!---        Iso-CO2 option  ---
!
            ELSEIF( M.EQ.IEQA+2 .AND. ISLC(45).EQ.0 ) THEN
              SN(M,N) = SN(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQDO+2 ) THEN
              SH(M,N) = SH(M,N) + DNR(IEQDO,N)
!
!---        Isobrine option  ---
!
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Aqueous or aqueous-gas, no hydrate, liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - liquid CO2 pressure
!         CH4 mass - CH4 vapor pressure
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
!
!---        Isothermal option  ---
!
            IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
!
!---        Iso-CO2 option  ---
!
            ELSEIF( M.EQ.IEQA+2 .AND. ISLC(45).EQ.0 ) THEN
              PN(M,N) = PN(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQDO+2 ) THEN
              PVO(M,N) = PVO(M,N) + DNR(IEQDO,N)
!
!---        Isobrine option  ---
!
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Aqueous or aqueous-gas, hydrate, liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - liquid CO2 pressure
!         CH4 mass - hydrate saturation
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
!
!---        Isothermal option  ---
!
            IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
!
!---        Iso-CO2 option  ---
!
            ELSEIF( M.EQ.IEQA+2 .AND. ISLC(45).EQ.0 ) THEN
              PN(M,N) = PN(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQDO+2 ) THEN
              SH(M,N) = SH(M,N) + DNR(IEQDO,N)
!
!---        Isobrine option  ---
!
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
          ENDIF
  200   CONTINUE
  300 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM37 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Compute total aqueous+hydrate saturation from gas-aqueous head.
!     Compute aqueous saturation from hydrate-aqueous head.
!
!     INDX = 0 : Hydrate saturation given, NAPL saturation computed.
!     INDX = 1 : Hydrate saturation given, NAPL saturation given.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 September 2004.
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
      REAL*8 RKLX(3),GX(8)
!
!----------------------Data Statements---------------------------------!
!
      SAVE TOL
      DATA TOL /1.D-12/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/KSP37'
!
!---  Residual aqueous saturation  ---
!
      SLRX = SCHR(4,IZN)
!
!---  van Genuchten saturation function w/o or w/ the
!     hydrate capillary pressure barrier  ---
!
      IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11 ) THEN
        IF( SHX.LT.TOL ) SHX = 0.D+0
        HCPBFX = 1.D+0
        IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
          HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &      SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
        ENDIF
        HNL = MAX( ZERO,BNL*(PNX-PLX)/RHORL/GRAV )/HCPBFX
        HGN = MAX( ZERO,BGN*(PGX-PNX)/RHORL/GRAV )/HCPBFX
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
        IF( INDX.NE.1 ) THEN
          ASNX = ASTX-ASLX
          IF( ASNX.LT.TOL ) ASNX = 0.D+0
        ENDIF
        HIL = MAX( 0.D+0,BIL*(PIX-PLX)/RHORL/GRAV )
!
!---    Potential ice saturation  ---
!
        IF( HIL.GT.EPSL ) THEN
!!
!!---     Newton-Raphson loop to find ice saturation  ---
!!
!  100     CONTINUE
!          DO 110 M = 1,2
!            SIY = SIX
!            IF( M.EQ.2 ) SIY = SIX + 1.D-6
!            SLX = (ASLX*(1.D+0-SLRX) + SLRX)*(1.D+0-SHX-SIY)
!            IF( SLX.LT.TOL ) SLX = 0.D+0
!            IF( INDX.NE.1 ) THEN
!              SNX = (ASNX*(1.D+0-SLRX))*(1.D+0-SHX-SIY)
!              IF( SNX.LT.TOL ) SNX = 0.D+0
!            ENDIF
!            SGX = MAX( 1.D+0-SLX-SNX-SHX-SIY,0.D+0 )
!            IF( SGX.LT.TOL ) SGX = 0.D+0
!            ASLHX = MAX( 1.D+0-SNX-SGX,0.D+0 )
!            HHL = (((1.D+0/ASLHX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))
!     &        /SCHR(1,IZN)
!            IF( HIL.GT.HHL ) THEN
!              ASLIX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HIL)**CN))**CM
!              ASLIX = MAX( ASLIX,1.D-9 )
!              SIZ = MAX( ASLHX-ASLIX,0.D+0 )
!            ELSE
!              SIZ = 0.D+0
!            ENDIF
!            GX(M) = SIY - SIZ
!  110     CONTINUE
!          FX = GX(1)
!          DFX = (GX(2)-GX(1))*1.D+6
!          DSIX = -FX/DFX
!          SIX = SIX + DSIX
!          IF( ABS(DSIX).GT.1.D-9 ) GOTO 100
!
!---     Newton-Raphson loop to find ice saturation  ---
!
  100     CONTINUE
          DO 110 M = 1,2
            SIY = SIX
            IF( M.EQ.2 ) SIY = SIX + 1.D-6
            SLX = (ASLX*(1.D+0-SLRX) + SLRX)*(1.D+0-SHX-SIY)
            IF( SLX.LT.TOL ) SLX = 0.D+0
            IF( INDX.NE.1 ) THEN
              SNX = (ASNX*(1.D+0-SLRX))*(1.D+0-SHX-SIY)
              IF( SNX.LT.TOL ) SNX = 0.D+0
            ENDIF
            SGX = MAX( 1.D+0-SLX-SNX-SHX-SIY,0.D+0 )
            IF( SGX.LT.TOL ) SGX = 0.D+0
            ASLHX = MAX( 1.D+0-SNX-SGX,0.D+0 )
            IF( ASLHX.LT.EPSL ) THEN
              HHL = HDOD
            ELSE
            HHL = (((1.D+0/ASLHX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))
     &        /SCHR(1,IZN)
            ENDIF
            IF( HIL.GT.HHL ) THEN
              ASLIX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HIL)**CN))**CM
              ASLIX = MAX( ASLIX,1.D-9 )
              SIZ = MAX( ASLHX-ASLIX,0.D+0 )
            ELSE
              SIZ = 0.D+0
            ENDIF
            GX(M) = SIY - SIZ
  110     CONTINUE
          FX = GX(1)
          DFX = (GX(2)-GX(1))*1.D+6
          DSIX = -FX/DFX
          SIX = SIX + DSIX
          IF( ABS(DSIX).GT.1.D-9 ) GOTO 100
          SLX = (ASLX*(1.D+0-SLRX) + SLRX)*(1.D+0-SHX-SIY)
          IF( SLX.LT.TOL ) SLX = 0.D+0
          IF( INDX.NE.1 ) THEN
            SNX = (ASNX*(1.D+0-SLRX))*(1.D+0-SHX-SIX)
            IF( SNX.LT.TOL ) SNX = 0.D+0
          ENDIF
          SGX = MAX( 1.D+0-SLX-SNX-SHX-SIX,0.D+0 )
          IF( SGX.LT.TOL ) SGX = 0.D+0
!
!---    Zero ice saturation  ---
!
        ELSE
          SIX = 0.D+0
          SLX = (ASLX*(1.D+0-SLRX) + SLRX)*(1.D+0-SHX-SIX)
          IF( SLX.LT.TOL ) SLX = 0.D+0
          IF( INDX.NE.1 ) THEN
            SNX = (ASNX*(1.D+0-SLRX))*(1.D+0-SHX-SIX)
            IF( SNX.LT.TOL ) SNX = 0.D+0
          ENDIF
          SGX = MAX( 1.D+0-SLX-SNX-SHX-SIX,0.D+0 )
          IF( SGX.LT.TOL ) SGX = 0.D+0
        ENDIF
!!
!!---    Converged solution on ice saturation  ---
!!
!        SLX = (ASLX*(1.D+0-SLRX) + SLRX)*(1.D+0-SHX-SIX)
!        IF( SLX.LT.TOL ) SLX = 0.D+0
!        IF( INDX.NE.1 ) THEN
!          SNX = (ASNX*(1.D+0-SLRX))*(1.D+0-SHX-SIX)
!          IF( SNX.LT.TOL ) SNX = 0.D+0
!        ENDIF
!        SGX = MAX( 1.D+0-SLX-SNX-SHX-SIX,0.D+0 )
!        IF( SGX.LT.TOL ) SGX = 0.D+0
!!
!!---    Effective saturations for relative permeability
!!       calculations  ---
!!
!        ASLX = MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 )
!        ASNX = SNX/(1.D+0-SLRX)
!        ASLNX = ASLX + ASNX
!        ASGX = SGX/(1.D+0-SLRX)
!        ASTX = 1.D+0-ASGX
!!!
!!!---    Converged solution on ice saturation  ---
!!!
!!        SLX = ASLX*(1.D+0-SLRX) + SLRX - SHX - SIX
!!        IF( SLX.LT.TOL ) SLX = 0.D+0
!!        IF( INDX.NE.1 ) THEN
!!          SNX = (ASNX*(1.D+0-SLRX))
!!          IF( SNX.LT.TOL ) SNX = 0.D+0
!!        ENDIF
!!        SGX = MAX( 1.D+0-SLX-SNX-SHX-SIX,0.D+0 )
!!        IF( SGX.LT.TOL ) SGX = 0.D+0
!
!---    Effective saturations for relative permeability
!       calculations  ---
!
        ASLX = MAX( (SLX-SHX-SIX-SLRX)/(1.D+0-SLRX),0.D+0 )
        ASNX = SNX/(1.D+0-SLRX)
        ASLNX = ASLX + ASNX
        ASGX = SGX/(1.D+0-SLRX)
        ASTX = 1.D+0-ASGX
!
!---    Constant aqueous relative permeability  ---
!
        IF( IRPL(IZN).EQ.0 ) THEN
          RKLX(1) = RPLC(2,IZN)
!
!---    Mualem aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.1 ) THEN
          CM = RPLC(2,IZN)
          RKLX(1) = SQRT(ASLX)*(1.D+0 -
     &      (1.D+0-ASLX**(1.D+0/CM))**CM)**2
!
!---    Burdine aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.2 ) THEN
          CM = RPLC(2,IZN)
          RKLX(1) = (ASLX**2)*(1.D+0 -
     &      (1.D+0-ASLX**(1.D+0/CM))**CM)
!
!---    Free Corey relative permeability function,
!       modified from STOMP-CO2 for the Moridis
!       code comparison problem  ---
!
        ELSEIF( IRPL(IZN).EQ.7 ) THEN
          SLRX = RPLC(3,IZN)
          SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          RKLX(1) = RPLC(1,IZN)*(SLPX**(RPLC(2,IZN)))
!
!---    Modified Mualem aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.22 ) THEN
          RKLX(1) = (ASLX**RPLC(1,IZN))*((1.D+0-(1.D+0-ASLX**
     &      (1.D+0/RPLC(2,IZN)))**RPLC(2,IZN))**2)
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
!---    Constant NAPL relative permeability  ---
!
        IF( IRPN(IZN).EQ.0 ) THEN
          RKNX = RPNC(1,IZN)
!
!---    Mualem NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.1 ) THEN
          CM = RPNC(1,IZN)
          RKNX = SQRT(ASNX)*((1.D+0 - (ASLNX-ASNX)**(1.D+0/CM))**CM -
     &      (1.D+0 - (ASLNX)**(1.D+0/CM))**CM)**2
!
!---    Burdine NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.2 ) THEN
          CM = RPNC(1,IZN)
          RKNX = (ASNX**2)*((1.D+0 - (ASLNX-ASNX)**(1.D+0/CM))**CM -
     &      (1.D+0 - (ASLNX)**(1.D+0/CM))**CM)
!
!---    Free-Corey relative permeability function,
!       modified from STOMP-CO2 for the Moridis
!       code comparison problem  ---
!
        ELSEIF( IRPN(IZN).EQ.7 ) THEN
          SLRX = RPNC(3,IZN)
          SNRX = RPNC(4,IZN)
          SNPX = MIN( MAX( (SNX-SNRX)/(1.D+0-SLRX),0.D+0 ),
     &      1.D+0 )
          RKNX = RPNC(1,IZN)*(SNPX**(RPNC(2,IZN)))
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
!---    Constant gas relative permeability  ---
!
        IF( IRPG(IZN).EQ.0 ) THEN
          RKGX = RPGC(3,IZN)
!
!---    Mualem gas relative permeability  ---
!
        ELSEIF( IRPG(IZN).EQ.1 ) THEN
          CM = RPGC(3,IZN)
          RKGX = SQRT(ASGX)*((1.D+0-ASTX**(1.D+0/CM))**CM)**2
!
!---    Burdine gas relative permeability  ---
!
        ELSEIF( IRPG(IZN).EQ.2 ) THEN
          CM = RPGC(3,IZN)
          RKGX = (ASGX**2)*((1.D+0-ASTX**(1.D+0/CM))**CM)
!
!---    Free-Corey relative permeability function,
!       modified from STOMP-CO2 for the Moridis
!       code comparison problem  ---
!
        ELSEIF( IRPG(IZN).EQ.7 ) THEN
          SLRX = RPGC(3,IZN)
          SGRX = RPGC(4,IZN)
          SGPX = MIN( MAX( (SGX-SGRX)/(1.D+0-SLRX),0.D+0 ),
     &      1.D+0 )
          RKGX = RPGC(1,IZN)*(SGPX**(RPGC(2,IZN)))
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
!---  Brooks and Corey saturation function w/o or w/ the
!     hydrate capillary pressure barrier  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12 ) THEN
        IF( SHX.LT.TOL ) SHX = 0.D+0
        HCPBFX = 1.D+0
        IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
          HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &      SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
        ENDIF
        HNL = MAX( ZERO,BNL*(PNX-PLX)/RHORL/GRAV )/HCPBFX
        HGN = MAX( ZERO,BGN*(PGX-PNX)/RHORL/GRAV )/HCPBFX
        CL = MAX( SCHR(3,IZN),SMALL )
        IF( HGN-SCHR(2,IZN).LE.1.D-12 ) THEN
          ASTX = 1.D+0
        ELSE
          ASTX = (SCHR(2,IZN)/HGN)**CL
        ENDIF
        IF( HNL-SCHR(2,IZN).LE.1.D-12 ) THEN
          ASLX = 1.D+0
        ELSE
          ASLX = (SCHR(2,IZN)/HNL)**CL
        ENDIF
        IF( INDX.NE.1 ) THEN
          ASNX = ASTX-ASLX
          IF( ASNX.LT.TOL ) ASNX = 0.D+0
        ENDIF
        HIL = MAX( 0.D+0,BIL*(PIX-PLX)/RHORL/GRAV )
!
!---    Potential ice saturation  ---
!
        IF( HIL.GT.EPSL ) THEN
!!
!!---     Newton-Raphson loop to find ice saturation  ---
!!
!  200     CONTINUE
!          DO 210 M = 1,2
!            SIY = SIX
!            IF( M.EQ.2 ) SIY = SIX + 1.D-6
!            SLX = (ASLX*(1.D+0-SLRX) + SLRX)*(1.D+0-SHX-SIY)
!            IF( SLX.LT.TOL ) SLX = 0.D+0
!            IF( INDX.NE.1 ) THEN
!              SNX = (ASNX*(1.D+0-SLRX))*(1.D+0-SHX-SIY)
!              IF( SNX.LT.TOL ) SNX = 0.D+0
!            ENDIF
!            SGX = MAX( 1.D+0-SLX-SNX-SHX-SIY,0.D+0 )
!            IF( SGX.LT.TOL ) SGX = 0.D+0
!            ASLHX = MAX( 1.D+0-SNX-SGX,0.D+0 )
!            HHL = SCHR(2,IZN)/(ASLHX**(1.D+0/CL))
!            IF( HIL.GT.HHL ) THEN
!              IF( HIL-SCHR(2,IZN).LE.1.D-12 ) THEN
!                ASLIX = 1.D+0
!              ELSE
!                ASLIX = (SCHR(2,IZN)/HIL)**CL
!              ENDIF
!              ASLIX = MAX( ASLIX,1.D-9 )
!              SIZ = MAX( ASLHX-ASLIX,0.D+0 )
!            ELSE
!              SIZ = 0.D+0
!            ENDIF
!            GX(M) = SIY - SIZ
!  210     CONTINUE
!          FX = GX(1)
!          DFX = (GX(2)-GX(1))*1.D+6
!          DSIX = -FX/DFX
!          SIX = SIX + DSIX
!          IF( ABS(DSIX).GT.1.D-9 ) GOTO 200
!
!---     Newton-Raphson loop to find ice saturation  ---
!
  200     CONTINUE
          DO 210 M = 1,2
            SIY = SIX
            IF( M.EQ.2 ) SIY = SIX + 1.D-6
            SLX = (ASLX*(1.D+0-SLRX) + SLRX)*(1.D+0-SHX-SIY)
            IF( SLX.LT.TOL ) SLX = 0.D+0
            IF( INDX.NE.1 ) THEN
              SNX = (ASNX*(1.D+0-SLRX))*(1.D+0-SHX-SIY)
              IF( SNX.LT.TOL ) SNX = 0.D+0
            ENDIF
            SGX = MAX( 1.D+0-SLX-SNX-SHX-SIY,0.D+0 )
            IF( SGX.LT.TOL ) SGX = 0.D+0
            ASLHX = MAX( 1.D+0-SNX-SGX,0.D+0 )
            IF( ASLHX.LT.EPSL ) THEN
              HHL = HDOD
            ELSE
            HHL = SCHR(2,IZN)/(ASLHX**(1.D+0/CL))
            ENDIF
            IF( HIL.GT.HHL ) THEN
              IF( HIL-SCHR(2,IZN).LE.1.D-12 ) THEN
                ASLIX = 1.D+0
              ELSE
                ASLIX = (SCHR(2,IZN)/HIL)**CL
              ENDIF
              ASLIX = MAX( ASLIX,1.D-9 )
              SIZ = MAX( ASLHX-ASLIX,0.D+0 )
            ELSE
              SIZ = 0.D+0
            ENDIF
            GX(M) = SIY - SIZ
  210     CONTINUE
          FX = GX(1)
          DFX = (GX(2)-GX(1))*1.D+6
          DSIX = -FX/DFX
          SIX = SIX + DSIX
          IF( ABS(DSIX).GT.1.D-9 ) GOTO 200
          SIX = 0.D+0
          SLX = (ASLX*(1.D+0-SLRX) + SLRX)*(1.D+0-SHX-SIY)
          IF( SLX.LT.TOL ) SLX = 0.D+0
          IF( INDX.NE.1 ) THEN
            SNX = (ASNX*(1.D+0-SLRX))*(1.D+0-SHX-SIX)
            IF( SNX.LT.TOL ) SNX = 0.D+0
          ENDIF
          SGX = MAX( 1.D+0-SLX-SNX-SHX-SIX,0.D+0 )
          IF( SGX.LT.TOL ) SGX = 0.D+0
!
!---    Zero ice saturation  ---
!
        ELSE
          SIX = 0.D+0
          SLX = (ASLX*(1.D+0-SLRX) + SLRX)*(1.D+0-SHX-SIX)
          IF( SLX.LT.TOL ) SLX = 0.D+0
          IF( INDX.NE.1 ) THEN
            SNX = (ASNX*(1.D+0-SLRX))*(1.D+0-SHX-SIX)
            IF( SNX.LT.TOL ) SNX = 0.D+0
          ENDIF
          SGX = MAX( 1.D+0-SLX-SNX-SHX-SIX,0.D+0 )
          IF( SGX.LT.TOL ) SGX = 0.D+0
        ENDIF
!!
!!---    Converged solution on ice saturation  ---
!!
!        SLX = (ASLX*(1.D+0-SLRX) + SLRX)*(1.D+0-SHX-SIX)
!        IF( SLX.LT.TOL ) SLX = 0.D+0
!        IF( INDX.NE.1 ) THEN
!          SNX = (ASNX*(1.D+0-SLRX))*(1.D+0-SHX-SIX)
!          IF( SNX.LT.TOL ) SNX = 0.D+0
!        ENDIF
!        SGX = MAX( 1.D+0-SLX-SNX-SHX-SIX,0.D+0 )
!        IF( SGX.LT.TOL ) SGX = 0.D+0
!!
!!---    Effective saturations for relative permeability
!!       calculations  ---
!!
!        ASLX = MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 )
!        ASNX = SNX/(1.D+0-SLRX)
!        ASLNX = ASLX + ASNX
!        ASGX = SGX/(1.D+0-SLRX)
!        ASTX = 1.D+0-ASGX
!!
!!---    Converged solution on ice saturation  ---
!!
!!        SLX = ASLX*(1.D+0-SLRX) + SLRX - SHX - SIX
!!        IF( SLX.LT.TOL ) SLX = 0.D+0
!!        IF( INDX.NE.1 ) THEN
!!          SNX = (ASNX*(1.D+0-SLRX))
!!          IF( SNX.LT.TOL ) SNX = 0.D+0
!!        ENDIF
!!        SGX = MAX( 1.D+0-SLX-SNX-SHX-SIX,0.D+0 )
!!        IF( SGX.LT.TOL ) SGX = 0.D+0
!
!---    Effective saturations for relative permeability
!       calculations  ---
!
        ASLX = MAX( (SLX-SHX-SIX-SLRX)/(1.D+0-SLRX),0.D+0 )
        ASNX = SNX/(1.D+0-SLRX)
        ASLNX = ASLX + ASNX
        ASGX = SGX/(1.D+0-SLRX)
        ASTX = 1.D+0-ASGX
!
!---    Constant aqueous relative permeability  ---
!
        IF( IRPL(IZN).EQ.0 ) THEN
          RKLX(1) = RPLC(2,IZN)
!
!---    Mualem aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.1 ) THEN
          CL = RPLC(2,IZN)
          RKLX(1) = ASLX**(2.5D+0 + 2.D+0/CL)
!
!---    Burdine aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.2 ) THEN
          CL = RPLC(2,IZN)
          RKLX(1) = ASLX**(3.0D+0 + 2.D+0/CL)
!
!---    Free Corey relative permeability function,
!       modified from STOMP-CO2 for the Moridis
!       code comparison problem  ---
!
        ELSEIF( IRPL(IZN).EQ.7 ) THEN
          SLRX = RPLC(3,IZN)
          SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          RKLX(1) = RPLC(1,IZN)*(SLPX**(RPLC(2,IZN)))
!
!---    Modified Mualem aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.22 ) THEN
          RKLX(1) = (ASLX**RPLC(1,IZN))*
     &      ASLX**(2.D+0 + 2.0D+0/RPLC(2,IZN))
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
!---    Constant NAPL relative permeability  ---
!
        IF( IRPN(IZN).EQ.0 ) THEN
          RKNX = RPNC(1,IZN)
!
!---    Mualem NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.1 ) THEN
          CL = RPNC(1,IZN)
          RKNX = SQRT(ASNX)*( ASTX**(2.0D+0 + 2.D+0/CL) -
     &        (ASTX-ASNX)**(2.0D+0 + 2.D+0/CL) )
!
!---    Burdine NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.2 ) THEN
          CL = RPNC(1,IZN)
          RKNX = (ASNX**2)*( ASTX**(1.0D+0 + 2.D+0/CL) -
     &        (ASTX-ASNX)**(1.0D+0 + 2.D+0/CL) )
!
!---    Free-Corey relative permeability function,
!       modified from STOMP-CO2 for the Moridis
!       code comparison problem  ---
!
        ELSEIF( IRPN(IZN).EQ.7 ) THEN
          SLRX = RPNC(3,IZN)
          SNRX = RPNC(4,IZN)
          SNPX = MIN( MAX( (SNX-SNRX)/(1.D+0-SLRX),0.D+0 ),
     &      1.D+0 )
          RKNX = RPNC(1,IZN)*(SNPX**(RPNC(2,IZN)))
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
!---    Constant gas relative permeability  ---
!
        IF( IRPG(IZN).EQ.0 ) THEN
          RKGX = RPGC(3,IZN)
!
!---    Mualem gas relative permeability  ---
!
        ELSEIF( IRPG(IZN).EQ.1 ) THEN
          CL = RPGC(3,IZN)
          RKGX = SQRT(ASGX)*(1.D+0-ASTX**(1.D+0+1.D+0/CL))**2
!
!---    Burdine gas relative permeability  ---
!
        ELSEIF( IRPG(IZN).EQ.2 ) THEN
          CL = RPGC(3,IZN)
          RKGX = (ASGX**2)*(1.D+0-ASTX**(1.D+0 + 2.D+0/CL))
!
!---    Free-Corey relative permeability function,
!       modified from STOMP-CO2 for the Moridis
!       code comparison problem  ---
!
        ELSEIF( IRPG(IZN).EQ.7 ) THEN
          SLRX = RPGC(3,IZN)
          SGRX = RPGC(4,IZN)
          SGPX = MIN( MAX( (SGX-SGRX)/(1.D+0-SLRX),0.D+0 ),
     &      1.D+0 )
          RKGX = RPGC(1,IZN)*(SGPX**(RPGC(2,IZN)))
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
      ENDIF
      RKLX(1) = MAX( RKLX(1),0.D+0 )
      RKLX(2) = MAX( RKLX(1),0.D+0  )
      RKLX(3) = MAX( RKLX(1),0.D+0  )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KSP37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Load the current time step values into the old time step
!     variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      USE FDVGC
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
      SUB_LOG(ISUB_LOG) = '/LDO37'
!
!---  Assign old time step values  ---
!
      DO 100 N = 1,NFLD
        ASLMIN(1,N) = MIN( ASL(N),ASLMIN(2,N) )
        ASLMIN(2,N) = ASLMIN(1,N)
        BTGL(1,N) = BTGL(2,N)
        DFGA(1,N) = DFGA(2,N)
        DFGO(1,N) = DFGO(2,N)
        DFLA(1,N) = DFLA(2,N)
        DFLO(1,N) = DFLO(2,N)
        DFNC(1,1,N) = DFNC(1,2,N)
        DFNC(2,1,N) = DFNC(2,2,N)
        DFNC(3,1,N) = DFNC(3,2,N)
        NPHAZ(1,N) = NPHAZ(2,N)
        PG(1,N) = PG(2,N)
        PH(1,N) = PH(2,N)
        PI(1,N) = PI(2,N)
        PL(1,N) = PL(2,N)
        PN(1,N) = PN(2,N)
        PORD(1,N) = PORD(2,N)
        PORT(1,N) = PORT(2,N)
        PVA(1,N) = PVA(2,N)
        PVO(1,N) = PVO(2,N)
        PVW(1,N) = PVW(2,N)
        RHOG(1,N) = RHOG(2,N)
        RHOH(1,N) = RHOH(2,N)
        RHOI(1,N) = RHOI(2,N)
        RHOL(1,N) = RHOL(2,N)
        RHON(1,N) = RHON(2,N)
        RKG(1,N) = RKG(2,N)
        RKL(1,1,N) = RKL(1,2,N)
        RKL(2,1,N) = RKL(2,2,N)
        RKL(3,1,N) = RKL(3,2,N)
        TMS(1,N) = TMS(2,N)
        SH(1,N) = SH(2,N)
        SI(1,N) = SI(2,N)
        SG(1,N) = SG(2,N)
        SL(1,N) = SL(2,N)
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
        XGW(1,N) = XGW(2,N)
        XHA(1,N) = XHA(2,N)
        XHO(1,N) = XHO(2,N)
        XHW(1,N) = XHW(2,N)
        XLA(1,N) = XLA(2,N)
        XLO(1,N) = XLO(2,N)
        XLS(1,N) = XLS(2,N)
        XLW(1,N) = XLW(2,N)
        XNC(1,1,N) = XNC(1,2,N)
        XNC(2,1,N) = XNC(2,2,N)
        XNC(3,1,N) = XNC(3,2,N)
        XNA(1,N) = XNA(2,N)
        XNO(1,N) = XNO(2,N)
        XNW(1,N) = XNW(2,N)
        XMGA(1,N) = XMGA(2,N)
        XMGO(1,N) = XMGO(2,N)
        XMGW(1,N) = XMGW(2,N)
        XMLA(1,N) = XMLA(2,N)
        XMLO(1,N) = XMLO(2,N)
        XMLW(1,N) = XMLW(2,N)
        XMNC(1,1,N) = XMNC(1,2,N)
        XMNC(2,1,N) = XMNC(2,2,N)
        XMNC(3,1,N) = XMNC(3,2,N)
        YLS(1,N) = YLS(2,N)
        YMGO(1,N) = YMGO(2,N)
#ifdef k1_hyd
        VISDG(2,N) = VISDG(1,N)
#endif
#ifdef k2_hyd
        VISDG(2,N) = VISDG(1,N)
        VISDG(4,N) = VISDG(3,N)
#endif
!
!---    Isothermal option  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          HGA(1,N) = HGA(2,N)
          HGO(1,N) = HGO(2,N)
          HGW(1,N) = HGW(2,N)
          HG(1,N) = HG(2,N)
          HLA(1,N) = HLA(2,N)
          HLO(1,N) = HLO(2,N)
          HLS(1,N) = HLS(2,N)
          HLW(1,N) = HLW(2,N)
          HL(1,N) = HL(2,N)
          HH(1,N) = HH(2,N)
          HI(1,N) = HI(2,N)
          HN(1,N) = HN(2,N)
          THKL(1,N) = THKL(2,N)
          THKG(1,N) = THKG(2,N)
          THKH(1,N) = THKH(2,N)
          THKI(1,N) = THKI(2,N)
          UEG(1,N) = UEG(2,N)
        ENDIF
        DO 90 NSL = 1,NSOLU
          CO(N,NSL) = C(N,NSL)
   90   CONTINUE
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO37 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Compute hydrologic, thermodynamic and physical properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 16 October 2007
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVN
      USE FDVI
      USE FDVH
      USE FDVGC
      USE FDVG
      USE FDVA
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Data Statements---------------------------------!
!
      SAVE HFLH
      DATA HFLH / 5.4D+7 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROP37'
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 200
        IZN = IZ(N)
        N_DB = N
!
!---    Minimum apparent aqueous saturation  ---
!
        ASLMINX = ASLMIN(1,N)
!
!---    Loop over increment indices  ---
!
        DO 100 M = 2,ISVC+2
!
!---      Saturated CO2 pressure  ---
!
          PLX = PL(M,N)+PATM
          PGX = PG(M,N)+PATM
          CALL SP_A( T(M,N),PSAX )
!
!---      Aqueous or aqueous-gas, no hydrate, no liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - CO2 vapor pressure
!         CH4 mass - CH4 vapor pressure
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          IF( NPHAZ(2,N).EQ.1 ) THEN
            PSBX = 1.D+0
            CALL FLH_PC1( PEQHX,PG(M,N),PH(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PSBX,PVA(M,N),PVO(M,N),PVW(M,N),RHOH(M,N),
     &        RKG(M,N),RKL(1,M,N),RKN(M,N),SH(M,N),SI(M,N),
     &        SG(M,N),SL(M,N),SN(M,N),T(M,N),TEQHX,XHA(M,N),XHO(M,N),
     &        XHW(M,N),XLS(M,N),YLS(M,N),YMGO(M,N),IZN,N )
!
!---      Aqueous or aqueous-gas, hydrate, no liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - Gas CH4 mole fraction of formers
!         CH4 mass - hydrate saturation
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
            PSBX = 1.D+0
#ifdef eq_hyd
            CALL FLH_PC2( PEQHX,PG(M,N),PH(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PSBX,PVA(M,N),PVO(M,N),PVW(M,N),RHOH(M,N),RKG(M,N),
     &        RKL(1,M,N),RKN(M,N),SH(M,N),SI(M,N),SG(M,N),
     &        SL(M,N),SN(M,N),T(M,N),TEQHX,XHA(M,N),XHO(M,N),XHW(M,N),
     &        XLS(M,N),YLS(M,N),YMGO(M,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC2( PEQHX,PG(M,N),PH(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PSBX,PVA(M,N),PVO(M,N),PVW(M,N),RHOH(M,N),RKG(M,N),
     &        RKL(1,M,N),RKN(M,N),SH(M,N),SI(M,N),SG(M,N),
     &        SL(M,N),SN(M,N),T(M,N),TEQHX,XHA(M,N),XHO(M,N),XHW(M,N),
     &        XLS(M,N),YLS(M,N),YMGO(M,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC2( PEQHX,PG(M,N),PH(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PSBX,PVA(M,N),PVAHX,PVAHOX,PVO(M,N),PVOHX,PVOHOX,
     &        PVW(M,N),RHOH(M,N),RKG(M,N),
     &        RKL(1,M,N),RKN(M,N),SH(M,N),SI(M,N),SG(M,N),
     &        SL(M,N),SN(M,N),T(M,N),TEQHX,XHA(M,N),XHO(M,N),XHW(M,N),
     &        XLS(M,N),YLS(M,N),YMGO(M,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
!
!---      Aqueous or aqueous-gas, no hydrate, low liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - liquid CO2 saturation
!         CH4 mass - CH4 vapor pressure
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
            PVA(M,N) = PSAX
            PSBX = 1.D+0
            CALL FLH_PC3( PEQHX,PG(M,N),PH(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PSBX,PVA(M,N),PVO(M,N),PVW(M,N),RHOH(M,N),
     &        RKG(M,N),RKL(1,M,N),RKN(M,N),SH(M,N),SI(M,N),
     &        SG(M,N),SL(M,N),SN(M,N),T(M,N),TEQHX,XHA(M,N),XHO(M,N),
     &        XHW(M,N),XLS(M,N),YLS(M,N),YMGO(M,N),IZN,N )
!
!---      Aqueous or aqueous-gas, hydrate, low liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - liquid CO2 saturation
!         CH4 mass - hydrate saturation
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
            PSBX = 1.D+0
            CALL FLH_PC4( PEQHX,PG(M,N),PH(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PSBX,PVA(M,N),PVO(M,N),PVW(M,N),
     &        RHOH(M,N),RKG(M,N),RKL(1,M,N),RKN(M,N),SH(M,N),
     &        SI(M,N),SG(M,N),SL(M,N),SN(M,N),T(M,N),TEQHX,XHA(M,N),
     &        XHO(M,N),XHW(M,N),XLS(M,N),YLS(M,N),YMGO(M,N),IZN,N )
!
!---      Aqueous or aqueous-gas, no hydrate, liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - liquid CO2 pressure
!         CH4 mass - CH4 vapor pressure
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
            PVA(M,N) = PSAX
            PSBX = 1.D+0
            CALL FLH_PC5( PEQHX,PG(M,N),PH(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PSBX,PVA(M,N),PVO(M,N),PVW(M,N),RHOH(M,N),
     &        RKG(M,N),RKL(1,M,N),RKN(M,N),SH(M,N),SI(M,N),
     &        SG(M,N),SL(M,N),SN(M,N),T(M,N),TEQHX,XHA(M,N),XHO(M,N),
     &        XHW(M,N),XLS(M,N),YLS(M,N),YMGO(M,N),IZN,N )
!
!---      Aqueous or aqueous-gas, hydrate, liquid CO2
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         CO2 mass - liquid CO2 pressure
!         CH4 mass - hydrate saturation
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
            PSBX = 1.D+0
#ifdef eq_hyd
            CALL FLH_PC6( PEQHX,PG(M,N),PH(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PSBX,PVA(M,N),PVO(M,N),PVW(M,N),RHOH(M,N),RKG(M,N),
     &        RKL(1,M,N),RKN(M,N),SH(M,N),SI(M,N),SG(M,N),
     &        SL(M,N),SN(M,N),T(M,N),TEQHX,XHA(M,N),XHO(M,N),XHW(M,N),
     &        XLS(M,N),YLS(M,N),YMGO(2,N),IZN,N )
#endif
#ifdef k1_hyd
            YMGOOX = VISDG(2,N)
            CALL FLH_PC6( PEQHX,PG(M,N),PH(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PSBX,PVA(M,N),PVO(M,N),PVW(M,N),RHOH(M,N),RKG(M,N),
     &        RKL(1,M,N),RKN(M,N),SH(M,N),SI(M,N),SG(M,N),
     &        SL(M,N),SN(M,N),T(M,N),TEQHX,XHA(M,N),XHO(M,N),XHW(M,N),
     &        XLS(M,N),YLS(M,N),YMGO(2,N),YMGOY,YMGOOX,IZN,N )
            VISDG(1,N) = YMGOY
#endif
#ifdef k2_hyd
            PVOHOX = VISDG(2,N)
            PVAHOX = VISDG(4,N)
            CALL FLH_PC6( PEQHX,PG(M,N),PH(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PSBX,PVA(M,N),PVAHX,PVAHOX,PVO(M,N),PVOHX,PVOHOX,
     &        PVW(M,N),RHOH(M,N),RKG(M,N),
     &        RKL(1,M,N),RKN(M,N),SH(M,N),SI(M,N),SG(M,N),
     &        SL(M,N),SN(M,N),T(M,N),TEQHX,XHA(M,N),XHO(M,N),XHW(M,N),
     &        XLS(M,N),YLS(M,N),YMGO(2,N),IZN,N )
            VISDG(1,N) = PVOHX
            VISDG(3,N) = PVAHX
#endif
          ENDIF
!
!---      Absolute system and gas pressure  ---
!
          PX = MAX( PG(M,N),PL(M,N),PN(M,N) ) + PATM
          PGX = PG(M,N) + PATM
!
!---      Aqueous component fractions and density  ---
!
          INCG = 1
          CALL SOL_LA( T(M,N),PVA(M,N),XLS(M,N),XMLA(M,N) )
          INCG = 2
          CALL SOL_LO( T(M,N),PVO(M,N),XLS(M,N),XMLO(M,N) )
          WTMLX = (XMLO(M,N)*WTMO + XMLA(M,N)*WTMA + 
     &      WTMW*(1.D+0-XMLO(M,N)-XMLA(M,N)))/
     &      (1.D+0-XLS(M,N)+WTMW*XLS(M,N)/WTMS)
          XLA(M,N) = XMLA(M,N)*WTMA/WTMLX
          XLO(M,N) = XMLO(M,N)*WTMO/WTMLX
          CALL DENS_B( T(M,N),PX,XLS(M,N),RHOBX )
          CALL DENS_L( T(M,N),RHOBX,XLA(M,N),XLO(M,N),RHOL(M,N) )
          XLW(M,N) = MAX( 1.D+0-XLA(M,N)-XLO(M,N)-XLS(M,N),0.D+0 )
          XMLS(M,N) = WTMLX*XLS(M,N)/WTMS
          XMLW(M,N) = WTMLX*XLW(M,N)/WTMW
!
!---      Aqueous molar density  ---
!
          RHOML(M,N) = RHOL(M,N)/WTMLX
!
!---      Ice density  ---
!
          CALL DENS_I( T(M,N),RHOI(M,N) )
!
!---      Porous-media porosity  ---
!
          CALL PORSTY( N,PX,PCMP(N),PORD(M,N),PORT(M,N) )
          PORD(M,N) = MAX( PORD(M,N),EPSL )
          PORT(M,N) = MAX( PORT(M,N),PORD(M,N) )
!
!---      Gas density and component fractions  ---
!
          ISRX = 2
          CALL DENS_W( T(M,N),PVW(M,N),RHOX,RHOGWX,ISRX )
          CALL DENS_A( T(M,N),PVA(M,N),RHOGAX,ISRX )
          CALL DENS_O( T(M,N),PVO(M,N),RHOGOX )
          RHOG(M,N) = RHOGWX+RHOGAX+RHOGOX
          XGA(M,N) = RHOGAX/RHOG(M,N)
          XGO(M,N) = RHOGOX/RHOG(M,N)
          XGW(M,N) = RHOGWX/RHOG(M,N)
          WTMGX = XGA(M,N)/WTMA + XGO(M,N)/WTMO + XGW(M,N)/WTMW
          XMGA(M,N) = XGA(M,N)/WTMA/WTMGX
          XMGO(M,N) = XGO(M,N)/WTMO/WTMGX
          XMGW(M,N) = XGW(M,N)/WTMW/WTMGX
          RHOMG(M,N) = RHOG(M,N)*WTMGX
!
!---      Gas viscosity  ---
!
          CALL VISC_A( T(M,N),RHOGAX,VISGAX )
          CALL VISC_O( T(M,N),RHOGOX,VISGOOX,VISGOEX,VISGOX )
          CALL VISC_W( T(M,N),PVW(M,N),RHOGWX,VISGWX )
          CALL VISC_G( VISGAX,VISGOX,VISGWX,XMGA(M,N),
     &      XMGO(M,N),XMGW(M,N),VISG(M,N) )
!
!---      Gas diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            DFGO(M,N) = DFGOC
            DFGA(M,N) = DFGAC
          ELSEIF( ISLC(2).EQ.2 ) THEN
            CALL BNDFAW( T(M,N),PGX,DFGAW )
            CALL BNDFAO( T(M,N),PGX,DFGAO )
            CALL BNDFOW( T(M,N),PGX,DFGOW )
            IF( (1.D+0-XMGO(M,N)).LT.EPSL ) THEN
              DFGO(M,N) = 2.D+0*DFGAO*DFGOW/(DFGAO+DFGOW)
            ELSE
              DFGO(M,N) = (1.D+0-XMGO(M,N))/
     &          ((XMGA(M,N)/DFGAO) + (XMGW(M,N)/DFGOW))
            ENDIF
            IF( (1.D+0-XMGA(M,N)).LT.EPSL ) THEN
              DFGA(M,N) = 2.D+0*DFGAO*DFGAW/(DFGAO+DFGAW)
            ELSE
              DFGA(M,N) = (1.D+0-XMGA(M,N))/
     &          ((XMGO(M,N)/DFGAO) + (XMGW(M,N)/DFGAW))
            ENDIF
          ENDIF
!
!---      Aqueous viscosity  ---
!
          ISRX = 1
          CALL DENS_W( T(M,N),PX,RHOLWX,RHOX,ISRX )
          CALL VISC_W( T(M,N),PX,RHOLWX,VISLWX )
          CALL VISC_B( T(M,N),XLS(M,N),VISLWX,VISBX )
          CALL VISC_L( XMLA(M,N),XMLO(M,N),VISBX,VISGAX,VISGOX,
     &      VISL(M,N) )
!
!---      Dissolved-CO2 and -NaCl diffusion coefficient  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLA(M,N) = DFLAC
            DFLO(M,N) = DFLOC
            DFLS(M,N) = DFLSC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            INCG = 1
            CALL DIFC_LA( T(M,N),VISL(M,N),VISGAX,DFLA(M,N) )
            INCG = 2
            CALL DIFC_LA( T(M,N),VISL(M,N),VISGOX,DFLO(M,N) )
            CALL DIFC_LS( T(M,N),XLS(M,N),VISL(M,N),DFLS(M,N) )
          ENDIF
!
!---      Liquid-CO2 density and viscosity  ---
!
          ISRX = 1
          CALL DENS_A( T(M,N),PX,RHON(M,N),ISRX )
          CALL VISC_A( T(M,N),RHON(M,N),VISN(M,N) )
          CALL SOL_N( PX,RHOG(M,N),RHOMG(M,N),RHON(M,N),RHOMN(M,N),
     &      T(M,N),WTMNX,XMGA(M,N),XMGO(M,N),XMGW(M,N),XNC(1,M,N),
     &      XNC(2,M,N),XNC(3,M,N) )
          XMNC(1,M,N) = XNC(1,M,N)/WTMA/WTMNX
          XMNC(2,M,N) = XNC(2,M,N)/WTMO/WTMNX
          XMNC(3,M,N) = XNC(3,M,N)/WTMW/WTMNX
          XNA(M,N) = XNC(1,M,N)
          XNO(M,N) = XNC(2,M,N)
          XNW(M,N) = XNC(3,M,N)
          CALL DIFC_N( DFNC(1,M,N),T(M,N),VISN(M,N),WTMNX )
!  
!---      Aqueous thermal conductivity  ---
!
          CALL THK_W( T(M,N),PX,RHOLWX,THKLWX )
          CALL THK_B( T(M,N),XLS(M,N),THKLWX,THKL(M,N) )
!
!---      Gas thermal conductivity  ---
!
          CALL THK_A( T(M,N),PGX,RHOGAX,THKGAX )
          CALL THK_O( T(M,N),RHOGOX,THKGOX,VISGOOX )
          CALL THK_W( T(M,N),PGX,RHOGWX,THKGWX )
          CALL THK_G( T(M,N),THKGAX,THKGOX,THKGWX,
     &      XMGA(M,N),XMGO(M,N),XMGW(M,N),THKG(M,N) )
!
!---      Liquid-CO2 thermal conductivity  ---
!
          CALL THK_A( T(M,N),PX,RHON(M,N),THKN(M,N) )
!
!---      Ice thermal conductivity  ---
!
          CALL THK_I( T(M,N),THKI(M,N) )
!
!---      Hydrate thermal conductivity  ---
!
          CALL THK_H( T(M,N),THKH(M,N) )
!
!---      Gas enthalpy and internal energy  ---
!
          ISRX = 2
          CALL ENTH_A( T(M,N),PVA(M,N),HGA(M,N),UEGA(M,N),ISRX )
          CALL ENTH_O( T(M,N),PVO(M,N),HGO(M,N),UGOX )
          CALL ENTH_W( T(M,N),PVW(M,N),HLWX,HGW(M,N),ISRX )
          UGWX = HGW(M,N) - PVW(M,N)/RHOGWX
          UEG(M,N) = XGA(M,N)*UEGA(M,N)+XGO(M,N)*UGOX+XGW(M,N)*UGWX
          HG(M,N) = XGA(M,N)*HGA(M,N) + XGO(M,N)*HGO(M,N) + 
     &      XGW(M,N)*HGW(M,N)
!
!---      Liquid water enthalpy  ---
!
          ISRX = 1
          CALL ENTH_W( T(M,N),PX,HLW(M,N),HGWX,ISRX )
!
!---      Liquid-CO2 enthalpy  ---
!
          ISRX = 1
          CALL ENTH_A( T(M,N),PX,HNAX,UNX,ISRX )
          FGWX = SG(M,N)*RHOG(M,N)*XGW(M,N) + SMALL
          FLWX = SL(M,N)*RHOL(M,N)*XLW(M,N) + SMALL
          HNWX = (HGW(M,N)*FGWX + HLW(M,N)*FLWX)/(FGWX + FLWX)
          HN(M,N) = XNA(M,N)*HNAX + XNW(M,N)*HNWX + 
     &      XNO(M,N)*HGO(M,N)
!
!---      Aqueous enthalpy  ---
!
          CALL ENTH_B( T(M,N),XLS(M,N),HLW(M,N),HBX )
          FGAX = SG(M,N)*RHOG(M,N)*XGA(M,N) + SMALL
          FNAX = SN(M,N)*RHON(M,N)*XNA(M,N) + SMALL
          HGAX = (HGA(M,N)*FGAX + HN(M,N)*FNAX)/(FGAX + FNAX)
          CALL ENTH_L( T(M,N),XLS(M,N),XLA(M,N),XLO(M,N),HBX,
     &      HGAX,HGO(M,N),HL(M,N) )
!
!---      Ice enthalpy  ---
!
          CALL ENTH_I( T(M,N),HI(M,N) )
!
!---      Hydrate enthalpy  ---
!
          CALL ENTH_H( T(M,N),XHA(M,N),XHO(M,N),XHW(M,N),
     &      HGA(M,N),HGO(M,N),HH(M,N),HLW(M,N) )
!
!---      Precipitated NaCl density, saturation, and enthalpy  ---
!
          CALL DENS_S( T(M,N),PX,RHOSP(M,N) )
          SS(M,N) = MAX(YLS(M,N)-XLS(M,N),0.D+0)*RHOL(M,N)*SL(M,N)/
     &      RHOSP(M,N)
!
!---      Permeability reduction factor  ---
!
          IF( IPRF(IZN).EQ.1 ) CALL PERM_R( SS(M,N),PERMRF(M,N),N )
!
!---      Change in porosity as a function of precipitated salt  ---
!
          PORD(M,N) = PORD(M,N)*(1.D+0-SS(M,N))
!
!---      Aqueous and gas tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL(M,N),SG(M,N),SN(M,N),
     &      PORD(M,N),TORL(M,N),TORG(M,N),TORN(M,N) )
!
!---      NaCl volumetric concentration  ---
!
          TMS(M,N) = YLS(M,N)*RHOL(M,N)*SL(M,N)*PORD(M,N)
  100   CONTINUE
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDBC37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Read input file for boundary condition information.
!
!     1 - Dirichlet
!     2 - Neumann
!     3 - Zero Flux
!     4 - Saturated
!     5 - Unit Gradient
!     6 - Free Gradient
!     7 - Outflow
!     8 - Aqueous Concentration
!     9 - Gas Concentration
!     10 - NAPL Concentration
!     11 - Hydraulic Gradient
!     12 - Initial Condition
!     13 - Inflow
!     14 - Inflow Aqueous-Phase
!     15 - Inflow Gas-Phase
!     16 - Inflow NAPL
!     17 - Seepage Face
!     18 - Convective
!     19 - Inflow-Outflow Volumetric
!     20 - Falling Head
!     21 - Falling Pond
!     22 - Free Boundary
!     23 - Inflow-Outflow Aqueous
!     24 - Potential Evaporation
!     25 - Fluctuating Water Table
!     26 - Dirichlet-Outflow
!     27 - Diode
!     28 - Convective-Radiative
!     29 - Convective Ground Surface
!     30 - Shuttleworth-Wallace
!     31 - Bare Shuttleworth-Wallace
!     32 - Relative Saturation
!     33 - Inflow Relative Saturation
!     34 - Aqu. Rel. Sat.
!     35 - Inflow Aqu. Rel. Sat.
!     36 - Aqu. Mass Frac.
!     37 - Inflow Aqu. Mass Frac.
!     38 - Vol. Conc.
!     39 - Inflow Vol. Conc.
!     40 - Aqu. Conc.
!     41 - Inflow Aqu. Conc.
!     42 - Dirichlet-Inflow
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      CHARACTER*32 CHTYP(42)
      CHARACTER*512 CHDUM
      REAL*8 VAR(LBTM,LBCV)
      INTEGER ITYP(LUK+LSOLU)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CHTYP
      DATA CHTYP /'Dirichlet','Neumann','Zero Flux','Saturated',
     &  'Unit Gradient','Free Gradient','Outflow',
     &  'Aqueous Concentration','Gas Concentration',
     &  'NAPL Concentration','Hydraulic Gradient',
     &  'Initial Condition','Inflow','Inflow Aqueous-Phase',
     &  'Inflow Gas-Phase','Inflow NAPL','Seepage Face','Convective',
     &  'Inflow-Outflow Volumetric','Falling Head','Falling Pond',
     &  'Free Boundary','Inflow-Outflow Aqueous',
     &  'Potential Evaporation','Fluctuating Water Table',
     &  'Dirichlet-Outflow','Diode',
     &  'Convective-Radiative','Convective Ground Surface',
     &  'Shuttleworth-Wallace','Bare Shuttleworth-Wallace',
     &  'Relative Saturation','Inflow Relative Saturation',
     &  'Aqu. Rel. Sat.','Inflow Aqu. Rel. Sat.','Aqu. Mass Frac.',
     &  'Inflow Aqu. Mass Frac.','Vol. Conc.','Inflow Vol. Conc.',
     &  'Aqu. Conc.','Inflow Aqu. Conc.','Dirichlet-Inflow'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDBC37'
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
!
!---    Isothermal option  ---
!
        IF( ISLC(30).EQ.0 )
     &    CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(2))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(3))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(4))
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 )
     &    CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(5))
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
!---    Isothermal option  ---
!
        IF( ISLC(30).EQ.0 ) THEN
!
!---      Read energy boundary condition type ---
!
          IF( INDEX(BDUM(1)(1:),'dirichlet').NE.0 ) THEN
            ITYP(1) = 1
          ELSEIF( INDEX(BDUM(1)(1:),'neumann').NE.0 ) THEN
            ITYP(1) = 2
          ELSEIF( INDEX(BDUM(1)(1:),'zero flux').NE.0 ) THEN
            ITYP(1) = 3
          ELSEIF( INDEX(BDUM(1)(1:),'outflow').NE.0 ) THEN
            ITYP(1) = 7
          ELSEIF( INDEX(BDUM(1)(1:),'initial cond').NE.0 ) THEN
            ITYP(1) = 12
          ELSEIF( INDEX(BDUM(1)(1:),'inflow').NE.0 ) THEN
            ITYP(1) = 13
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Energy Boundary Condition: ' //
     &        BDUM(1)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  Read aqueous boundary condition type ---
!
        IF( INDEX(BDUM(2)(1:),'dirichlet-outflow').NE.0 ) THEN
          ITYP(2) = 26
        ELSEIF( INDEX(BDUM(2)(1:),'dirichlet-inflow').NE.0 ) THEN
          ITYP(2) = 42
        ELSEIF( INDEX(BDUM(2)(1:),'dirichlet').NE.0 ) THEN
          IF( INDEX(BDUM(2)(1:),'inflow').NE.0 .OR.
     &      INDEX(BDUM(2)(1:),'outflow').NE.0 ) THEN
            INDX = 4
            CHMSG = 'Ambiguous Aqueous Boundary Condition: ' //
     &        BDUM(2)
            CALL WRMSGS( INDX )
          ENDIF
          ITYP(2) = 1
        ELSEIF( INDEX(BDUM(2)(1:),'neumann').NE.0 ) THEN
          ITYP(2) = 2
        ELSEIF( INDEX(BDUM(2)(1:),'zero flux').NE.0 ) THEN
          ITYP(2) = 3
        ELSEIF( INDEX(BDUM(2)(1:),'unit gradient').NE.0 ) THEN
          ITYP(2) = 5
        ELSEIF( INDEX(BDUM(2)(1:),'hydraulic gradient').NE.0 ) THEN
          ITYP(2) = 11
        ELSEIF( INDEX(BDUM(2)(1:),'initial cond').NE.0 ) THEN
          ITYP(2) = 12
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Aqueous Boundary Condition: ' //
     &      BDUM(2)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Read gas boundary condition type ---
!
        IF( INDEX(BDUM(3)(1:),'dirichlet-outflow').NE.0 ) THEN
          ITYP(3) = 26
        ELSEIF( INDEX(BDUM(3)(1:),'dirichlet-inflow').NE.0 ) THEN
          ITYP(3) = 42
        ELSEIF( INDEX(BDUM(3)(1:),'dirichlet').NE.0 ) THEN
          IF( INDEX(BDUM(2)(1:),'inflow').NE.0 .OR.
     &      INDEX(BDUM(2)(1:),'outflow').NE.0 ) THEN
            INDX = 4
            CHMSG = 'Ambiguous Gas Boundary Condition: ' //
     &        BDUM(2)
            CALL WRMSGS( INDX )
          ENDIF
          ITYP(3) = 1
        ELSEIF( INDEX(BDUM(3)(1:),'neumann').NE.0 ) THEN
          ITYP(3) = 2
        ELSEIF( INDEX(BDUM(3)(1:),'zero flux').NE.0 ) THEN
          ITYP(3) = 3
        ELSEIF( INDEX(BDUM(3)(1:),'hydraulic gradient').NE.0 ) THEN
          ITYP(3) = 11
        ELSEIF( INDEX(BDUM(3)(1:),'initial cond').NE.0 ) THEN
          ITYP(3) = 12
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Gas Boundary Condition: ' //
     &      BDUM(3)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Read liquid-CO2 boundary condition type ---
!
        IF( INDEX(BDUM(4)(1:),'dirichlet-outflow').NE.0 ) THEN
          ITYP(4) = 26
        ELSEIF( INDEX(BDUM(4)(1:),'dirichlet-inflow').NE.0 ) THEN
          ITYP(4) = 42
        ELSEIF( INDEX(BDUM(4)(1:),'dirichlet').NE.0 ) THEN
          IF( INDEX(BDUM(2)(1:),'inflow').NE.0 .OR.
     &      INDEX(BDUM(2)(1:),'outflow').NE.0 ) THEN
            INDX = 4
            CHMSG = 'Ambiguous Liquid-CO2 Boundary Condition: ' //
     &        BDUM(2)
            CALL WRMSGS( INDX )
          ENDIF
          ITYP(4) = 1
        ELSEIF( INDEX(BDUM(4)(1:),'neumann').NE.0 ) THEN
          ITYP(4) = 2
        ELSEIF( INDEX(BDUM(4)(1:),'zero flux').NE.0 ) THEN
          ITYP(4) = 3
        ELSEIF( INDEX(BDUM(4)(1:),'hydraulic gradient').NE.0 ) THEN
          ITYP(4) = 11
        ELSEIF( INDEX(BDUM(4)(1:),'initial cond').NE.0 ) THEN
          ITYP(4) = 12
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Liquid-CO2 Boundary Condition: ' //
     &      BDUM(4)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read salt boundary condition type ---
!
        IF( ISLC(32).EQ.0 ) THEN
          IF( INDEX(BDUM(5)(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(5)(1:),'rel').NE.0 .AND.
     &      INDEX(BDUM(5)(1:),'sat').NE.0 ) THEN
            ITYP(IEQS) = 34
          ELSEIF( INDEX(BDUM(5)(1:),'inflow').NE.0 .AND.
     &      (INDEX(BDUM(5)(1:),'rel').NE.0 .OR.
     &      INDEX(BDUM(5)(1:),'sat').NE.0) ) THEN
            ITYP(IEQS) = 35
          ELSEIF( INDEX(BDUM(5)(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(5)(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(5)(1:),'frac').NE.0 ) THEN
             ITYP(IEQS) = 36
          ELSEIF( INDEX(BDUM(5)(1:),'inflow').NE.0 .AND.
     &      (INDEX(BDUM(5)(1:),'mass').NE.0 .OR.
     &      INDEX(BDUM(5)(1:),'frac').NE.0) ) THEN
            ITYP(IEQS) = 37
          ELSEIF( INDEX(BDUM(5)(1:),'zero flux').NE.0 ) THEN
             ITYP(IEQS) = 3
          ELSEIF( INDEX(BDUM(5)(1:),'outflow').NE.0 ) THEN
             ITYP(IEQS) = 7
          ELSEIF( INDEX(BDUM(5)(1:),'initial cond').NE.0 ) THEN
             ITYP(IEQS) = 12
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Salt Boundary Condition: '//BDUM(5)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  Read solute boundary condition type(s) ---
!
        IF( IEQC.GT.0 ) THEN
          DO 25 NSL = 1,NSOLU
            IF(INDEX(BDUM(NSL+LUK)(1:),'volumetric conc').NE.0 ) THEN
               ITYP(NSL+LUK) = 1
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'aqueous conc').NE.0 )THEN
               ITYP(NSL+LUK) = 8
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'gas conc').NE.0 ) THEN
               ITYP(NSL+LUK) = 9
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'outflow').NE.0 ) THEN
               ITYP(NSL+LUK) = 7
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'zero flux').NE.0 ) THEN
               ITYP(NSL+LUK) = 3
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'initial cond').NE.0 ) THEN
               ITYP(NSL+LUK) = 12
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Solute Boundary Condition: ' //
     &          BDUM(NSL+LUK)
              CALL WRMSGS( INDX )
            ENDIF
   25     CONTINUE
        ENDIF
!
!---  Write boundary condition type(s) ---
!
        WRITE(IWR,'(A)') 'Boundary Condition Type: '
!
!---    Isothermal option  ---
!
        IF( ISLC(30).EQ.0 )
     &    WRITE(IWR,'(2X,2A)') 'Energy: ',CHTYP(ITYP(1))
        WRITE(IWR,'(2X,2A)') 'Aqueous: ',CHTYP(ITYP(2))
        WRITE(IWR,'(2X,2A)') 'Gas: ',CHTYP(ITYP(3))
        WRITE(IWR,'(2X,2A)') 'Liquid-CO2: ',CHTYP(ITYP(4))
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 )
     &    WRITE(IWR,'(2X,2A)') 'Salt: ',CHTYP(ITYP(IEQS))
        IF( IEQC.GT.0 ) THEN
          DO 30 NSL = 1,NSOLU
            IDB = INDEX( SOLUT(NSL)(1:),'  ') - 1
            WRITE(IWR,'(2X,2A)') SOLUT(NSL)(1:IDB),CHTYP(ITYP(NSL+LUK))
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
!---  Read, write, and convert boundary condition time, variables,
!     and units  ---
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
!---      Isothermal option  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Dirichlet, outflow, or inflow energy ---
!
            IF( ITYP(1).EQ.1 .OR. ITYP(1).EQ.7 .OR.
     &        ITYP(1).EQ.13 ) THEN
              VARB = 'Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
!
!---        Neumann energy ---
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
!---        Zero flux energy ---
!
            ELSEIF( ITYP(1).EQ.3 ) THEN
              VARB = 'Zero Flux'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
!
!---        Initial condition energy ---
!
            ELSE
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            ENDIF
          ENDIF
!
!---      Dirichlet aqueous, Dirichlet-inflow aqueous, or
!         Dirichlet-outflow aqueous  ---
!
          IF( ITYP(2).EQ.1 .OR. ITYP(2).EQ.26
     &      .OR. ITYP(2).EQ.42 ) THEN
            VARB = 'Aqueous Pressure'
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
!---      Neumann aqueous  ---
!
          ELSEIF( ITYP(2).EQ.2 ) THEN
            VARB = 'Volumetric Aqueous Flux'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
!
!---      Zero-flux aqueous  ---
!
          ELSEIF( ITYP(2).EQ.3 ) THEN
            VARB = 'Aqueous Pressure'
            IDFLT = 1
            VAR(NTM,3) = -1.D+20
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
!---      Hydraulic-gradient aqueous  ---
!
          ELSEIF( ITYP(2).EQ.11 ) THEN
            VARB = 'Base Aqueous Pressure'
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
!---      Initial-condition aqueous, saturation aqueous, 
!         or unit-gradient aqueous  ---
!
          ELSE
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          ENDIF
!
!---      Read dissolved-CO2 aqueous relative saturation ---
!
          VARB = 'Dissolved-CO2 Aqueous Relative Saturation'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &      ': ',VAR(NTM,7)
!
!---      Read dissolved-CH4 aqueous relative saturation ---
!
          VARB = 'Dissolved-CH4 Aqueous Relative Saturation'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &      ': ',VAR(NTM,8)
          
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Read dissolved-NaCl aqueous relative saturation ---
!
            IF( ITYP(IEQS).EQ.34 .OR. ITYP(IEQS).EQ.35 ) THEN
              VARB = 'Dissolved-NaCl Aqueous Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,6)
!
!---        Read dissolved-NaCl aqueous mass fraction ---
!
            ELSEIF( ITYP(IEQS).EQ.36 .OR. ITYP(IEQS).EQ.37 ) THEN
              VARB = 'Dissolved-NaCl Aqueous Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,6)
!
!---        Initial condition, zero flux or outflow  ---
!
            ELSE
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            ENDIF
          ENDIF
!
!---      Read CO2-slurry aqueous volumetric fraction ---
!
          VARB = 'CO2-Slurry Aqueous Volumetric Fraction'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,10))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &      ': ',VAR(NTM,10)
!
!---      Set dissolved-CO2 aqueous relative saturation
!         for non-zero CO2-slurry fractions  ---
!
          IF( VAR(NTM,10).GT.EPSL ) THEN
            VAR(NTM,7) = 1.D+0
            WRITE(IWR,'(2X,A,1PE11.4)') 
     &        'Overwrite Dissolved-CO2 Aqueous Relative Saturation: ',
     &        VAR(NTM,7)
          ENDIF
!
!---      Dirichlet gas, Dirichlet-inflow gas, or
!         Dirichlet-outflow gas  ---
!
          IF( ITYP(3).EQ.1 .OR. ITYP(3).EQ.26
     &       .OR. ITYP(3).EQ.42 ) THEN
            VARB = 'Gas Pressure'
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
!
!---      Neumann gas  ---
!
          ELSEIF( ITYP(3).EQ.2 ) THEN
            VARB = 'Volumetric Gas Flux'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
!
!---      Zero-flux gas  ---
!
          ELSEIF( ITYP(3).EQ.3 ) THEN
            VARB = 'Gas Pressure'
            IDFLT = 1
            VAR(NTM,4) = -1.D+20
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
!
!---      Hydraulic-gradient gas  ---
!
          ELSEIF( ITYP(3).EQ.11 ) THEN
            VARB = 'Base Gas Pressure'
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
          ELSE
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          ENDIF
!
!---      CO2 gas mass fraction of hydrate formers ---
!
          VARB = 'CO2 Gas Mass Fraction of Hydrate Formers'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,9))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &      ': ',VAR(NTM,9)
!
!---      Read water-vapor relative saturation ---
!
          VARB = 'Water-Vapor Relative Saturation, '
          IDFLT = 1
          VAR(NTM,5) = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &      ': ',VAR(NTM,5)
!
!---      Dirichlet liquid-CO2, Dirichlet-inflow liquid-CO2, or
!         Dirichlet-outflow liquid-CO2  ---
!
          IF( ITYP(4).EQ.1 .OR. ITYP(4).EQ.26
     &      .OR. ITYP(4).EQ.42 ) THEN
            VARB = 'Liquid-CO2 Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,11))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,11)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,11),INDX)
            VAR(NTM,11) = VAR(NTM,11) - PATM
!
!---      Neumann liquid-CO2  ---
!
          ELSEIF( ITYP(4).EQ.2 ) THEN
            VARB = 'Volumetric Liquid-CO2 Flux'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,11))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,11)
            INDX = 0
            IUNM = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,11),INDX)
!
!---      Zero-flux liquid-CO2  ---
!
          ELSEIF( ITYP(4).EQ.3 ) THEN
            VARB = 'Liquid-CO2 Pressure'
            IDFLT = 1
            VAR(NTM,11) = -1.D+20
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,11))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,11)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,11),INDX)
            VAR(NTM,11) = VAR(NTM,11) - PATM
!
!---      Hydraulic gradient liquid-CO2  ---
!
          ELSEIF( ITYP(4).EQ.11 ) THEN
            VARB = 'Base Liquid-CO2 Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,11))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,11)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,11),INDX)
            VAR(NTM,11) = VAR(NTM,11) - PATM
          ELSE
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          ENDIF
!
!---      Hydrate saturation  ---
!
          VARB = 'Hydrate Saturation, '
          IDFLT = 1
          VAR(NTM,12) = 0.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,12))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &      ': ',VAR(NTM,12)
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
              IF( ISLC(30).EQ.0 ) IBCT(1,NBC) = ITYP(1)
              IBCT(2,NBC) = ITYP(2)
              IBCT(3,NBC) = ITYP(3)
              IBCT(4,NBC) = ITYP(4)
              IF( ISLC(32).EQ.0 ) IBCT(IEQS,NBC) = ITYP(IEQS)
              IF( IEQC.GT.0 ) THEN
                DO 110 NSL = 1,NSOLU
                  IBCT(NSL+LUK,NBC) = ITYP(NSL+LUK)
  110           CONTINUE
              ENDIF
              IBCM(NBC) = IBCMX
              IBCIN(NBC) = NB
!
!---  Assign hydraulic gradient type boundary condition
!     for the base surface  ---
!
              IF( ITYP(2).EQ.11 ) THEN
                IF( NBCL.EQ.1 ) THEN
                  IBCT(2,NBC) = -11
                  NBHG(2,NBC) = NBC
                ELSE
                  NBHG(2,NBC) = NBHG(2,NBC-1)
                ENDIF
              ENDIF
              IF( ITYP(3).EQ.11 ) THEN
                IF( NBCL.EQ.1 ) THEN
                  IBCT(3,NBC) = -11
                  NBHG(3,NBC) = NBC
                ELSE
                  NBHG(3,NBC) = NBHG(3,NBC-1)
                ENDIF
              ENDIF
              IF( ITYP(4).EQ.11 ) THEN
                IF( NBCL.EQ.1 ) THEN
                  IBCT(4,NBC) = -11
                  NBHG(4,NBC) = NBC
                ELSE
                  NBHG(4,NBC) = NBHG(4,NBC-1)
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
        IF( INDEX(ADUM(1:),'file').NE.0 ) CLOSE(UNIT=26)
  400 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDBC37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIC37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Read input file for initial conditions information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE NCG_PT
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
      SUB_LOG(ISUB_LOG) = '/RDIC37'
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
        ELSEIF( INDEX( ADUM(1:),'liquid co2 pres' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'liquid-co2 pres' ).NE.0 ) THEN
          VARB = 'Initial Liquid-CO2 Pressure'
          IUNM = -1
          IUNKG = 1
          IUNS = -2
        ELSEIF( INDEX( ADUM(1:),'temperature' ).NE.0 ) THEN
          VARB = 'Initial Temperature'
          IUNK = 1
        ELSEIF( INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Hydrate Saturation'
        ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          VARB = 'Initial Hydrate CO2 Mass Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Hydrate CO2 Mole Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'gas' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          VARB = 'Initial Gas CO2 Mass Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'gas' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Gas CO2 Mole Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'partial' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'press' ).NE.0 ) THEN
          VARB = 'Initial CO2 Partial Pressure'
          IUNM = -1
          IUNKG = 1
          IUNS = -2
        ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'rel' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-CO2 Relative Saturation'
        ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-CO2 Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'conc' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-CO2 Concentration'
          IUNKG = 1
          IUNM = -3
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          VARB = 'Initial Hydrate CH4 Mass Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Hydrate CH4 Mole Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'gas' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          VARB = 'Initial Gas CH4 Mass Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'gas' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Gas CH4 Mole Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'partial' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'press' ).NE.0 ) THEN
          VARB = 'Initial CH4 Partial Pressure'
          IUNM = -1
          IUNKG = 1
          IUNS = -2
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-CH4 Relative Saturation'
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-CH4 Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'conc' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-CH4 Concentration'
          IUNKG = 1
          IUNM = -3
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
          CHMSG = 'Unrecognized Rock/Soil Type: '//FDUM
          CALL WRMSGS( INDX )
          GOTO 1000
   30     CONTINUE
          WRITE(IWR,'(2X,3A,1PE11.4,2A)') ADUM(1:NCHA),UNTS(1:NCHU),
     &      ': ',VAR(1),' Rock/Soil Type: ',FDUM(1:NCHF)
          INDX = 0
          CALL RDUNIT( UNTS,VAR(1),INDX )
!
!---  Read initial condtions input from the input file  ---
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
              CALL RDINBS( T,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( T,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( T,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
            DO 202 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 202
              IC_NCG(1,N) = IC_NCG(1,N) + 100
              IC_NCG(2,N) = IC_NCG(2,N) + 100
  202       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( T,VAR(1),ADDER,IROCK,INDX )
            DO 204 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 204
              IF( IZ(N).EQ.IROCK ) IC_NCG(1,N) = IC_NCG(1,N) + 100
              IF( IZ(N).EQ.IROCK ) IC_NCG(2,N) = IC_NCG(2,N) + 100
  204       CONTINUE
          ELSE
            CALL RDINIS( T,VAR,ADDER,IDOM,INDX )
            DO 206 K = IDOM(5),IDOM(6)
            DO 206 J = IDOM(3),IDOM(4)
            DO 206 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 206
              IC_NCG(1,N) = IC_NCG(1,N) + 100
              IC_NCG(2,N) = IC_NCG(2,N) + 100
  206       CONTINUE
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          ADDER = -PATM
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( PL,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( PL,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( PL,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PL,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( PL,VAR,ADDER,IDOM,INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas pres').NE.0 ) THEN
          ADDER = -PATM
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( PG,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( PG,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( PG,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
            DO 208 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 208
              IC_NCG(1,N) = IC_NCG(1,N) + 10
              IC_NCG(2,N) = IC_NCG(2,N) + 10
  208       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PG,VAR(1),ADDER,IROCK,INDX )
            DO 210 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 210
              IF( IZ(N).EQ.IROCK ) IC_NCG(1,N) = IC_NCG(1,N) + 10
              IF( IZ(N).EQ.IROCK ) IC_NCG(2,N) = IC_NCG(2,N) + 10
  210       CONTINUE
          ELSE
            CALL RDINIS( PG,VAR,ADDER,IDOM,INDX )
            DO 212 K = IDOM(5),IDOM(6)
            DO 212 J = IDOM(3),IDOM(4)
            DO 212 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 212
              IC_NCG(1,N) = IC_NCG(1,N) + 10
              IC_NCG(2,N) = IC_NCG(2,N) + 10
  212       CONTINUE
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'sat').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( SH,ADDER,UNTS,INDX )
              DO 214 N = 1,NFLD
                IC_OPT(2,N) = 1
  214         CONTINUE
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( SH,ADDER,UNTS,INDX )
              DO 216 N = 1,NFLD
                IC_OPT(2,N) = 1
  216         CONTINUE
            ELSE
              CALL RDINFS( SH,VAR,ADDER,UNTS,INDX )
              DO 218 N = 1,NFLD
                IC_OPT(2,N) = 1
  218         CONTINUE
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( SH,VAR(1),ADDER,IROCK,INDX )
            DO 220 N = 1,NFLD
              IF( IZ(N).EQ.IROCK ) IC_OPT(2,N) = 1
  220       CONTINUE
          ELSE
            CALL RDINIS( SH,VAR,ADDER,IDOM,INDX )
            DO 222 K = IDOM(5),IDOM(6)
            DO 222 J = IDOM(3),IDOM(4)
            DO 222 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 222
              IC_OPT(2,N) = 1
  222       CONTINUE
          ENDIF
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
            DO 224 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 224
              ICBRN(N) = IVAR
  224       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( TMS,VAR(1),ADDER,IROCK,INDX )
            DO 226 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 226
              IF( IZ(N).EQ.IROCK ) ICBRN(N) = IVAR
  226       CONTINUE
          ELSE
            CALL RDINIS( TMS,VAR,ADDER,IDOM,INDX )
            DO 228 K = IDOM(5),IDOM(6)
            DO 228 J = IDOM(3),IDOM(4)
            DO 228 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 228
              ICBRN(N) = IVAR
  228       CONTINUE
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'co2').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'liquid').NE.0 .AND.
     &      INDEX(ADUM(1:),'press').NE.0 ) THEN
            IVAR = 9
          ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &      INDEX(ADUM(1:),'mole').NE.0 .AND.
     &      INDEX(ADUM(1:),'form').NE.0 ) THEN
            IVAR = 8
          ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &      INDEX(ADUM(1:),'mass').NE.0 .AND.
     &      INDEX(ADUM(1:),'form').NE.0 ) THEN
            IVAR = 7
          ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &      INDEX(ADUM(1:),'mole').NE.0 .AND.
     &      INDEX(ADUM(1:),'form').NE.0 ) THEN
            IVAR = 6
          ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &      INDEX(ADUM(1:),'mass').NE.0 .AND.
     &      INDEX(ADUM(1:),'form').NE.0 ) THEN
            IVAR = 5
          ELSEIF( INDEX(ADUM(1:),'partial').NE.0 .AND.
     &      INDEX(ADUM(1:),'press').NE.0 ) THEN
            IVAR = 4
          ELSEIF( INDEX(ADUM(1:),'mass').NE.0 .AND.
     &      INDEX(ADUM(1:),'frac').NE.0 ) THEN
            IVAR = 3
          ELSEIF( INDEX(ADUM(1:),'rel').NE.0 .AND.
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
            ELSEIF( IVAR.EQ.4 ) THEN
              IUNM = -1
              IUNKG = 1
              IUNS = -2
            ENDIF
            CALL RDINFS( PVA,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
            DO 306 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 306
              IF( MOD(IC_NCG(1,N),10).GT.0 ) THEN
                INDX = 7
                IMSG = N
                CHMSG = 'Mulitple Initial Condition Specification'
     &            // 'of CO2 at Node: '
                CALL WRMSGS( INDX )
              ENDIF
              IC_NCG(1,N) = IC_NCG(1,N) + IVAR
  306       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PVA,VAR(1),ADDER,IROCK,INDX )
            DO 308 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 308
              IF( IZ(N).EQ.IROCK ) THEN
                IF( MOD(IC_NCG(1,N),10).GT.0 ) THEN
                  INDX = 7
                  IMSG = N
                  CHMSG = 'Mulitple Initial Condition Specification'
     &              // 'of CO2 at Node: '
                  CALL WRMSGS( INDX )
                ENDIF
                IC_NCG(1,N) = IC_NCG(1,N) + IVAR
              ENDIF
  308       CONTINUE
          ELSE
            CALL RDINIS( PVA,VAR,ADDER,IDOM,INDX )
            DO 310 K = IDOM(5),IDOM(6)
            DO 310 J = IDOM(3),IDOM(4)
            DO 310 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 310
              IF( MOD(IC_NCG(1,N),10).GT.0 ) THEN
                INDX = 7
                IMSG = N
                CHMSG = 'Mulitple Initial Condition Specification'
     &            // 'of CO2 at Node: '
                CALL WRMSGS( INDX )
              ENDIF
              IC_NCG(1,N) = IC_NCG(1,N) + IVAR
  310       CONTINUE
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'ch4').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &      INDEX(ADUM(1:),'mole').NE.0 .AND.
     &      INDEX(ADUM(1:),'form').NE.0 ) THEN
            IVAR = 8
          ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &      INDEX(ADUM(1:),'mass').NE.0 .AND.
     &      INDEX(ADUM(1:),'form').NE.0 ) THEN
            IVAR = 7
          ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &      INDEX(ADUM(1:),'mole').NE.0 .AND.
     &      INDEX(ADUM(1:),'form').NE.0 ) THEN
            IVAR = 6
          ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &      INDEX(ADUM(1:),'mass').NE.0 .AND.
     &      INDEX(ADUM(1:),'form').NE.0 ) THEN
            IVAR = 5
          ELSEIF( INDEX(ADUM(1:),'partial').NE.0 .AND.
     &      INDEX(ADUM(1:),'press').NE.0 ) THEN
            IVAR = 4
          ELSEIF( INDEX(ADUM(1:),'mass').NE.0 .AND.
     &      INDEX(ADUM(1:),'frac').NE.0 ) THEN
            IVAR = 3
          ELSEIF( INDEX(ADUM(1:),'rel').NE.0 .AND.
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
            ELSEIF( IVAR.EQ.4 ) THEN
              IUNM = -1
              IUNKG = 1
              IUNS = -2
            ENDIF
            CALL RDINFS( PVO,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
            DO 326 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 326
              IF( MOD(IC_NCG(2,N),10).GT.0 ) THEN
                INDX = 7
                IMSG = N
                CHMSG = 'Mulitple Initial Condition Specification'
     &            // 'of CH4 at Node: '
                CALL WRMSGS( INDX )
              ENDIF
              IC_NCG(2,N) = IC_NCG(2,N) + IVAR
  326       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PVO,VAR(1),ADDER,IROCK,INDX )
            DO 328 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 328
              IF( IZ(N).EQ.IROCK ) THEN
                IF( MOD(IC_NCG(2,N),10).GT.0 ) THEN
                  INDX = 7
                  IMSG = N
                  CHMSG = 'Mulitple Initial Condition Specification'
     &              // 'of CH4 at Node: '
                  CALL WRMSGS( INDX )
                ENDIF
                IC_NCG(2,N) = IC_NCG(2,N) + IVAR
              ENDIF
  328       CONTINUE
          ELSE
            CALL RDINIS( PVO,VAR,ADDER,IDOM,INDX )
            DO 330 K = IDOM(5),IDOM(6)
            DO 330 J = IDOM(3),IDOM(4)
            DO 330 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 330
              IF( MOD(IC_NCG(2,N),10).GT.0 ) THEN
                INDX = 7
                IMSG = N
                CHMSG = 'Mulitple Initial Condition Specification'
     &            // 'of CH4 at Node: '
                CALL WRMSGS( INDX )
              ENDIF
              IC_NCG(2,N) = IC_NCG(2,N) + IVAR
  330       CONTINUE
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
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIC37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINPT37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Read input file cards.
!     Direct control to card reader subroutines.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/RDINPT37'
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
        CALL RDSP37
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
!---  Search input file for liquid-CO2 relative permeability card  ---
!
 1000 CONTINUE
 1009 READ(IRD,'(A)', END=1010) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 1009
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'liquid').NE.0 .AND.
     &  INDEX(CHDUM(2:),'co2').NE.0 .AND.
     &  INDEX(CHDUM(2:),'rel').NE.0 ) THEN
        CALL RDNRP
        REWIND(IRD)
        GOTO 1100
      ELSE
        GOTO 1000
      ENDIF
 1010 CONTINUE
      INDX = 4
      CHMSG = 'Missing Liquid-CO2 Relative Permeability Card'
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
        CALL RDTF37
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
        CALL RDTP37
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
        CALL RDST37
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
        CALL RDIC37
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
        CALL RDBC37
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
        CALL RDSR37
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
        CALL RDOU37
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
        CALL RDSF37
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
!---  End of RDINPT37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDOU37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Read input file for output information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/RDOU37'
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
        IF( INDEX(ADUM(1:),'aqueous pressure').NE.0 ) THEN
          IREF(NV) = 1
        ELSEIF( INDEX(ADUM(1:),'gas pressure').NE.0 ) THEN
          IREF(NV) = 2
        ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 ) THEN
          IREF(NV) = 4
        ELSEIF( INDEX(ADUM(1:),'phase condition').NE.0 ) THEN
          IREF(NV) = 5
        ELSEIF( INDEX(ADUM(1:),'aqueous gauge pressure').NE.0 ) THEN
          IREF(NV) = 6
        ELSEIF( INDEX(ADUM(1:),'gas gauge pressure').NE.0 ) THEN
          IREF(NV) = 7
        ELSEIF( INDEX(ADUM(1:),'apparent aqueous sat').NE.0 ) THEN
          IREF(NV) = 9
        ELSEIF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
          IREF(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'gas saturation').NE.0 ) THEN
          IREF(NV) = 12
        ELSEIF( INDEX(ADUM(1:),'liquid co2 sat').NE.0 .OR.
     &    INDEX(ADUM(1:),'liquid-co2 sat').NE.0 ) THEN
          IREF(NV) = 13
        ELSEIF( INDEX(ADUM(1:),'aqueous moisture cont').NE.0 ) THEN
          IREF(NV) = 15
        ELSEIF( INDEX(ADUM(1:),'co2 aqueous mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous co2 mole frac').NE.0 ) THEN
          CHREF(204) = 'XMLA'
          IREF(NV) = 204
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mole frac').NE.0 ) THEN
          CHREF(205) = 'XMLS'
          IREF(NV) = 205
        ELSEIF( INDEX(ADUM(1:),'effective trapped gas').NE.0 ) THEN
          IREF(NV) = 19
        ELSEIF( INDEX(ADUM(1:),'diffusive porosity').NE.0 ) THEN
          IREF(NV) = 20
        ELSEIF( INDEX(ADUM(1:),'h2o gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o mass frac').NE.0 ) THEN
          IREF(NV) = 21
        ELSEIF( INDEX(ADUM(1:),'co2 gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas co2 mass frac').NE.0 ) THEN
          IREF(NV) = 22
        ELSEIF( INDEX(ADUM(1:),'ch4 gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas ch4 mass frac').NE.0 ) THEN
          IREF(NV) = 23
        ELSEIF( INDEX(ADUM(1:),'h2o aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous h2o mass frac').NE.0 ) THEN
          IREF(NV) = 24
        ELSEIF( INDEX(ADUM(1:),'co2 aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous co2 mass frac').NE.0 ) THEN
          IREF(NV) = 25
        ELSEIF( INDEX(ADUM(1:),'ch4 aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous ch4 mass frac').NE.0 ) THEN
          IREF(NV) = 26
        ELSEIF( INDEX(ADUM(1:),'h2o liquid-co2 mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'liquid-co2 h2o mass frac').NE.0 ) THEN
          IREF(NV) = 163
        ELSEIF( INDEX(ADUM(1:),'co2 liquid-co2 mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'liquid-co2 co2 mass frac').NE.0 ) THEN
          IREF(NV) = 223
        ELSEIF( INDEX(ADUM(1:),'ch4 liquid-co2 mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'liquid-co2 ch4 mass frac').NE.0 ) THEN
          IREF(NV) = 166
        ELSEIF( INDEX(ADUM(1:),'aqueous hydraulic head').NE.0 ) THEN
          IREF(NV) = 27
        ELSEIF( INDEX(ADUM(1:),'gas hydraulic head').NE.0 ) THEN
          IREF(NV) = 28
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
        ELSEIF( INDEX(ADUM(1:),'liquid-co2 density').NE.0 ) THEN
          IREF(NV) = 36
        ELSEIF( INDEX(ADUM(1:),'total h2o mass').NE.0 ) THEN
          IREF(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'total co2 mass').NE.0 ) THEN
          IREF(NV) = 38
        ELSEIF( INDEX(ADUM(1:),'total ch4 mass').NE.0 ) THEN
          IREF(NV) = 39
        ELSEIF( INDEX(ADUM(1:),'h2o mass source int').NE.0 ) THEN
          IREF(NV) = 40
        ELSEIF( INDEX(ADUM(1:),'co2 mass source int').NE.0 ) THEN
          IREF(NV) = 41
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
        ELSEIF ( INDEX(ADUM(1:),'x salt flux').NE.0 ) THEN
          IREF(NV) = 64
        ELSEIF ( INDEX(ADUM(1:),'y salt flux').NE.0 ) THEN
          IREF(NV) = 65
        ELSEIF ( INDEX(ADUM(1:),'z salt flux').NE.0 ) THEN
          IREF(NV) = 66
        ELSEIF ( INDEX(ADUM(1:),'xnc salt flux').NE.0 ) THEN
          IREF(NV) = 67
        ELSEIF ( INDEX(ADUM(1:),'ync salt flux').NE.0 ) THEN
          IREF(NV) = 68
        ELSEIF ( INDEX(ADUM(1:),'znc salt flux').NE.0 ) THEN
          IREF(NV) = 69
        ELSEIF( INDEX(ADUM(1:),'h2o gas mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o mole').NE.0 ) THEN
          IREF(NV) = 70
        ELSEIF( INDEX(ADUM(1:),'co2 gas mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas co2 mole').NE.0 ) THEN
          IREF(NV) = 71
        ELSEIF( INDEX(ADUM(1:),'ch4 gas mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas ch4 mole').NE.0 ) THEN
          IREF(NV) = 72
        ELSEIF( INDEX(ADUM(1:),'h2o gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o conc').NE.0 ) THEN
          IREF(NV) = 73
        ELSEIF( INDEX(ADUM(1:),'co2 gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas co2 conc').NE.0 ) THEN
          IREF(NV) = 74
        ELSEIF( INDEX(ADUM(1:),'ch4 gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas ch4 conc').NE.0 ) THEN
          IREF(NV) = 75
        ELSEIF( INDEX(ADUM(1:),'h2o aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous h2o conc').NE.0 ) THEN
          IREF(NV) = 76
        ELSEIF( INDEX(ADUM(1:),'co2 aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous co2 conc').NE.0 ) THEN
          IREF(NV) = 77
        ELSEIF( INDEX(ADUM(1:),'ch4 aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous ch4 conc').NE.0 ) THEN
          IREF(NV) = 78
        ELSEIF( INDEX(ADUM(1:),'gas courant').NE.0 ) THEN
          ICRNT = 1
          IREF(NV) = 79
        ELSEIF( INDEX(ADUM(1:),'ice pressure').NE.0 ) THEN
          IREF(NV) = 80
        ELSEIF( INDEX(ADUM(1:),'ice saturation').NE.0 ) THEN
          IREF(NV) = 81
        ELSEIF( INDEX(ADUM(1:),'ice density').NE.0 ) THEN
          IREF(NV) = 82
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
        ELSEIF( INDEX(ADUM(1:),'osmotic pressure').NE.0 ) THEN
          IREF(NV) = 101
        ELSEIF( INDEX(ADUM(1:),'osmotic eff').NE.0 ) THEN
          IREF(NV) = 102
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mass').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mass').NE.0 ) THEN
          IREF(NV) = 110
        ELSEIF( INDEX(ADUM(1:),'gas-aqueous scaling').NE.0 ) THEN
          IREF(NV) = 131
        ELSEIF( INDEX(ADUM(1:),'h2o mass source rate').NE.0 ) THEN
          IREF(NV) = 140
        ELSEIF( INDEX(ADUM(1:),'co2 mass source rate').NE.0 ) THEN
          IREF(NV) = 141
        ELSEIF( INDEX(ADUM(1:),'salt mass source rate').NE.0 ) THEN
          IREF(NV) = 147
        ELSEIF( INDEX(ADUM(1:),'salt mass source int').NE.0 ) THEN
          IREF(NV) = 148
        ELSEIF( INDEX(ADUM(1:),'aqueous viscosity').NE.0 ) THEN
          IREF(NV) = 176
        ELSEIF( INDEX(ADUM(1:),'gas viscosity').NE.0 ) THEN
          IREF(NV) = 289
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    INDEX(ADUM(1:),'h2o').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 277
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 278
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    INDEX(ADUM(1:),'ch4').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 279
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    INDEX(ADUM(1:),'ch4').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
          IREF(NV) = 280
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    INDEX(ADUM(1:),'ch4').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IREF(NV) = 281
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    INDEX(ADUM(1:),'h2o').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 282
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 283
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    INDEX(ADUM(1:),'ch4').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 284
        ELSEIF( INDEX(ADUM(1:),'integrated h2o mass').NE.0 ) THEN
          IREF(NV) = 191
        ELSEIF( INDEX(ADUM(1:),'integrated co2 mass').NE.0 ) THEN
          IREF(NV) = 192
        ELSEIF( INDEX(ADUM(1:),'integrated ch4 mass').NE.0 ) THEN
          IREF(NV) = 193
        ELSEIF( INDEX(ADUM(1:),'integrated h2o aqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated aqueous h2o').NE.0 ) THEN
          IREF(NV) = 194
        ELSEIF( INDEX(ADUM(1:),'integrated co2 aqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated aqueous co2').NE.0 ) THEN
          IREF(NV) = 195
        ELSEIF( INDEX(ADUM(1:),'integrated ch4 aqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated aqueous ch4').NE.0 ) THEN
          IREF(NV) = 259
        ELSEIF( INDEX(ADUM(1:),'integrated h2o gas').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated gas h2o').NE.0 ) THEN
          IREF(NV) = 197
        ELSEIF( INDEX(ADUM(1:),'integrated co2 gas').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated gas co2').NE.0 ) THEN
          IREF(NV) = 198
        ELSEIF( INDEX(ADUM(1:),'integrated ch4 gas').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated gas ch4').NE.0 ) THEN
          IREF(NV) = 260
        ELSEIF( INDEX(ADUM(1:),'h2o hydrate mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate h2o mass frac').NE.0 ) THEN
          IREF(NV) = 250
        ELSEIF( INDEX(ADUM(1:),'co2 hydrate mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate co2 mass frac').NE.0 ) THEN
          IREF(NV) = 251
        ELSEIF( INDEX(ADUM(1:),'ch4 hydrate mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate ch4 mass frac').NE.0 ) THEN
          IREF(NV) = 252
        ELSEIF( INDEX(ADUM(1:),'hydrate dens').NE.0 ) THEN
          IREF(NV) = 253
        ELSEIF( INDEX(ADUM(1:),'hydrate satur').NE.0 ) THEN
          IREF(NV) = 254
        ELSEIF( INDEX(ADUM(1:),'hydrate press').NE.0 ) THEN
          IREF(NV) = 255
        ELSEIF( INDEX(ADUM(1:),'integrated h2o hydrate').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated hydrate h2o').NE.0 ) THEN
          IREF(NV) = 256
        ELSEIF( INDEX(ADUM(1:),'integrated co2 hydrate').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated hydrate co2').NE.0 ) THEN
          IREF(NV) = 257
        ELSEIF( INDEX(ADUM(1:),'integrated ch4 hydrate').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated hydrate ch4').NE.0 ) THEN
          IREF(NV) = 258
        ELSEIF( INDEX(ADUM(1:),'salt saturation').NE.0 ) THEN
          IREF(NV) = 264
        ELSEIF( INDEX(ADUM(1:),'h2o hydrate mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate h2o mole frac').NE.0 ) THEN
          IREF(NV) = 265
        ELSEIF( INDEX(ADUM(1:),'co2 hydrate mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate co2 mole frac').NE.0 ) THEN
          IREF(NV) = 266
        ELSEIF( INDEX(ADUM(1:),'ch4 hydrate mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate ch4 mole frac').NE.0 ) THEN
          IREF(NV) = 267
        ELSEIF( INDEX(ADUM(1:),'solute volumetric conc').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 1
        ELSEIF( INDEX(ADUM(1:),'solute aqueous conc').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 2
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
        ELSEIF( INDEX(ADUM(1:),'solute source').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 11
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
      IF( NPRTM.GT.LPTM ) THEN
        INDX = 5
        CHMSG = 'Number of Output Times > Parameter LPTM'
        CALL WRMSGS( INDX )
      ENDIF
      WRITE(IWR,'(/,A)') ' Plot File Output Times:'
      DO 300 N = 1, NPRTM
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Plot File Output Time'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PRTM(N))
        VARB = 'Plot File Output Time Units'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,1PE11.4,1X,A)') PRTM(N),UNTS(1:NCH)
        INDX = 0
        IUNS = 1
        CALL RDUNIT(UNTS,PRTM(N),INDX)
        TMPR = MIN( TMPR,PRTM(N) )
 300  CONTINUE
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
        ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 ) THEN
          IPLOT(NV) = 4
        ELSEIF( INDEX(ADUM(1:),'phase condition').NE.0 ) THEN
          IPLOT(NV) = 5
        ELSEIF( INDEX(ADUM(1:),'aqueous gauge pressure').NE.0 ) THEN
          IPLOT(NV) = 6
        ELSEIF( INDEX(ADUM(1:),'gas gauge pressure').NE.0 ) THEN
          IPLOT(NV) = 7
        ELSEIF( INDEX(ADUM(1:),'apparent aqueous sat').NE.0 ) THEN
          IPLOT(NV) = 9
        ELSEIF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
          IPLOT(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'gas saturation').NE.0 ) THEN
          IPLOT(NV) = 12
        ELSEIF( INDEX(ADUM(1:),'liquid co2 sat').NE.0 .OR.
     &    INDEX(ADUM(1:),'liquid-co2 sat').NE.0 ) THEN
          IPLOT(NV) = 13
        ELSEIF( INDEX(ADUM(1:),'aqueous moisture cont').NE.0 ) THEN
          IPLOT(NV) = 15
        ELSEIF( INDEX(ADUM(1:),'co2 aqueous mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous co2 mole frac').NE.0 ) THEN
          IPLOT(NV) = 204
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mole frac').NE.0 ) THEN
          IPLOT(NV) = 205
        ELSEIF( INDEX(ADUM(1:),'effective trapped gas').NE.0 ) THEN
          IPLOT(NV) = 19
        ELSEIF( INDEX(ADUM(1:),'diffusive porosity').NE.0 ) THEN
          IPLOT(NV) = 20
        ELSEIF( INDEX(ADUM(1:),'h2o gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o mass frac').NE.0 ) THEN
          IPLOT(NV) = 21
        ELSEIF( INDEX(ADUM(1:),'co2 gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas co2 mass frac').NE.0 ) THEN
          IPLOT(NV) = 22
        ELSEIF( INDEX(ADUM(1:),'ch4 gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas ch4 mass frac').NE.0 ) THEN
          IPLOT(NV) = 23
        ELSEIF( INDEX(ADUM(1:),'h2o aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous h2o mass frac').NE.0 ) THEN
          IPLOT(NV) = 24
        ELSEIF( INDEX(ADUM(1:),'co2 aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous co2 mass frac').NE.0 ) THEN
          IPLOT(NV) = 25
        ELSEIF( INDEX(ADUM(1:),'ch4 aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous ch4 mass frac').NE.0 ) THEN
          IPLOT(NV) = 26
        ELSEIF( INDEX(ADUM(1:),'h2o liquid-co2 mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'liquid-co2 h2o mass frac').NE.0 ) THEN
          IPLOT(NV) = 163
        ELSEIF( INDEX(ADUM(1:),'co2 liquid-co2 mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'liquid-co2 co2 mass frac').NE.0 ) THEN
          IPLOT(NV) = 223
        ELSEIF( INDEX(ADUM(1:),'ch4 liquid-co2 mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'liquid-co2 ch4 mass frac').NE.0 ) THEN
          IPLOT(NV) = 166
        ELSEIF( INDEX(ADUM(1:),'aqueous hydraulic head').NE.0 ) THEN
          IPLOT(NV) = 27
        ELSEIF( INDEX(ADUM(1:),'gas hydraulic head').NE.0 ) THEN
          IPLOT(NV) = 28
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
        ELSEIF( INDEX(ADUM(1:),'total h2o mass').NE.0 ) THEN
          IPLOT(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'total co2 mass').NE.0 ) THEN
          IPLOT(NV) = 38
        ELSEIF( INDEX(ADUM(1:),'total ch4 mass').NE.0 ) THEN
          IPLOT(NV) = 39
        ELSEIF( INDEX(ADUM(1:),'h2o mass source int').NE.0 ) THEN
          IPLOT(NV) = 40
        ELSEIF( INDEX(ADUM(1:),'co2 mass source int').NE.0 ) THEN
          IPLOT(NV) = 41
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
        ELSEIF ( INDEX(ADUM(1:),'x salt flux').NE.0 ) THEN
          IPLOT(NV) = 64
        ELSEIF ( INDEX(ADUM(1:),'y salt flux').NE.0 ) THEN
          IPLOT(NV) = 65
        ELSEIF ( INDEX(ADUM(1:),'z salt flux').NE.0 ) THEN
          IPLOT(NV) = 66
        ELSEIF ( INDEX(ADUM(1:),'xnc salt flux').NE.0 ) THEN
          IPLOT(NV) = 67
        ELSEIF ( INDEX(ADUM(1:),'ync salt flux').NE.0 ) THEN
          IPLOT(NV) = 68
        ELSEIF ( INDEX(ADUM(1:),'znc salt flux').NE.0 ) THEN
          IPLOT(NV) = 69
        ELSEIF( INDEX(ADUM(1:),'h2o gas mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o mole').NE.0 ) THEN
          IPLOT(NV) = 70
        ELSEIF( INDEX(ADUM(1:),'co2 gas mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas co2 mole').NE.0 ) THEN
          IPLOT(NV) = 71
        ELSEIF( INDEX(ADUM(1:),'h2o gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o conc').NE.0 ) THEN
          IPLOT(NV) = 73
        ELSEIF( INDEX(ADUM(1:),'co2 gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas co2 conc').NE.0 ) THEN
          IPLOT(NV) = 74
        ELSEIF( INDEX(ADUM(1:),'h2o aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous h2o conc').NE.0 ) THEN
          IPLOT(NV) = 76
        ELSEIF( INDEX(ADUM(1:),'co2 aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous co2 conc').NE.0 ) THEN
          IPLOT(NV) = 77
        ELSEIF( INDEX(ADUM(1:),'gas courant').NE.0 ) THEN
          ICRNT = 1
          IPLOT(NV) = 79
        ELSEIF( INDEX(ADUM(1:),'ice pressure').NE.0 ) THEN
          IPLOT(NV) = 80
        ELSEIF( INDEX(ADUM(1:),'ice saturation').NE.0 ) THEN
          IPLOT(NV) = 81
        ELSEIF( INDEX(ADUM(1:),'ice density').NE.0 ) THEN
          IPLOT(NV) = 82
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
        ELSEIF( INDEX(ADUM(1:),'osmotic press').NE.0 ) THEN
          IPLOT(NV) = 101
        ELSEIF( INDEX(ADUM(1:),'osmotic eff').NE.0 ) THEN
          IPLOT(NV) = 102
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mass').NE.0 ) THEN
          IPLOT(NV) = 110
        ELSEIF( INDEX(ADUM(1:),'gas-aqueous scaling').NE.0 ) THEN
          IPLOT(NV) = 131
        ELSEIF( INDEX(ADUM(1:),'h2o mass source rate').NE.0 ) THEN
          IPLOT(NV) = 140
        ELSEIF( INDEX(ADUM(1:),'co2 mass source rate').NE.0 ) THEN
          IPLOT(NV) = 141
        ELSEIF( INDEX(ADUM(1:),'salt mass source rate').NE.0 ) THEN
          IPLOT(NV) = 147
        ELSEIF( INDEX(ADUM(1:),'salt mass source int').NE.0 ) THEN
          IPLOT(NV) = 148
        ELSEIF( INDEX(ADUM(1:),'aqueous viscosity').NE.0 ) THEN
          IPLOT(NV) = 176
        ELSEIF( INDEX(ADUM(1:),'gas viscosity').NE.0 ) THEN
          IPLOT(NV) = 289
        ELSEIF( INDEX(ADUM(1:),'h2o hydrate mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate h2o mass frac').NE.0 ) THEN
          IPLOT(NV) = 250
        ELSEIF( INDEX(ADUM(1:),'co2 hydrate mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate co2 mass frac').NE.0 ) THEN
          IPLOT(NV) = 251
        ELSEIF( INDEX(ADUM(1:),'ch4 hydrate mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate ch4 mass frac').NE.0 ) THEN
          IPLOT(NV) = 252
        ELSEIF( INDEX(ADUM(1:),'hydrate dens').NE.0 ) THEN
          IPLOT(NV) = 253
        ELSEIF( INDEX(ADUM(1:),'hydrate satur').NE.0 ) THEN
          IPLOT(NV) = 254
          IF( INDEX(ADUM(1:),'equil').NE.0 ) IPLOT(NV) = 256
        ELSEIF( INDEX(ADUM(1:),'hydrate press').NE.0 ) THEN
          IPLOT(NV) = 255
        ELSEIF( INDEX(ADUM(1:),'salt saturation').NE.0 ) THEN
          IPLOT(NV) = 264
        ELSEIF( INDEX(ADUM(1:),'h2o hydrate mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate h2o mole frac').NE.0 ) THEN
          IPLOT(NV) = 265
        ELSEIF( INDEX(ADUM(1:),'co2 hydrate mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate co2 mole frac').NE.0 ) THEN
          IPLOT(NV) = 266
        ELSEIF( INDEX(ADUM(1:),'ch4 hydrate mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'hydrate ch4 mole frac').NE.0 ) THEN
          IPLOT(NV) = 267
        ELSEIF( INDEX(ADUM(1:),'solute volumetric conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 1
        ELSEIF( INDEX(ADUM(1:),'solute aqueous conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 2
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
        ELSEIF( INDEX(ADUM(1:),'solute source').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 11
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDOU37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSF37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Read input file surface flux information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/RDSF37'
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
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0) THEN
          IF( INDEX(ADUM(1:),'volum').NE.0) THEN
            ISFT(NS) = 2
            WRITE(IWR,'(A)') 'Aqueous Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ELSEIF( INDEX(ADUM(1:),'co2').NE.0) THEN
            ISFT(NS) = 29
            WRITE(IWR,'(A)') 'Aqueous CO2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'ch4').NE.0) THEN
            ISFT(NS) = 12
            WRITE(IWR,'(A)') 'Aqueous CH4 Mass Flux Surface'
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
        ELSEIF( INDEX(ADUM(1:),'liquid-co2').NE.0 .OR.
     &    INDEX(ADUM(1:),'liquid co2').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'volum').NE.0) THEN
            ISFT(NS) = 4
            WRITE(IWR,'(A)') 'Liquid-CO2 Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ELSEIF( INDEX(ADUM(1:),'co2').NE.0) THEN
            ISFT(NS) = 55
            WRITE(IWR,'(A)') 'Liquid-CO2 CO2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'ch4').NE.0) THEN
            ISFT(NS) = 56
            WRITE(IWR,'(A)') 'Liquid-CO2 CH4 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'water').NE.0 .OR.
     &      INDEX(ADUM(1:),'h2o').NE.0 ) THEN
            ISFT(NS) = 57
            WRITE(IWR,'(A)') 'Liquid-CO2 Water Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            ISFT(NS) = 7
            WRITE(IWR,'(A)') 'Liquid-CO2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'volum').NE.0) THEN
            ISFT(NS) = 3
            WRITE(IWR,'(A)') 'Gas Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ELSEIF( INDEX(ADUM(1:),'co2').NE.0) THEN
            ISFT(NS) = 28
            WRITE(IWR,'(A)') 'Gas CO2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'ch4').NE.0) THEN
            ISFT(NS) = 11
            WRITE(IWR,'(A)') 'Gas CH4 Mass Flux Surface'
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
        ELSEIF( INDEX(ADUM(1:),'co2').NE.0) THEN
          ISFT(NS) = 30
          WRITE(IWR,'(A)') 'CO2-Mass Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'ch4').NE.0) THEN
          ISFT(NS) = 13
          WRITE(IWR,'(A)') 'CH4-Mass Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'water').NE.0) THEN
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
!---  End of RDSF37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSP37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Read input file for rock/soil saturation function information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/RDSP37'
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
      VARB ='Liquid CO2-Aqueous Interfacial Tension'
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
      VARB ='Gas-Liquid CO2 Interfacial Tension'
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
     &    'Liquid CO2-Aqueous Surface Tension'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Check for zero gas-NAPL interfacial tension  ---
!
      IF( SCHR(7,1).LT.EPSL ) THEN
        INDX = 4
        CHMSG = 'Zero Gas-Liquid CO2 Interfacial Tension'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Write spreading coefficient  ---
!
      C_SPRD = MIN( (SCHR(5,1)-SCHR(6,1)-SCHR(7,1)),0.D+0 )
      WRITE(IWR,'(2X,A,1PE11.4)') 
     &  'Liquid CO2 Spreading Coefficient: ',C_SPRD
      IF( C_SPRD.LT.-EPSL ) THEN
        CA_NL = ACOS( ((SCHR(6,1)**2) - (SCHR(7,1)**2) + (SCHR(5,1)**2))
     &    /(2.D+0*SCHR(5,1)*SCHR(6,1)) )
        CA_GN = ACOS( ((SCHR(7,1)**2) - (SCHR(6,1)**2) + (SCHR(5,1)**2))
     &    /(2.D+0*SCHR(5,1)*SCHR(7,1)) )
!
!---    Write gas-NAPL and NAPL-aqueous contact angles  ---
!

        WRITE(IWR,'(2X,A,1PE11.4)') 'Gas-Liquid CO2 Contact ' //
     &    'Angle, deg: ',(CA_GN*1.8D+2/GPI)
        WRITE(IWR,'(2X,A,1PE11.4)') 'Liquid CO2-Aqueous ' //
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
          GOTO 20
        ENDIF
        INDX = 2
        CHMSG = 'Unrecognized Rock/Soil Type: ' // RDUM(1:NCH)
        CALL WRMSGS( INDX )
        GOTO 20
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
!---  Read saturation/capillary pressure function  ---
!
        VARB = 'Saturation Function Type: '
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        ISMX = 0
        IF( INDEX(ADUM(1:),'extended').NE.0 ) ISMX = 1
        IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &      INDEX(ADUM(1:),'barrier').NE.0 ) THEN
            ISCHRX = 11
          ELSE
            ISCHRX = 1
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &      INDEX(ADUM(1:),'corey').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &      INDEX(ADUM(1:),'barrier').NE.0 ) THEN
            ISCHRX = 12
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
  230     CONTINUE
        ELSE
          ISCHR(IROCK) = ISCHRX
          ISM(IROCK) = ISMX
        ENDIF
!
!---    van Genuchten Function  ---
!
        IF( ISCHRX.EQ.1 .OR. ISCHRX.EQ.11 ) THEN
          IF( ISCHRX.EQ.1 ) THEN
            WRITE(IWR,'(A)') 'van Genuchten Function'
          ELSEIF( ISCHRX.EQ.11 ) THEN
           WRITE(IWR,'(A)')'van Genuchten  Function w/ Hydrate Barrier'
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
            INDX = 16
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
     &        UNTS(1:NCH),': ',SCHR(8,IROCK)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHR(8,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(8,IROCK),', N/m)'
          ENDIF
!
!---      Hydrate barrier model  ---
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
     &          SCHR(10,IROCK)
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
     &          SCHR(11,IROCK)
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
     &          SCHR(12,IROCK)
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
     &          SCHR(13,IROCK)
            ENDIF
          ENDIF
!
!---    Brooks and Corey Function  ---
!
        ELSEIF( ISCHRX.EQ.2 .OR. ISCHRX.EQ.12 ) THEN
          IF( ISCHRX.EQ.2 ) THEN
            WRITE(IWR,'(A)') 'Brooks and Corey Function'
          ELSEIF( ISCHRX.EQ.12 ) THEN
           WRITE(IWR,'(A)')'Brooks and Corey Function ' // 
     &       'w/ Hydrate Barrier'
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
     &        ': ',SCHR(2,IROCK)
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
            INDX = 16
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
     &        UNTS(1:NCH),': ',SCHR(8,IROCK)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHR(8,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(8,IROCK),', N/m)'
          ENDIF
!
!---      Hydrate barrier model  ---
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
     &          SCHR(10,IROCK)
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
     &          SCHR(11,IROCK)
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
     &          SCHR(12,IROCK)
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
     &          SCHR(13,IROCK)
            ENDIF
          ENDIF
        ENDIF
!
!---    Bound-water aqueous saturation  ---
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
     &      SCHR(15,IROCK)
        ENDIF
!
!---    Loop over remaining rock/soils within scaling group  ---
!
        IF( ISLC(19).EQ.1 .AND. IROCK.LT.NROCK ) THEN
          DO 490 M = IROCK+1,NROCK
            IF( ISCALE(M).EQ.ISGRP ) THEN
              N = N+1
              ISM(M) = ISM(IROCK)
              ISCHR(M) = ISCHR(IROCK)
              DO 480 L = 1,LSCHR
                SCHR(L,M) = SCHR(L,IROCK)
  480         CONTINUE
              DO 482 L = 1,2
                ISLTBL(L,M) = ISLTBL(L,IROCK)
  482         CONTINUE
            ENDIF
  490     CONTINUE
        ENDIF
!
!---  Read next rock/soil type or scaling group  ---
!
        IF( N.LT.NROCK ) WRITE(IWR,'(/)')
        GOTO 20
 500  CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSP37 group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSR37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Read input file for source information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      CHARACTER*64 ADUM,BDUM,CDUM,UNTS
      CHARACTER*512 CHDUM
      REAL*8 VAR(LSTM,8+LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSR37'
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
        ELSEIF( INDEX(ADUM(1:),'fluid').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'ic').NE.0 .OR.
     &    INDEX(ADUM(1:),'initial').NE.0 ) ) THEN
          WRITE(IWR,'(2X,A)') 'Fluid-Mass Source'
          ISRTX = 4
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
              GOTO 40
            ENDIF
   30     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute Name: '//ADUM
          CALL WRMSGS( INDX )
   40     CONTINUE
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Source Type: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Read source domain indices  ---
!
        VARB = 'Source Domain Index'
        I1X = ISTART
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
          CHMSG = 'Invalid Source Domain: ' // CHDUM(I1X:ICX)
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
!---  Check for z-direction well sources  ---
!
        IF( ISRTX.EQ.4 ) THEN
          IF( ((I2X-I1X).GT.1) .OR. ((J2X-J1X).GT.1) ) THEN
            INDX = 4
            CHMSG = 'Invalid Well-Domain Indices'
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  Read number of source times  ---
!
        VARB = 'Number of Source Times'
        CALL RDINT(ISTART,ICOMMA,CHDUM,ISRM(NS))
        IF( ISRM(NS).GT.LSTM ) THEN
          INDX = 5
          CHMSG = 'Number of Source Times > Parameter LSTM'
          CALL WRMSGS( INDX )
        ENDIF
        SRTMO = -SMALL
        DO 100 NTM = 1,ISRM(NS)
          DO 60 M = 1,8
            VAR(NTM,M) = 0.D+0
   60     CONTINUE
!
!---      Read and write source values and units  ---
!
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
          VARB = 'Source: Time'
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
!---      Initial Condition Fluid Mass Sink  ---
!
          ELSEIF( ISRTX.EQ.4 ) THEN
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
!---  End of RDSR37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDST37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Reads the salt transport card dispersivities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/RDST37'
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDST37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTF37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Reads solute/fluid interaction card for diffusion and partition
!     coefficients, and internodal diffusion term averaging scheme for
!     single phase (aqueous) solute transport equation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/RDTF37'
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
!---  End of RDTF37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTP37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Reads the solute/porous media interaction card for the
!     dispersivities, half-lives, and partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/RDTP37'
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
!---  End of RDTP37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Compute the maximum relative residuals
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE GRID
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
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
      CHARACTER*64 PH_CND(6)
      INTEGER IRSDX(LFD)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Aqueous-Gas w/o Hydrate',
     &  'Aqueous-Gas w/ Hydrate',
     &  'Aqueous-Low Liquid CO2-Gas w/o Hydrate',
     &  'Aqueous-Low Liquid CO2-Gas w Hydrate',
     &  'Aqueous-Liquid CO2-Gas w/o Hydrate',
     &  'Aqueous-Liquid CO2-Gas w Hydrate'/
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL37'
!
!---  Zero maximum residuals  ---
!
      DO 100 M = 1,ISVC
        RSD(M) = 0.D+0
        NSD(M) = 0
  100 CONTINUE
!
!---  Loop over all active nodes  ---
!
      DO 200 N = 1,NFLD
        IRSDX(N) = 0
        IF( IXP(N).EQ.0 ) GOTO 200
        IZN = IZ(N)
        N_DB = N
        NMD = IXP(N)
!
!---    Isothermal option  ---
!
        IF( ISLC(30).EQ.0 ) MPT = IM(IEQT,NMD)
        MPW = IM(IEQW,NMD)
!
!---    Iso-CO2 option  ---
!
        IF( ISLC(45).EQ.0 ) MPA = IM(IEQA,NMD)
        MPO = IM(IEQDO,NMD)
        DPO = BLU(MPO)
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) MPS = IM(IEQS,NMD)
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
          RSDX = MIN( (ABS(BLU(MPS))/XLSMX),
     &      ABS(RSDL(IEQS,N)/(ACP+SMALL)) )
          RSDX = RSDX*1.D-1
          IF( RSDX.GT.RSD(IEQS) ) THEN
            RSD(IEQS) = RSDX
            NSD(IEQS) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
        ENDIF
!
!---    Aqueous or aqueous-gas, no hydrate, no liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - CO2 vapor pressure
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
!          IF( T(2,N).GT.43.D+0 .AND. T(2,N).LT.47.D+0 ) GOTO 110
!
!---      Energy equation, isothermal option  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            RSDX = ABS(BLU(MPT))/TABS
            IF( RSDX.GT.RSD(IEQT) ) THEN
              RSD(IEQT) = RSDX
              NSD(IEQT) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Water residual  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N) +
     &      RHOH(2,N)*SH(2,N)*XHW(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPW))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      CO2 residual, iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            PX = MAX( PL(2,N),0.D+0 ) + PATM
            IF( PVA(2,N).GT.(1.D-6*PX) ) THEN
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLA(2,N) +
     &          RHOH(2,N)*SH(2,N)*XHA(2,N))*DTI*VOL(N)
              INCG = 1
              RSDX = MIN( ABS(BLU(MPA))/PX,
     &          ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD(IEQA) ) THEN
                RSD(IEQA) = RSDX
                NSD(IEQA) = N
              ENDIF
              IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
            ENDIF
          ENDIF
!
!---      CH4 residual  ---
!
          PX = MAX( PL(2,N),0.D+0 ) + PATM
          IF( PVO(2,N).GT.(1.D-6*PX) ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGO(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLO(2,N) + 
     &        RHOH(2,N)*SH(2,N)*XHO(2,N))*DTI*VOL(N)
            INCG = 2
            RSDX = MIN( ABS(DPO)/PX,
     &        ABS(RSDL(IEQDO,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQDO) ) THEN
              RSD(IEQDO) = RSDX
              NSD(IEQDO) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
  110     CONTINUE
!
!---    Aqueous or aqueous-gas, hydrate, no liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - Gas CH4 mole fraction of formers
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
!          IF( RKL(1,2,N).GT.EPSL ) THEN
!
!---        Energy equation, isothermal option  ---
!
            IF( ISLC(30).EQ.0 ) THEN
              RSDX = ABS(BLU(MPT))/TABS
              IF( RSDX.GT.RSD(IEQT) ) THEN
                RSD(IEQT) = RSDX
                NSD(IEQT) = N
              ENDIF
              IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
            ENDIF
!
!---        CO2 residual, iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLA(2,N) +
     &          RHOH(2,N)*SH(2,N)*XHA(2,N))*DTI*VOL(N)
              INCG = 1
              RSDX = MIN( ABS(BLU(MPA)),
     &          ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD(IEQA) ) THEN
                RSD(IEQA) = RSDX
                NSD(IEQA) = N
              ENDIF
              IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
            ENDIF
!
!---        CH4 residual  ---
!
            IF( 1.D+0-SH(2,N).GT.EPSL ) THEN
              IF( SL(2,N).LT.SCHR(15,IZN) ) DPO = MIN( DPO,0.D+0 )
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGO(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &          RHOH(2,N)*SH(2,N)*XHO(2,N))*DTI*VOL(N)
              INCG = 1
              RSDX = MIN( ABS(DPO),
     &          ABS(RSDL(IEQDO,N)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD(IEQDO) ) THEN
                RSD(IEQDO) = RSDX
                NSD(IEQDO) = N
              ENDIF
              IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
            ENDIF
!
!---        Water residual  ---
!
            ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &        RHOG(2,N)*SG(2,N)*XGW(2,N) + 
     &        RHOH(2,N)*SH(2,N)*XHW(2,N))*DTI*VOL(N)
            RSDX = MIN( ABS(BLU(MPW))/(ABS(PL(2,N))+PATM),
     &        ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQW) ) THEN
              RSD(IEQW) = RSDX
              NSD(IEQW) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!          ENDIF
!
!---    Aqueous or aqueous-gas, no hydrate, low liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 saturation
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
!
!---      Energy equation, isothermal option  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            RSDX = ABS(BLU(MPT))/TABS
            IF( RSDX.GT.RSD(IEQT) ) THEN
              RSD(IEQT) = RSDX
              NSD(IEQT) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Water residual  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N) +
     &      RHOH(2,N)*SH(2,N)*XHW(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPW))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      CO2 residual, iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLA(2,N) +
     &        RHOH(2,N)*SH(2,N)*XHA(2,N))*DTI*VOL(N)
            INCG = 1
            RSDX = MIN( ABS(BLU(MPA)),
     &        ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQA) ) THEN
              RSD(IEQA) = RSDX
              NSD(IEQA) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      CH4 residual  ---
!
          PX = MIN( PL(2,N),0.D+0 ) + PATM
          IF( PVO(2,N).GT.(1.D-9*PX) ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGO(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLO(2,N) + 
     &        RHOH(2,N)*SH(2,N)*XHO(2,N))*DTI*VOL(N)
            INCG = 2
            RSDX = MIN( ABS(DPO)/PX,
     &        ABS(RSDL(IEQDO,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQDO) ) THEN
              RSD(IEQDO) = RSDX
              NSD(IEQDO) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---    Aqueous or aqueous-gas, hydrate, low liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 saturation
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
!
!---      Energy equation, isothermal option  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            RSDX = ABS(BLU(MPT))/TABS
            IF( RSDX.GT.RSD(IEQT) ) THEN
              RSD(IEQT) = RSDX
              NSD(IEQT) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Water residual  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N) +
     &      RHOH(2,N)*SH(2,N)*XHW(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPW))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      CO2 residual, iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLA(2,N) +
     &        RHOH(2,N)*SH(2,N)*XHA(2,N))*DTI*VOL(N)
            INCG = 1
            RSDX = MIN( ABS(BLU(MPA)),
     &        ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQA) ) THEN
              RSD(IEQA) = RSDX
              NSD(IEQA) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      CH4 residual  ---
!
          IF( 1.D+0-SH(2,N).GT.EPSL ) THEN
            IF( SL(2,N).LT.SCHR(15,IZN) ) DPO = MIN( DPO,0.D+0 )
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGO(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &        RHOH(2,N)*SH(2,N)*XHO(2,N))*DTI*VOL(N)
            INCG = 1
            RSDX = MIN( ABS(DPO),
     &        ABS(RSDL(IEQDO,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQDO) ) THEN
              RSD(IEQDO) = RSDX
              NSD(IEQDO) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---    Aqueous or aqueous-gas, no hydrate, liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 pressure
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
!
!---      Energy equation, isothermal option  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            RSDX = ABS(BLU(MPT))/TABS
            IF( RSDX.GT.RSD(IEQT) ) THEN
              RSD(IEQT) = RSDX
              NSD(IEQT) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Water residual  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N) +
     &      RHOH(2,N)*SH(2,N)*XHW(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPW))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      CO2 residual, iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLA(2,N) +
     &        RHOH(2,N)*SH(2,N)*XHA(2,N))*DTI*VOL(N)
            INCG = 1
            RSDX = MIN( ABS(BLU(MPA))/(ABS(PN(2,N))+PATM),
     &        ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQA) ) THEN
              RSD(IEQA) = RSDX
              NSD(IEQA) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      CH4 residual  ---
!
          PX = MIN( PL(2,N),0.D+0 ) + PATM
          IF( PVO(2,N).GT.(1.D-9*PX) ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGO(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLO(2,N) + 
     &        RHOH(2,N)*SH(2,N)*XHO(2,N))*DTI*VOL(N)
            INCG = 2
            RSDX = MIN( ABS(DPO)/PX,
     &        ABS(RSDL(IEQDO,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQDO) ) THEN
              RSD(IEQDO) = RSDX
              NSD(IEQDO) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---    Aqueous or aqueous-gas, hydrate, liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 pressure
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
!
!---      Energy equation, isothermal option  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            RSDX = ABS(BLU(MPT))/TABS
            IF( RSDX.GT.RSD(IEQT) ) THEN
              RSD(IEQT) = RSDX
              NSD(IEQT) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Water residual  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N) +
     &      RHOH(2,N)*SH(2,N)*XHW(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPW))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      CO2 residual, iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLA(2,N) +
     &        RHOH(2,N)*SH(2,N)*XHA(2,N))*DTI*VOL(N)
            INCG = 1
            RSDX = MIN( ABS(BLU(MPA))/(ABS(PN(2,N))+PATM),
     &        ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQA) ) THEN
              RSD(IEQA) = RSDX
              NSD(IEQA) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      CH4 residual  ---
!
          IF( 1.D+0-SH(2,N).GT.EPSL ) THEN
            IF( SL(2,N).LT.SCHR(15,IZN) ) DPO = MIN( DPO,0.D+0 )
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGO(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &        RHOH(2,N)*SH(2,N)*XHO(2,N))*DTI*VOL(N)
            INCG = 1
            RSDX = MIN( ABS(DPO),
     &        ABS(RSDL(IEQDO,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQDO) ) THEN
              RSD(IEQDO) = RSDX
              NSD(IEQDO) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  Assign a convergence index  ---
!
!      RSDHX = 1.D+0
!      RSDAX = 0.D+0
!      NC = 0
      DO 300 M = 1,ISVC
        IF( RSD(M).GT.RSDMX ) ICNV = 2
!        IF( RSD(M).GT.EPSL ) THEN
!          RSDHX = RSDHX*RSD(M)
!          RSDAX = RSDAX + RSD(M)
!          NC = NC + 1
!        ENDIF
  300 CONTINUE
!      RSDHX = RSDHX**(1.D+0/REAL(NC))
!      RSDAX = RSDAX/REAL(NC)
!      IF( RSDHX.LT.RSDMX ) ICNV = 3
      IF( ICNV.EQ.2 .AND. NITER.GE.NRIMX ) ICNV = 1
!
!---  Unconverged solution Newton-Raphson iteration limit exceeded  ---
!
      IF( ICNV.EQ.1 ) THEN
        WRITE(ISC,'(10X,A)') '---  Convergence Failure  ---'
        WRITE(IWR,'(10X,A)') '---  Convergence Failure  ---'
!
!---    Isothermal option  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          NT = NSD(IEQT)
          IF( NT.GT.0 ) THEN
            NPT = MOD(NPHAZ(2,NT),10)
            NCHT = INDEX( PH_CND(NPT),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Energy Equation Maximum Residual = ',RSD(IEQT),
     &        ': Node = ',NT,
     &        ': Phase Condition = ',PH_CND(NPT)(1:NCHT)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Energy Equation Maximum Residual = ',RSD(IEQT),
     &        ': Node = ',NT,
     &        ': Phase Condition = ',PH_CND(NPT)(1:NCHT)
          ENDIF
        ENDIF
        NW = NSD(IEQW)
!
!---    Iso-CO2 option  ---
!
        IF( ISLC(45).EQ.0 ) NA = NSD(IEQA)
        NDO = NSD(IEQDO)
        IF( NW.GT.0 ) THEN
          NPW = MOD(NPHAZ(2,NW),10)
          NCHW = INDEX( PH_CND(NPW),'  ') - 1
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Equation Maximum Residual = ',RSD(IEQW),
     &      ': Node = ',NW,
     &      ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Equation Maximum Residual = ',RSD(IEQW),
     &      ': Node = ',NW,
     &      ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
        ENDIF
        IF( ISLC(45).EQ.0 ) THEN
          IF( NA.GT.0 ) THEN
            NPA = MOD(NPHAZ(2,NA),10)
            NCHA = INDEX( PH_CND(NPA),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'CO2 Equation Maximum Residual = ',RSD(IEQA),
     &        ': Node = ',NA,
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'CO2 Equation Maximum Residual = ',RSD(IEQA),
     &        ': Node = ',NA,
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
          ENDIF
        ENDIF
        IF( NDO.GT.0 ) THEN
          NPO = MOD(NPHAZ(2,NDO),10)
          NCHO = INDEX( PH_CND(NPO),'  ') - 1
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &      'CH4 Equation Maximum Residual = ',RSD(IEQDO),
     &      ': Node = ',NDO,
     &      ': Phase Condition = ',PH_CND(NPO)(1:NCHO)
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &      'CH4 Equation Maximum Residual = ',RSD(IEQDO),
     &      ': Node = ',NDO,
     &      ': Phase Condition = ',PH_CND(NPO)(1:NCHO)
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          NS = NSD(IEQS)
          IF( NS.GT.0 ) THEN
            NPS = MOD(NPHAZ(2,NS),10)
            NCHS = INDEX( PH_CND(NPS),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD(IEQS),
     &        ': Node = ',NS,
     &        ': Phase Condition = ',PH_CND(NPS)(1:NCHS)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD(IEQS),
     &         ': Node = ',NS,
     &        ': Phase Condition = ',PH_CND(NPS)(1:NCHS)
          ENDIF
        ENDIF
!
!---  Reduce time step  ---
!
        IF( NTSR.LT.4 .OR. (DTCF*DT).GT.DTMN ) THEN
          NTSR = NTSR + 1
          DTX = DT
!          DTCFX = 1.D+1**(1.D+0/(DTCF-1.D+0))
!          DTCFX = 1.D+0 + 1.D+0/(LOG10(MIN(DTCFX,RSDAX)))
!          TM = TM - (1.D+0-DTCFX)*DT
!          DT = DTCFX*DT
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
            XLA(2,N) = XLA(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            SH(2,N) = SH(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            PVA(2,N) = PVA(1,N)
            PVO(2,N) = PVO(1,N)
            YMGO(2,N) = YMGO(1,N)
#ifdef k1_eq
            VISDG(1,N) = VISDG(2,N)
#endif
#ifdef k2_eq
            VISDG(1,N) = VISDG(2,N)
            VISDG(3,N) = VISDG(4,N)
#endif
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
!---  End of RSDL37 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SBND37( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Modify the Jacobian matrix for the solute transport equation
!     to incorporate boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/SBND37'
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
!---  End of SBND37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORC37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Compute source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE JACOB
      USE GRID
      USE FDVT
      USE FDVS
      USE FDVP
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
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORC37'
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
                  XFOX = XFOGX + XFOLX
                  XFSLX = XLS(M,N)*RHOL(M,N)*SL(M,N)
                  XFSX = XFSLX
                  XFWGX = XGW(M,N)*RHOG(M,N)*SG(M,N)
                  XFWLX = XLW(M,N)*RHOL(M,N)*SL(M,N)
                  XFWX = XFWGX + XFWLX
                  XFTX = XFAX+XFOX+XFSX+XFWX
                  XFAX = XFAX/XFTX
                  XFOX = XFOX/XFTX
                  XFSX = XFSX/XFTX
                  XFWX = XFWX/XFTX
                  XFGX = (XFAGX+XFOGX+XFWGX)/XFTX
                  XFLX = (XFALX+XFOLX+XFSLX+XFWLX)/XFTX
                  XFNX = (XFANX)/XFTX
                  SRCA(M,N) = SRCA(M,N) + SRX(4)*XFAX
                  SRCO(M,N) = SRCO(M,N) + SRX(4)*XFOX
                  SRCW(M,N) = SRCW(M,N) + SRX(4)*XFWX
                  SRCS(M,N) = SRCS(M,N) + SRX(4)*XFSX
                  SRCT(M,N) = SRCT(M,N) + SRX(4)*
     &              (XFGX*HG(M,N) + XFLX*HL(M,N) + XFNX*HN(M,N))
!
!---            Initial-condition fluid mass rate  ---
!
                ELSEIF( ISRT(NS).EQ.4 ) THEN
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
!---  End of SORC37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORIT37( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Compute solute transport source integrals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/SORIT37'
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
!---  End of SORIT37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT37( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Compute solute transport source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      PetscErrorCode :: IERR
#endif
      REAL*8 SRX(8)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORT37'
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
!---  End of SORT37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPRP37( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Calculates the aqueous- and gas-phase solute
!     mole fractions from user-specified partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      SUB_LOG(ISUB_LOG) = '/SPRP37'
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
!---  End of SPRP37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT37
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     HYD (Gas Hydrate w/ Equilibrium Exchange) Mode
!
!     Update the primary variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 August 2004
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
      USE FDVH
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
      CHARACTER*64 PH_CND(6)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Aqueous-Gas w/o Hydrate',
     &  'Aqueous-Gas w/ Hydrate',
     &  'Aqueous-Low Liquid CO2-Gas w/o Hydrate',
     &  'Aqueous-Low Liquid CO2-Gas w Hydrate',
     &  'Aqueous-Liquid CO2-Gas w/o Hydrate',
     &  'Aqueous-Liquid CO2-Gas w Hydrate'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT37'
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
!
!---    Isothermal option  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          MPT = IM(IEQT,NMD)
          DPT = BLU(MPT)
        ELSE
          DPT = 0.D+0
        ENDIF
        MPW = IM(IEQW,NMD)
        DPW = BLU(MPW)
!
!---    Iso-CO2 option  ---
!
        IF( ISLC(45).EQ.0 ) THEN
          MPA = IM(IEQA,NMD)
          DPA = BLU(MPA)
        ELSE
          DPA = 0.D+0
        ENDIF
        MPO = IM(IEQDO,NMD)
        DPO = BLU(MPO)
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
!---    Limit changes in temperature to 1.0 C without hydrate, 
!       and to 0.25 C with hydrate, isothermal option   ---
!
        IF( ISLC(30).EQ.0 ) THEN
          IF( SH(2,N).GT.EPSL ) THEN
            DPT = SIGN( MIN( 2.5D-1,ABS(DPT) ),DPT )
          ELSE
            DPT = SIGN( MIN( 1.D+0,ABS(DPT) ),DPT )
          ENDIF
          IF( T(2,N).LT.0.D+0 ) DPT = 5.D-1*DPT
          T(2,N) = T(2,N) + DPT
        ENDIF
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
          IF( YLS(2,N)/EPSL.LT.EPSL .AND. BLU(MPS)/EPSL.LT.EPSL ) THEN
            BLU(MPS) = 0.D+0
            DPS = 0.D+0
          ENDIF
          YLS(2,N) = YLS(2,N)+DPS
          IF( YLS(2,N).LT.1.D-12 ) YLS(2,N) = 0.D+0
          XLS(2,N) = MIN( YLS(2,N),XLSMX )
        ENDIF
!
!---    Interfacial surface tensions  ---
!
        SIG_GL = SCHR(5,IZN)
        SIG_NL = SCHR(6,IZN)
        SIG_GN = SCHR(7,IZN)
        SIG_HL = SCHR(9,IZN)
        SIG_IL = SCHR(9,IZN)
        CALL SFT_L( T(2,N),XLS(2,N),SIG_GL )
        CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &    /(2.D+0*SIG_GL*SIG_NL)
        CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &    /(2.D+0*SIG_GL*SIG_GN)
!
!---    Interphase scaling factors  ---
!
        BGL = SCHR(8,IZN)/SIG_GL
        BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
        BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
        BHL = SCHR(8,IZN)/SCHR(9,IZN)
        BIL = SCHR(8,IZN)/SCHR(9,IZN)
!
!---    Entry pressure with hydrate barrier adjustment  ---
!
        HCPBFX = 1.D+0
        IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
          HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &      SIN(SCHR(12,IZN)*SH(2,N) + SCHR(13,IZN)) )
        ENDIF
        ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---    Aqueous or aqueous-gas, no hydrate, no liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - CO2 vapor pressure
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
!!
!!---      Limit changes in aqueous pressure  ---
!!
!          DPX = 2.5D-1*MAX( PL(2,N)+PATM,1.D+6 )
!          DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
!          PL(2,N) = PL(2,N) + DPW
!
!---      Zero negative corrections for zero CO2 partial pressure,
!         iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            IF( PVA(2,N)/EPSL.LT.EPSL .AND. 
     &        BLU(MPA)/EPSL.LT.EPSL ) THEN
              BLU(MPA) = 0.D+0
              DPA = 0.D+0
            ENDIF
            DPX = MAX( 1.D-1*PVA(2,N),1.D+4 )
            DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
            PVA(2,N) = MAX( (PVA(2,N)+DPA),0.D+0 )
            IF( PVA(2,N).LT.1.D-6 ) PVA(2,N) = 0.D+0
          ENDIF
!
!---      Zero negative corrections for zero CH4 partial pressure  ---
!
          IF( PVO(2,N)/EPSL.LT.EPSL .AND. BLU(MPO)/EPSL.LT.EPSL ) THEN
            BLU(MPO) = 0.D+0
            DPO = 0.D+0
          ENDIF
          DPX = MAX( 1.D-1*PVO(2,N),1.D+4 )
          DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
          PVO(2,N) = MAX( (PVO(2,N)+DPO),0.D+0 )
          IF( PVO(2,N).LT.1.D-6 ) PVO(2,N) = 0.D+0
!
!---      Limit changes in aqueous pressure  ---
!
          DPX = 5.D-2*(PG(2,N)+PATM)
          DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
          PL(2,N) = PL(2,N) + DPW
!!
!!---      Limit changes in aqueous pressure  ---
!!
!          IF( SG(1,N).LT.EPSL .AND. SG(2,N).LT.EPSL ) THEN
!            CALL SP_B( T(2,N),XLS(2,N),PSBX )
!            PPLX = MAX( PL(2,N)+PATM,PSBX )
!            CALL DENS_B( TX,PPLX,XLSX,RHOBX )
!            IF( ISLC(44).EQ.1 ) THEN
!             PVWX = PSBX
!            ELSE
!             CALL VPL_BX( TEQHX,PL(2,N),PVA(2,N),PVO(2,N),
!     &         PVWX,PSBX,RHOBX,XLS(2,N) )
!            ENDIF
!            PGX = PVA(2,N)+PVO(2,N)+PVWX-PATM
!!
!!---        Transition to gas  ---
!!
!            IF( (PL(2,N)+DPW).LT.(PGX-(ENPR/BGL)-1.D+3) ) THEN
!              PL(2,N) = PGX-(ENPR/BGL)-(1.D+2/(REAL(NITER)**3))
!!
!!---        No gas  ---
!!
!            ELSE
!              IF( ABS(DPW).GT.1.D+6 ) THEN
!                DPX = MIN( 6.D-1*ABS(DPW),1.D+6 )
!              ELSE
!                DPX = DPW
!              ENDIF
!              DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
!              PL(2,N) = PL(2,N) + DPW
!              PL(2,N) = MAX( PGX-(ENPR/BGL),PL(2,N) )
!            ENDIF
!!
!!---      Gas  ---
!!
!          ELSE
!             DPX = MAX( 1.D+5,1.D-1*(PG(2,N)-PL(2,N)) )
!             DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
!             PL(2,N) = PL(2,N) + DPW
!          ENDIF
!
!---    Aqueous or aqueous-gas, hydrate, no liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - Gas CH4 mole fraction of formers
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
!
!---      Limit changes in gas CH4 mole fraction of formers,
!         iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            DPX = 1.D-1
            DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
            YMGO(2,N) = MIN( MAX( YMGO(2,N)+DPA,0.D+0 ),1.D+0 )
          ENDIF
!
!---      Limit changes in hydrate saturation  ---
!
          DPX = 1.D-1
          DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
          IF( SL(2,N).LT.SCHR(15,IZN) ) DPO = MIN( DPO,0.D+0 )
          SH(2,N) = MIN( MAX( SH(2,N)+DPO,0.D+0 ),1.D+0 )
!
!---      Limit changes in aqueous pressure  ---
!
          DPX = 5.D-2*(PG(2,N)+PATM)
          DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
          PL(2,N) = PL(2,N) + DPW
!!
!!---      Limit changes in aqueous pressure  ---
!!
!          IF( SG(1,N).LT.EPSL .AND. SG(2,N).LT.EPSL ) THEN
!            TEQHX = T(2,N)
!            YMGOX = YMGO(2,N)
!            XLSX = XLS(2,N)
!            CALL HYD_P( TEQHX,YMGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
!            PVOX = YMGOX*PEQHX
!            PVAX = PEQHX - PVOX
!            CALL SP_B( TEQHX,XLSX,PSBX )
!            PPLX = MAX( PL(2,N)+PATM,PSBX )
!            CALL DENS_B( TX,PPLX,XLSX,RHOBX )
!            IF( ISLC(44).EQ.1 ) THEN
!             PVWX = PSBX
!            ELSE
!             CALL VPL_BX( TEQHX,PL(2,N),PVAX,PVOX,PVWX,PSBX,RHOBX,XLSX )
!            ENDIF
!            PGX = PVAX+PVOX+PVWX-PATM
!!
!!---        Transition to gas  ---
!!
!            IF( (PL(2,N)+DPW).LT.(PGX-(ENPR/BGL)-1.D+3) ) THEN
!              PL(2,N) = PGX-(ENPR/BGL)-(1.D+2/(REAL(NITER)**3))
!!
!!---        No gas  ---
!!
!            ELSE
!              IF( ABS(DPW).GT.1.D+6 ) THEN
!                DPX = MIN( 6.D-1*ABS(DPW),1.D+6 )
!              ELSE
!                DPX = DPW
!              ENDIF
!              DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
!              PL(2,N) = PL(2,N) + DPW
!              PL(2,N) = MAX( PGX-(ENPR/BGL),PL(2,N) )
!            ENDIF
!!
!!---      Gas  ---
!!
!          ELSE
!             DPX = MAX( 1.D+5,1.D-1*(PG(2,N)-PL(2,N)) )
!             DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
!             PL(2,N) = PL(2,N) + DPW
!          ENDIF
!
!---    Aqueous or aqueous-gas, no hydrate, low liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 saturation
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
!
!---      Limit changes in aqueous pressure  ---
!
          IF( ABS(DPW).GT.1.D+6 ) THEN
            DPX = MIN( 6.D-1*ABS(DPW),1.D+6 )
          ELSE
            DPX = DPW
          ENDIF
          DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
          PL(2,N) = PL(2,N) + DPW
!
!---      Limit changes in liquid-CO2 saturation,
!         iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            DPX = 1.D-1
            DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
            SN(2,N) = MIN( MAX( SN(2,N)+DPA,0.D+0 ),1.D+0 )
          ENDIF
!
!---      Zero negative corrections for zero CH4 partial pressure  ---
!
          IF( PVO(2,N)/EPSL.LT.EPSL .AND. BLU(MPO)/EPSL.LT.EPSL ) THEN
            BLU(MPO) = 0.D+0
            DPO = 0.D+0
          ENDIF
          DPX = MAX( 2.5D-1*PVO(2,N),1.D+6 )
          DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
          PVO(2,N) = MAX( (PVO(2,N)+DPO),0.D+0 )
          IF( PVO(2,N).LT.1.D-6 ) PVO(2,N) = 0.D+0
!
!---    Aqueous or aqueous-gas, hydrate, low liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 saturation
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
!
!---      Limit changes in aqueous pressure  ---
!
          IF( ABS(DPW).GT.1.D+6 ) THEN
            DPX = MIN( 6.D-1*ABS(DPW),1.D+6 )
          ELSE
            DPX = DPW
          ENDIF
          DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
          PL(2,N) = PL(2,N) + DPW
!
!---      Limit changes in liquid-CO2 saturation,
!         iso-CO2 option  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            DPX = 1.D-1
            DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
            SN(2,N) = MIN( MAX( SN(2,N)+DPA,0.D+0 ),1.D+0 )
          ENDIF
!
!---      Limit changes in hydrate saturation  ---
!
          DPX = 1.D-1
          DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
          IF( SL(2,N).LT.SCHR(15,IZN) ) DPO = MIN( DPO,0.D+0 )
          SH(2,N) = MIN( MAX( SH(2,N)+DPO,0.D+0 ),1.D+0 )
!
!---    Aqueous or aqueous-gas, no hydrate, liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 pressure
!       CH4 mass - CH4 vapor pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
!
!---      Limit changes in pressure  ---
!
          DPX = 5.D-2*(PG(2,N)+PATM)
!
!---      Limit changes in aqueous pressure  ---
!
!          IF( ABS(DPW).GT.1.D+6 ) THEN
!            DPX = MIN( 6.D-1*ABS(DPW),1.D+6 )
!          ELSE
!            DPX = DPW
!          ENDIF
          DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
          PL(2,N) = PL(2,N) + DPW
!
!---      Limit changes in liquid-CO2 pressure  ---
!
!          IF( ABS(DPA).GT.1.D+6 ) THEN
!            DPX = MIN( 6.D-1*ABS(DPA),1.D+6 )
!          ELSE
!            DPX = DPA
!          ENDIF
          DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
          PN(2,N) = PN(2,N) + DPA
!
!---      Zero negative corrections for zero CH4 partial pressure  ---
!
          IF( PVO(2,N)/EPSL.LT.EPSL .AND. BLU(MPO)/EPSL.LT.EPSL ) THEN
            BLU(MPO) = 0.D+0
            DPO = 0.D+0
          ENDIF
          DPX = MAX( 1.D-1*PVO(2,N),1.D+4 )
          DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
          PVO(2,N) = MAX( (PVO(2,N)+DPO),0.D+0 )
          IF( PVO(2,N).LT.1.D-6 ) PVO(2,N) = 0.D+0
!
!---    Aqueous or aqueous-gas, hydrate, liquid CO2
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       CO2 mass - liquid CO2 pressure
!       CH4 mass - hydrate saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
!
!---      Limit changes in pressure  ---
!
          DPX = 5.D-2*(PG(2,N)+PATM)
!
!---      Limit changes in aqueous pressure  ---
!
!          IF( ABS(DPW).GT.1.D+6 ) THEN
!            DPX = MIN( 6.D-1*ABS(DPW),1.D+6 )
!          ELSE
!            DPX = DPW
!          ENDIF
          DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
          PL(2,N) = PL(2,N) + DPW
!
!---      Limit changes in liquid-CO2 pressure  ---
!
!          IF( ABS(DPA).GT.1.D+6 ) THEN
!            DPX = MIN( 6.D-1*ABS(DPA),1.D+6 )
!          ELSE
!            DPX = DPA
!          ENDIF
          DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
          PN(2,N) = PN(2,N) + DPA
!
!---      Limit changes in hydrate saturation  ---
!
          DPX = 1.D-1
          DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
          IF( SL(2,N).LT.SCHR(15,IZN) ) DPO = MIN( DPO,0.D+0 )
          SH(2,N) = MIN( MAX( SH(2,N)+DPO,0.D+0 ),1.D+0 )
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
          WRITE(ISC,'(4X,A,I6)') 'Node = ',N
          WRITE(IWR,'(4X,A,I6)') 'Node = ',N
          WRITE(ISC,'(4X,2A)') 'Phase Condition = ',PH_CND(NPHAZ(2,N))
          WRITE(IWR,'(4X,2A)') 'Phase Condition = ',PH_CND(NPHAZ(2,N))
          WRITE(ISC,'(4X,A,1PE12.5)') 'Temperature = ',T(2,N)
          WRITE(IWR,'(4X,A,1PE12.5)') 'Temperature = ',T(2,N)
          WRITE(ISC,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',PL(2,N)+PATM
          WRITE(IWR,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',PL(2,N)+PATM
!
!---      Aqueous or aqueous-gas, no hydrate, no liquid CO2
!
          IF( NPHAZ(2,N).EQ.1 ) THEN
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              WRITE(ISC,'(4X,A,1PE12.5)')
     &          'CO2 Partial Pressure = ',PVA(2,N)
              WRITE(IWR,'(4X,A,1PE12.5)')
     &          'CO2 Partial Pressure = ',PVA(2,N)
            ENDIF
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'CH4 Partial Pressure = ',PVO(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'CH4 Partial Pressure = ',PVO(2,N)
!
!---      Aqueous or aqueous-gas, hydrate, no liquid CO2
!
          ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              WRITE(ISC,'(4X,A,1PE12.5)')
     &          'Gas CH4 Mole Fraction of Formers = ',YMGO(2,N)
              WRITE(IWR,'(4X,A,1PE12.5)')
     &          'Gas CH4 Mole Fraction of Formers = ',YMGO(2,N)
            ENDIF
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Hydrate Saturation = ',SH(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Hydrate Saturation = ',SH(2,N)
!
!---      Aqueous or aqueous-gas, no hydrate, low liquid CO2
!
          ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              WRITE(ISC,'(4X,A,1PE12.5)')
     &          'Liquid-CO2 Saturation = ',SN(2,N)
              WRITE(IWR,'(4X,A,1PE12.5)')
     &          'Liquid-CO2 Saturation = ',SN(2,N)
            ENDIF
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'CH4 Partial Pressure = ',PVO(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'CH4 Partial Pressure = ',PVO(2,N)
!
!---      Aqueous or aqueous-gas, hydrate, low liquid CO2
!
          ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              WRITE(ISC,'(4X,A,1PE12.5)')
     &          'Liquid-CO2 Saturation = ',SN(2,N)
              WRITE(IWR,'(4X,A,1PE12.5)')
     &          'Liquid-CO2 Saturation = ',SN(2,N)
            ENDIF
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Hydrate Saturation = ',SH(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Hydrate Saturation = ',SH(2,N)
!
!---      Aqueous or aqueous-gas, no hydrate, liquid CO2
!
          ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              WRITE(ISC,'(4X,A,1PE12.5)')
     &          'Liquid-CO2 Pressure = ',PN(2,N)
              WRITE(IWR,'(4X,A,1PE12.5)')
     &          'Liquid-CO2 Pressure = ',PN(2,N)
            ENDIF
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'CH4 Partial Pressure = ',PVO(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'CH4 Partial Pressure = ',PVO(2,N)
!
!---      Aqueous or aqueous-gas, hydrate, liquid CO2
!
          ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
!
!---        Iso-CO2 option  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              WRITE(ISC,'(4X,A,1PE12.5)')
     &          'Liquid-CO2 Pressure = ',PN(2,N)
              WRITE(IWR,'(4X,A,1PE12.5)')
     &          'Liquid-CO2 Pressure = ',PN(2,N)
            ENDIF
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Hydrate Saturation = ',SH(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Hydrate Saturation = ',SH(2,N)
          ENDIF
          WRITE(ISC,'(4X,A,1PE12.5,A,I6)')
     &      'Total-Salt Aqu. Mass Fraction = ',YLS(2,N)
          WRITE(IWR,'(4X,A,1PE12.5,A,I6)')
     &      'Total-Salt Aqu. Mass Fraction = ',YLS(2,N)
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
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            XLA(2,N) = XLA(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            SH(2,N) = SH(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            PVA(2,N) = PVA(1,N)
            PVO(2,N) = PVO(1,N)
            YMGO(2,N) = YMGO(1,N)
#ifdef k1_hyd
            VISDG(1,N) = VISDG(2,N)
#endif
#ifdef k2_hyd
            VISDG(1,N) = VISDG(2,N)
            VISDG(3,N) = VISDG(4,N)
#endif
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
            XLA(2,N) = XLA(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            SH(2,N) = SH(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            PVA(2,N) = PVA(1,N)
            PVO(2,N) = PVO(1,N)
            YMGO(2,N) = YMGO(1,N)
#ifdef k1_hyd
            VISDG(1,N) = VISDG(2,N)
#endif
#ifdef k2_hyd
            VISDG(1,N) = VISDG(2,N)
            VISDG(3,N) = VISDG(4,N)
#endif
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
!---  End of UPDT37 group
!
      RETURN
      END

!----------------------Program-----------------------------------------!
!
      PROGRAM STOMP4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     STOMP: Subsurface Transport Over Multiple Phases
!
!     Water-Oil Mode
!
!     This operational mode of the STOMP simulator models the flow
!     of aqueous and non-aqueous phase liquid through variably
!     saturated porous media.  Gas phase flow is not considered.
!     The aqueous phase is considered to comprise liquid water and
!     dissolved-oil, whereas water is considered immiscible in the
!     nonaqueous phase.  Solute transport is modeled through the
!     aqueous and nonaqueous phases with sorption to the geologic
!     media.
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
      USE REACT
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
#ifdef awtap
      REAL*8 HCFL(LFZ,LSR),HCFN(LFZ,LSR)
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
      SUB_LOG(1) = 'STOMP4'
      ICODE = 4
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
      CALL RDINPT4
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
!---  Define surface spill grid  ---
!
      IF( ISLC(49).EQ.1 ) CALL SP_GRID
!
!---  Check thermodynamic and hydrologic initial states  ---
!
      CALL CHK4
!
!---  For geomechanics set k iterate value of pore pressure  ---
!
      IF( ISLC(50).NE.0 ) THEN
        INDX = 2
        CALL PRESS_GM( INDX )
      ENDIF
#ifdef ecke
!
!---  Sequence reaction equations  ---
!
      IF( ISLC(40).EQ.1 ) CALL SEQEQ
#endif
!
!---  Initialize well states  ---
!
      IF( LWELL.EQ.1 ) CALL CHKW4
!
!---  Compute Jacobian matrix pointers
!     with surface spills  ---
!
      IF( ISLC(49).EQ.1 ) THEN
        CALL JCBPSP
      
!---  Compute Jacobian matrix pointers 
!     without surface spills  ---
!
      ELSE
        CALL JCBP
      ENDIF
!
!---  For geomechanics simulations compute Jacobian matrix pointers  --
!
      IF( ISLC(50).NE.0 .AND. ISLC(67).EQ.0 ) CALL JCBP_GM
!
!---  Compute primary variable increments  ---
!
      CALL INCRM4
!
!---  Surface spill head increments  ---
!
      IF( ISLC(49).EQ.1 ) CALL INCRMSP
!
!---  Saturation, relative permeability, porosity, and tortuosity  ---
!
      CALL PROP4( INDX )
      IF( INDX.LT.0 ) THEN
        WRITE(IWR,'(/,A)') '---  End of STOMP4 Simulation  ---'
        WRITE(ISC,'(/,A)') '---  End of STOMP4 Simulation  ---'
        STOP
      ENDIF
!
!---  Thermodynamic properties and equations of state  ---
!
      CALL TMPR4
      IF( LWELL.EQ.1 ) CALL PROPW4
!
!---  Surface spill properties  ---
!
      IF( ISLC(49).EQ.1 ) CALL PROPSP
!
!---  Compute boundary saturation, relative permeability, and
!     thermodynamic properties  ---
!
      CALL BCP4
!
!---  Compute initial solute concentrations  ---
!
      CALL CISC4
#ifdef ecke
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
#endif
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
!---  Compute initial aqueous-phase volumetric flux,
!     NAPL volumetric flux, and dissolved-oil flux  ---
!
      ISVF = 1
      CALL DRCVL
      CALL DRCVN
      IF( ISLC(12).EQ.1 ) THEN
        CALL DDFLO
      ELSE
        CALL DFFLO
      ENDIF
      CALL BCF4
!
!---  Initial NAPL surface spill flux and infiltration  ---
!
      IF( ISLC(49).EQ.1 ) CALL FLXSPN
!
!---  Initial aqueous surface spill flux and infiltration  ---
!
      IF( ISLC(49).EQ.1 ) CALL FLXSPL
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
      CALL LDO4
!
!---  Load old time step well arrays  ---
!
      IF( LWELL.EQ.1 ) CALL LDOW4
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
      CALL TRPNL4
!
!---  Reference node(s) output  ---
!
      IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 .OR.
     &  MOD( (NSTEP-NRST),IFQO ).EQ.0) CALL REFNOD
#ifdef ecke
!
!---  Normalize mineral species concentrations after initial
!     output for normal simulations  ---
!
      IF( (NSTEP-NRST).EQ.0 ) CALL NMNSP
#endif
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
        CALL BCP4
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
      CALL BCP4
!
!---  Compute well boundary thermodynamic properties  ---
!
      IF( LWELL.EQ.1 ) CALL BCPW4
!
!---  Compute NAPL surface spill source contributions  ---
!
      IF( ISLC(49).EQ.1 ) CALL SORCSP4
!
!---  Compute source contributions  ---
!
#ifdef awtap
      CALL SORC4( HCFL,HCFN )
#else
      CALL SORC4
#endif
!
!---  Compute aqueous-phase volumetric flux (interior surfaces)  ---
!
      CALL DRCVL
!
!---  Compute NAPL volumetric flux (interior surfaces)  ---
!
      CALL DRCVN
!
!---  Courant number limit check  ---
!
      IF( CRNTMXC.GT.EPSL ) THEN
        ICNV = 3
        CALL CRN_L
        CALL CRN_N
!
!---    Cut time step  ---
!
        IF( ICNV.EQ.1 ) GOTO 200
      ENDIF
!
!---  Compute aqueous-phase volumetric
!     and NAPL volumetric well fluxes
!     (interior surfaces)  ---
!
      IF( LWELL.EQ.1 ) CALL FLUXW4
!
!---  Compute oil diffusion flux through the aqueous phase
!     (interior surfaces)  ---
!
      IF( ISLC(12).EQ.1 ) THEN
        CALL DDFLO
      ELSE
        CALL DFFLO
      ENDIF
!
!---  Compute well water vapor diffusive flux through the gas phase
!     and well air diffusive flux through the aqueous phase
!     (interior surfaces)  ---
!
      IF( LWELL.EQ.1 ) CALL DFLUXW4
!
!---  Compute aqueous-phase, and NAPL volumetric flux
!     (boundary surfaces)  ---
!
      CALL BCF4
!
!---  Compute aqueous-phase and NAPL well fluxes
!     (boundary surfaces)  ---
!
      IF( LWELL.EQ.1 ) CALL BCFW4
!
!---  Initial NAPL surface spill flux and infiltration  ---
!
      IF( ISLC(49).EQ.1 ) CALL FLXSPN
!
!---  Initial aqueous surface spill flux and infiltration  ---
!
      IF( ISLC(49).EQ.1 ) CALL FLXSPL
!
!---  Compute well source contributions  ---
!
      IF( LWELL.EQ.1 ) CALL SORCW4
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
      CALL JCBWLN
!
!---  Load Jacobian matrix for the oil equation
!     (zero flux boundary)  ---
!
      CALL JCBOLN
!
!---  Load Jacobian matrix for the well-water equation and
!     modify the Jacobian matrix for the field-water for
!     well flow  ---
!
      IF( LWELL.EQ.1 ) CALL JCBWW4
!
!---  Load Jacobian matrix for the well-oil equation and
!     modify the Jacobian matrix for the field-oil for
!     well flow  ---
!
      IF( LWELL.EQ.1 ) CALL JCBOW4
!
!---  Load Jacobian matrix for water equation
!     spill surface nodes  ---
!
      IF( ISLC(49).EQ.1 ) CALL JCSPW
!
!---  Load Jacobian matrix for oil equation
!     spill surface nodes  ---
!
      IF( ISLC(49).EQ.1 ) CALL JCSPO
!
!---  Modify the Jacobian matrix for the boundary conditions  ---
!
      CALL BCJ4
!
!---  Modify the Jacobian matrix for the well
!     boundary conditions  ---
!
      IF( LWELL.EQ.1 ) CALL BCJW4
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
      CALL UPDT4
!
!---  Update well primary variables  ---
!
      IF( LWELL.EQ.1 ) CALL UPDTW4
!
!---  Update surface spill heads  ---
!
      IF( ISLC(49).EQ.1 ) CALL UPDTSP
!
!---  Compute convergence from maximum relative residuals  ---
!
      CALL RSDL4
!
!---  Compute well equation residuals  ---
!
      IF( LWELL.EQ.1 ) CALL RSDLW4
!
!---  Compute surface spill residuals  ---
!
      IF( ISLC(49).EQ.1 ) CALL RSDSP
!
!---  Compute primary variable increments, saturation,
!     relative permeability, porosity, tortuosity,
!     thermodynamic properties for interior nodes,
!     except immediately after a new time step  ---
!
  400 CONTINUE
      CALL INCRM4
      CALL PROP4( INDX )
      IF( INDX.LT.0 ) GOTO 400
      CALL TMPR4
!
!---  Compute well properties  ---
!
      IF( LWELL.EQ.1 ) CALL PROPW4
!
!---  Surface spill head increments  ---
!
      IF( ISLC(49).EQ.1 ) CALL INCRMSP
!
!---  Surface spill properties  ---
!
      IF( ISLC(49).EQ.1 ) CALL PROPSP
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
!---  Compute aqueous-phase, and napl volumetric flux
!     (interior surfaces)  ---
!
      ISVF = 1
      CALL DRCVL
      CALL DRCVN
!
!---  Compute aqueous-phase, and napl volumetric flux
!     (boundary surfaces)  ---
!
      CALL BCF4
!
!---  Initial NAPL surface spill flux and infiltration  ---
!
      IF( ISLC(49).EQ.1 ) CALL FLXSPN
!
!---  Initial aqueous surface spill flux and infiltration  ---
!
      IF( ISLC(49).EQ.1 ) CALL FLXSPL
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
          CALL SPRP4( NSL )
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
#ifdef awtap
          CALL SORT4( HCFL,HCFN,NSL )
#else
          CALL SORT4( NSL )
#endif
!
!---      Zero solute transport fluxes  ---
!
          CALL SFXZ( NSL )
!
!---      Load Jacobian matrix (aqueous-phase transport)  ---
!
          CALL SJCBL( NSL )
!
!---      Load Jacobian matrix (NAPL transport)  ---
!
          CALL SJCBN( NSL )
!
!---      Modify Jacobian matrix for boundary conditions ---
!
          CALL SBND4( NSL )
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
!---      Compute solute NAPL fluxes (interior nodes)  ---
!
          CALL SFXN( NSL )
!
!---      Compute solute NAPL fluxes (boundary surfaces)  ---
!
          CALL SFXNB( NSL )
!
!---      Integrate solute sources  ---
!
#ifdef awtap
          CALL SORIT4( HCFL,HCFN,NSL )
#else
          CALL SORIT4( NSL )
#endif
!
!---      Load old sub-time-step concentrations  ---
!
          IF( ISLC(17).NE.0 ) CALL UPDTCO( NSL)
!
!---    Bottom of sub-time step loop  ---
!
  690   CONTINUE
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
#ifdef ecke
!
!---  Reactive transport  ---
!
      IF( ISLC(40).EQ.1 ) THEN
!
!---    Temporarily store time stepping  ---
!
        DT_RST = DT
        DTI_RST = DTI
        TM_RST = TM
        TM = TM - DT
        N_RST = 1
  710   CONTINUE
!
!---    Zero linked sources  ---
!
        CALL ZLKSRC
!
!---    Sub-time step reduction limit exceeded  ---
!
        IF( N_RST.GT.16 ) THEN
          WRITE(ISC,'(A)') '          ---  ECKEChem ' // 
     &      'Sub-Time Step Reduction Limit Exceeded  ---'
          WRITE(IWR,'(A)') '          ---  ECKEChem ' // 
     &      'Sub-Time Step Reduction Limit Exceeded  ---'
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
!---    Sub-time step loop  ---
!
        DO 790 NC = 1,N_RST
          TM = TM + DT
!
!---      Loop over number of conservation component species  ---
!
          DO 730 NEQ = 1,NEQC
            NSL = NEQ + NSOLU
!
!---        Skip transport for linked aqueous oil   ---
!
            IF( ISPLK(12).EQ.NSL ) GOTO 720
!
!---        Mobile conservation component fractions   ---
!
            CALL MOBCF( NEQ )
!
!---        Solute transport ---
!
            CALL TPORT32( NSL )
!
!---        Add immobile conservation component fractions   ---
!
  720       CONTINUE
            CALL IMOBCF( NEQ )
!
!---      End of conservation component species transport  ---
!
  730     CONTINUE
!
!---      Loop over number of kinetic component species  ---
!
          DO 750 NEQ = 1,NEQK
            NSL = NEQ + NEQC + NSOLU
!
!---        Skip transport for linked aqueous oil   ---
!
            IF( ISPLK(12).EQ.NSL ) GOTO 740
! 
!---        Mobile kinetic component fractions   ---
!
            CALL MOBKF( NEQ )
!
!---        Solute transport ---
!
            CALL TPORT32( NSL )
! 
!---        Add immobile kinetic component fractions   ---
!
  740       CONTINUE
            CALL IMOBKF( NEQ )
!
!---      End of conservation component species transport  ---
!
  750     CONTINUE
!
!---      Equilibrium-conservation-kinetic reaction chemistry   ---
!
          CALL ECKECHEM
          IF( ECKE_ER ) GOTO 710
!
!---      Load old sub-time-step reactive species
!         concentrations and component species concentrations  ---
!
          IF( ISLC(17).NE.0 ) CALL UPDTCHEM
!
!---    Bottom of sub-time step loop  ---
!
  790   CONTINUE
!
!---    Reset time stepping  ---
!
        IF( N_RST.GT.1 ) THEN
          DT = DT_RST
          DTI = DTI_RST
          TM = TM_RST
        ENDIF
      ENDIF
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
!---  Integrate well pumping rates  ---
!
      IF( LWELL.EQ.1 ) CALL SFINW4
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
      STOP
!
!---  End of STOMP program  ---
!
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCF4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Compute boundary surface fluxes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/BCF4'
!
!---  Zero boundary fluxes  ---
!
      DO 70 NB = 1,NBC
        N = IBCN(NB)
        IF( IBCD(NB).EQ.-3 ) THEN
          DO 10 M = 1,ISVF
            WL(M,NSZ(N)) = 0.D+0
            WN(M,NSZ(N)) = 0.D+0
            WDLO(M,NSZ(N)) = 0.D+0
            WLO(M,NSZ(N)) = 0.D+0
   10     CONTINUE
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          DO 20 M = 1,ISVF
            VL(M,NSY(N)) = 0.D+0
            VN(M,NSY(N)) = 0.D+0
            VDLO(M,NSY(N)) = 0.D+0
            VLO(M,NSY(N)) = 0.D+0
   20     CONTINUE
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          DO 30 M = 1,ISVF
            UL(M,NSX(N)) = 0.D+0
            UN(M,NSX(N)) = 0.D+0
            UDLO(M,NSX(N)) = 0.D+0
            ULO(M,NSX(N)) = 0.D+0
   30     CONTINUE
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          DO 40 M = 1,ISVF
            UL(M,NSX(N)+1) = 0.D+0
            UN(M,NSX(N)+1) = 0.D+0
            UDLO(M,NSX(N)+1) = 0.D+0
            ULO(M,NSX(N)+1) = 0.D+0
   40     CONTINUE
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          DO 50 M = 1,ISVF
            VL(M,NSY(N)+IFLD) = 0.D+0
            VN(M,NSY(N)+IFLD) = 0.D+0
            VDLO(M,NSY(N)+IFLD) = 0.D+0
            VLO(M,NSY(N)+IFLD) = 0.D+0
   50     CONTINUE
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          DO 60 M = 1,ISVF
            WL(M,NSZ(N)+IJFLD) = 0.D+0
            WN(M,NSZ(N)+IJFLD) = 0.D+0
            WDLO(M,NSZ(N)+IJFLD) = 0.D+0
            WLO(M,NSZ(N)+IJFLD) = 0.D+0
   60     CONTINUE
        ENDIF
   70 CONTINUE
!
!---  Loop over boundary conditions  ---
!
      DO 200 NB = 1,NBC
!
!---  Special decreasing pond boundary condition for the MOFHYS
!     verification problem  ---
!
!        MB = IBCIN(NB)
!        IF( NB.EQ.2 .AND. BC(3,1,MB).LE.SMALL ) GOTO 200
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
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 110 M = 1,ISVF
              WL(M,NSZ(N)) = BCX(2)
  110       CONTINUE
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLOB( N,NB )
            ELSE
              CALL DFFLOB( N,NB )
            ENDIF
          ELSEIF( IBCT(IEQW,NB).EQ.7 ) THEN
            CALL DRCVLB( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLB( N,NB )
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLOB( N,NB )
            ELSE
              CALL DFFLOB( N,NB )
            ENDIF
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQO,NB).EQ.2 ) THEN
            DO 112 M = 1,ISVF
              WN(M,NSZ(N)) = BCX(3)
  112       CONTINUE
          ELSEIF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL DRCVNB( N,NB )
          ENDIF
!
!---  South boundary
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---  Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 120 M = 1,ISVF
              VL(M,NSY(N)) = BCX(2)
  120       CONTINUE
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLOS( N,NB )
            ELSE
              CALL DFFLOS( N,NB )
            ENDIF
          ELSEIF( IBCT(IEQW,NB).EQ.7 ) THEN
            CALL DRCVLS( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLS( N,NB )
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLOS( N,NB )
            ELSE
              CALL DFFLOS( N,NB )
            ENDIF
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQO,NB).EQ.2 ) THEN
            DO 122 M = 1,ISVF
              VN(M,NSY(N)) = BCX(3)
  122       CONTINUE
          ELSEIF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL DRCVNS( N,NB )
          ENDIF
!
!---  West boundary
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---  Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 130 M = 1,ISVF
              UL(M,NSX(N)) = BCX(2)
  130       CONTINUE
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLOW( N,NB )
            ELSE
              CALL DFFLOW( N,NB )
            ENDIF
          ELSEIF( IBCT(IEQW,NB).EQ.7 ) THEN
            CALL DRCVLW( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLW( N,NB )
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLOW( N,NB )
            ELSE
              CALL DFFLOW( N,NB )
            ENDIF
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQO,NB).EQ.2 ) THEN
            DO 132 M = 1,ISVF
              UN(M,NSX(N)) = BCX(3)
  132       CONTINUE
          ELSEIF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL DRCVNW( N,NB )
          ENDIF
!
!---  East boundary
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---  Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 140 M = 1,ISVF
              UL(M,NSX(N)+1) = BCX(2)
  140       CONTINUE
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLOE( N,NB )
            ELSE
              CALL DFFLOE( N,NB )
            ENDIF
          ELSEIF( IBCT(IEQW,NB).EQ.7 ) THEN
            CALL DRCVLE( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLE( N,NB )
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLOE( N,NB )
            ELSE
              CALL DFFLOE( N,NB )
            ENDIF
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQO,NB).EQ.2 ) THEN
            DO 142 M = 1,ISVF
              UN(M,NSX(N)+1) = BCX(3)
  142       CONTINUE
          ELSEIF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL DRCVNE( N,NB )
          ENDIF
!
!---  North boundary
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---  Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 150 M = 1,ISVF
              VL(M,NSY(N)+IFLD) = BCX(2)
  150       CONTINUE
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLON( N,NB )
            ELSE
              CALL DFFLON( N,NB )
            ENDIF
          ELSEIF( IBCT(IEQW,NB).EQ.7 ) THEN
            CALL DRCVLN( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLN( N,NB )
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLON( N,NB )
            ELSE
              CALL DFFLON( N,NB )
            ENDIF
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQO,NB).EQ.2 ) THEN
            DO 152 M = 1,ISVF
              VN(M,NSY(N)+IFLD) = BCX(3)
  152       CONTINUE
          ELSEIF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL DRCVNN( N,NB )
          ENDIF
!
!---  Top boundary
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---  Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 160 M = 1,ISVF
              WL(M,NSZ(N)+IJFLD) = BCX(2)
  160       CONTINUE
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLOT( N,NB )
            ELSE
              CALL DFFLOT( N,NB )
            ENDIF
          ELSEIF( IBCT(IEQW,NB).EQ.7 ) THEN
            CALL DRCVLT( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLT( N,NB )
            IF( ISLC(12).EQ.1 ) THEN
              CALL DDFLOT( N,NB )
            ELSE
              CALL DFFLOT( N,NB )
            ENDIF
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQO,NB).EQ.2 ) THEN
            DO 162 M = 1,ISVF
              WN(M,NSZ(N)+IJFLD) = BCX(3)
  162       CONTINUE
          ELSEIF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL DRCVNT( N,NB )
          ENDIF
        ENDIF
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCF4 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCJ4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Modify the Jacobian matrix for boundary conditions
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/BCJ4'
!
!---  Loop over boundary conditions  ---
!
      DO 100 NB = 1,NBC
!
!---  Special decreasing pond boundary condition for the MOFHYS
!     verification problem  ---
!
!        MB = IBCIN(NB)
!        IF( NB.EQ.2 .AND. BC(3,1,MB).LE.SMALL ) GOTO 100
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
!---  Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---  Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWB( N,NB,NPZ )
            CALL JCBLOB( N,NB,NPZ )
          ENDIF
!
!---  NAPL  ---
!
          IF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL JCBNOB( N,NB,NPZ )
          ENDIF
!
!---  South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---  Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWS( N,NB,NPY )
            CALL JCBLOS( N,NB,NPY )
          ENDIF
!
!---  NAPL
!
          IF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL JCBNOS( N,NB,NPY )
          ENDIF
!
!---  West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---  Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWW( N,NB,NPX )
            CALL JCBLOW( N,NB,NPX )
          ENDIF
!
!---  NAPL  ---
!
          IF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL JCBNOW( N,NB,NPX )
          ENDIF
!
!---  East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---  Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWE( N,NB,NQX )
            CALL JCBLOE( N,NB,NQX )
          ENDIF
!
!---  NAPL  ---
!
          IF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL JCBNOE( N,NB,NQX )
          ENDIF
!
!---  North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---  Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWN( N,NB,NQY )
            CALL JCBLON( N,NB,NQY )
          ENDIF
!
!---  NAPL  ---
!
          IF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL JCBNON( N,NB,NQY )
          ENDIF
!
!---  Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---  Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWT( N,NB,NQZ )
            CALL JCBLOT( N,NB,NQZ )
          ENDIF
!
!---  NAPL  ---
!
          IF( IBCT(IEQO,NB).NE.3 ) THEN
            CALL JCBNOT( N,NB,NQZ )
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCJ4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCP4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Compute saturation, relative permeability and thermodynamic
!     properties for boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     Last Modified by MD White, Battelle, November 8, 1999.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE FDVD
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/BCP4'
      INDX = 1
!
!---  Assign values for initial condition type boundary conditions  ---
!
      IF( NSTEP-NRST.LE.1 .AND. NITER.LE.1 ) THEN
       DO 50 NB = 1,NBC
        IF( IBCT(IEQW,NB).EQ.12 .OR. IBCT(IEQO,NB).EQ.12 ) THEN
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
        IF( IBCT(IEQW,NB).EQ.12 ) PLB(1,NB) = PL(2,N) + RHOL(2,N)*GB
        IF( IBCT(IEQO,NB).EQ.12 ) PNB(1,NB) = PN(2,N) + RHON(2,N)*GB
   50  CONTINUE
      ENDIF
!
!---  Loop over boundary conditions  ---
!
      DO 400 NB = 1,NBC
!!
!!---    Special decreasing pond boundary condition for the MOFHYS
!!       verification problem  ---
!!
!        MB = IBCIN(NB)
!        IF( NB.EQ.2 .AND. BC(3,1,MB).LE.SMALL ) GOTO 400
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 400
!
!---    Assign local boundary condition variables  ---
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
             IF( IBCT(IEQW,NB).EQ.2 ) THEN
               BCX(2) = BCX(2)-5.D-1*DTBC*(BC(2,M,MB)-BC(2,M-1,MB))/TDBC
             ENDIF
             IF( IBCT(IEQO,NB).EQ.2 ) THEN
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
        IF( IBCT(IEQW,NB).EQ.12 ) BCX(2) = PLB(1,NB)
        IF( IBCT(IEQO,NB).EQ.12 ) BCX(3) = PNB(1,NB)
        N = IBCN(NB)
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX+1
        NQY = NPY+IFLD
        NQZ = NPZ+IJFLD
        IZN = IZ(N)
        IBD = ABS(IBCD(NB))
#ifdef ecke
        POR0(1,N) = POR0(1,N)
        POR0(2,N) = POR0(2,N)
#endif
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
!---      Loop over secondary variable indices  ---
!
        DO 300 M = 2,ISVC+2
          TX = T(2,N)
          PLX = PL(M,N)
          PGX = PG(2,N)
          PNX = PN(M,N)
!
!---      Aqueous Dirichlet  ---
!
          IF( IBCT(IEQW,NB).EQ.1 .OR. IBCT(IEQW,NB).EQ.12 ) THEN
            PLX = BCX(2)
!
!---      Aqueous Neumann  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.2 ) THEN
            PLX = PLX + BCX(2)*DB*VISL(M,N)/PERM(IBD,IZN)
     &          + RHOL(M,N)*GB
!
!---      Aqueous Zero Flux  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.3 ) THEN
            IF( ABS(BCX(2)+PATM).GT.EPSL ) THEN
              PLX = BCX(2)
            ELSE
              PLX = PLX + RHOL(M,N)*GB
            ENDIF
!
!---      Aqueous Unit Gradient  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.5 ) THEN
            PLX = PLX
!
!---      Aqueous Outflow  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.7 ) THEN
            PLX = MIN( PLX+RHOL(M,N)*GB,PGX )
!
!---      Aqueous Hydraulic Gradient Base Node  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.-11 ) THEN
            PLX = BCX(2)
            NHG = N
!
!---      Aqueous Hydraulic Gradient  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.11 ) THEN
            PLX = FNHGBL( NBHG(1,NB),NB,M )
          ENDIF
!
!---      NAPL Dirichlet  ---
!
          IF( IBCT(IEQO,NB).EQ.1 .OR. IBCT(IEQO,NB).EQ.12 ) THEN
            PNX = BCX(3)
!
!---      NAPL Neumann; Logic prevents NAPL flow 
!         from two-phase nodes  ---
!
          ELSEIF( IBCT(IEQO,NB).EQ.2 ) THEN
            IF( NPHAZ(2,N).LE.4 .AND. BCX(3)*DB.LT.ZERO ) THEN
              INDX = 12
              CHMSG = 'NAPL Boundary Flux: Boundary Node: '
              IMSG = N
              CALL WRMSGS( INDX )
            ELSE
              PNX = PNX + BCX(3)*DB*VISN(M,N)/PERM(IBD,IZN)
     &          + RHON(M,N)*GB
            ENDIF
!
!---      NAPL Zero Flux  ---
!
          ELSEIF( IBCT(IEQO,NB).EQ.3 ) THEN
            IF( ABS(BCX(3)+PATM).GT.EPSL ) THEN
              PNX = BCX(3)
            ELSE
              PNX = PNX + RHON(M,N)*GB
            ENDIF
!
!---      NAPL Unit Gradient  ---
!
          ELSEIF( IBCT(IEQO,NB).EQ.5 ) THEN
            PNX = PNX
!
!---      NAPL Outflow  ---
!
          ELSEIF( IBCT(IEQO,NB).EQ.7 ) THEN
            PNX = MIN( PNX+RHON(M,N)*GB,PGX )
!
!---      NAPL Hydraulic Gradient Base Node  ---
!
          ELSEIF( IBCT(IEQO,NB).EQ.-11 ) THEN
            PNX = BCX(3)
            NHG = N
!
!---      NAPL Hydraulic Gradient  ---
!
          ELSEIF( IBCT(IEQO,NB).EQ.11 ) THEN
            PNX = FNHGBN( NBHG(3,NB),NB,M )
          ENDIF
!
!---      Absolute system pressure  ---
!
          PX = MAX( PLX,PGX,PNX ) + PATM
!
!---      Aqueous density, saturated water-vapor pressure and
!         aqueous viscosity  ---
!
          CALL WATLQD( TX,PX,RHOLB(M,NB) )
          CALL WATSP( TX,PSWB(2,NB) )
          CALL WATLQV( TX,PX,PSWB(2,NB),VISLB(M,NB) )
!
!---      Saturated oil-vapor pressure, aqueous-phase mole and mass
!         fractions  ---
!
          CALL VOCSP( TX,PSOB(2,NB) )
          PVOB(M,NB) = PSOB(2,NB)
!
!---      Compute the aqueous-phase mole and mass fractions  ---
!
          XMLOB(M,NB) = BCX(4)*PVOB(M,NB)/HCOW
          IF( XMLOB(M,NB).LT.EPSL ) XMLOB(M,NB) = 0.D+0
          XMLWB(M,NB) = MAX( 1.D+0-XMLOB(M,NB),0.D+0 )
          WTMX = XMLWB(M,NB)*WTMW + XMLOB(M,NB)*WTMO
          RHOMLB(M,NB) = RHOLB(M,NB)/WTMX
          XLWB(M,NB) = XMLWB(M,NB)*WTMW/WTMX
          XLOB(M,NB) = XMLOB(M,NB)*WTMO/WTMX
          IF( XLOB(M,NB).LT.EPSL ) XLOB(M,NB) = 0.D+0
          IF( XLWB(M,NB).LT.EPSL ) XLWB(M,NB) = 0.D+0
!
!---      Fluid-pair interfacial tension and scaling factor   ---
!
          CALL SFT4( XLOB(M,NB),IZN )
          ENPR = SCHR(2,IZN)*RHORL*GRAV
          PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
          PNX = MAX( PNC,PNX )
!
!---      Compute phase saturation and relative permeability
!         assuming no entrapment on boundary surfaces  ---
!
          INDX = 2
          TRPNLX = SCHR(9,IZN)/(SCHR(9,IZN)+TRPNL(2,N))
          ASNTMX = SCHR(11,IZN)*TRPNLX
          ASNRMX = SCHR(10,IZN)
          CALL KSP4( PGX,PNX,PLX,SLB(M,NB),SNB(M,NB),SGB(M,NB),
     &      RKLB(1,M,NB),RKNB(M,NB),ASTX,ASLX,ASLMINX,
     &      ASTMAXX,ASNTX,ASNRX,ASNTMX,ASNRMX,INDX,IZN,NB )
!
!---      Compute porosity and tortuosity  ---
!
          CALL PORSTY( N,PX,PCMP(N),PORDB(M,NB),PORTB(M,NB))
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SLB(M,NB),SGB(M,NB),
     &      SNB(M,NB),PORDB(M,NB),TORLB(M,NB),TORGX,TORNB(M,NB) )
!
!---      Compute the NAPL density and viscosity  ---
!
          CALL VOCLQD( TX,PX,PSOB(M,NB),RHONB(M,NB) )
          CALL VOCLQV( TX,VISNB(M,NB),VISCO,IVISO )
!
!---      Compute dissolved-oil diffusion coefficients  ---
!
          CALL VOCDFL( TX,VISLB(M,NB),DFLOB(M,NB) )
!
!---      Assign boundary primary variables  ---
!
          TB(M,NB) = TX
          PLB(M,NB) = PLX
          PNB(M,NB) = PNX
          PGB(M,NB) = PGX
  300   CONTINUE
  400 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCP4 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CAP4( ASLMINX,PGX,PLX,PNX,ASNTMX,ASNX,INDX,IZN,NX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Compute the NAPL pressure from the NAPL saturation and
!     aqueous pressure.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1999.
!     Last Modified by MD White on January 11, 1999.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/CAP4'
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
          WRITE(ISC,'(A)') 'Execution Note: CAP4 Convergence Failure'
          WRITE(IWR,'(A)') 'Execution Note: CAP4 Convergence Failure'
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
          CALL UPDT4
          INDX = -INDX
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
        IF( ABS(DPNZ).GT.1.D-6 ) GOTO 150
  151   CONTINUE
        IF( NC.LT.1 .AND. ASNX.GT.1.D-6 ) THEN
          ICNV = 1
          WRITE(ISC,'(A)') 'Execution Note: CAP4 Convergence Failure'
          WRITE(IWR,'(A)') 'Execution Note: CAP4 Convergence Failure'
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
          CALL UPDT4
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
          WRITE(ISC,'(A)') 'Execution Note: CAP4 Convergence Failure'
          WRITE(IWR,'(A)') 'Execution Note: CAP4 Convergence Failure'
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
          CALL UPDT4
          INDX = -INDX
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
        IF( ABS(DPNZ).GT.1.D-6 ) GOTO 160
  161   CONTINUE
        IF( NC.LT.1 .AND. ASNX.GT.1.D-6 ) THEN
          ICNV = 1
          WRITE(ISC,'(A)') 'Execution Note: CAP4 Convergence Failure'
          WRITE(IWR,'(A)') 'Execution Note: CAP4 Convergence Failure'
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
          CALL UPDT4
          INDX = -INDX
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
        PNX = MAX( PNC,PNZ )
  162   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Water-Oil Mode
!
!     Check the thermodynamic and hydrologic states declared through
!     user inputs.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     Last Modified by MD White, PNNL, November 2, 1999.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
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
      SUB_LOG(ISUB_LOG) = '/CHK4'
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
        IF( ISCHR(N).GE.31 .AND. ISCHR(N).LE.34 ) THEN
          IF( SCHR(14,N).LE.ZERO ) THEN
            IF( IRPL(N).EQ.2 ) THEN
              SCHR(14,N) = 1.D+0 - 2.D+0/SCHR(3,N)
            ELSE
              SCHR(14,N) = 1.D+0 - 1.D+0/SCHR(3,N)
            ENDIF
          ENDIF
          IF( SCHR(14,N).LE.ZERO ) THEN
            NCH = INDEX( ROCK(N),'  ' )-1
            CHMSG = 'Negative or Zero van Genuchten m Parameter: '//
     &ROCK(N)(1:NCH)
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
#ifdef ecke
!
!---    Load reactive transport total and diffusive porosity  ---
!
        POR0(1,N) = POR(1,IZ(N))
        POR0(2,N) = POR(2,IZ(N))
#endif
  100 CONTINUE
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 110 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 110
        PX = MAX( PG(2,N),PL(2,N),PN(2,N) ) + PATM
        PX = MAX( PX,PSW(2,N)+PSO(2,N) )
        CALL WATLQD( T(2,N),PX,RHOLWX )
        RHOL(2,N) = RHOLWX
!
!---    Dissolved-oil initial condition provided as concentration  ---
!
        IF( XLO(2,N).LT.-EPSL .AND. XMLO(2,N).LT.-EPSL ) THEN
          RHOLOX = -XLO(2,N)
          RHOL(2,N) = RHOLWX + RHOLOX
          XLO(2,N) = RHOLOX/RHOL(2,N)
          XLW(2,N) = RHOLWX/RHOL(2,N)
          WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO
          XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
          XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---    Dissolved-oil initial condition provided as mass fraction  ---
!
        ELSEIF( XLO(2,N).GT.EPSL .AND. XMLO(2,N).LT.-EPSL ) THEN
          RHOL(2,N) = RHOLWX/(1.D+0-XLO(2,N))
          XLW(2,N) = RHOLWX/RHOL(2,N)
          WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO
          XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
          XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---    Dissolved-oil initial condition provided as mole fraction  ---
!
        ELSEIF( XLO(2,N).LT.-EPSL .AND. XMLO(2,N).GT.EPSL ) THEN
          RHOLOX = (XMLO(2,N)*RHOLWX*WTMO)/((1.D+0-XMLO(2,N))*WTMW)
          RHOL(2,N) = RHOLOX+RHOLWX
          XLO(2,N) = RHOLOX/RHOL(2,N)
          XLW(2,N) = RHOLWX/RHOL(2,N)
          WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO
          XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
          XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
        ELSE
          XLO(2,N) = 0.D+0
          XLW(2,N) = 0.D+0
          XMLO(2,N) = 0.D+0
          XMLW(2,N) = 0.D+0
        ENDIF
        IF( XMLO(2,N).LT.EPSL ) XMLO(2,N) = 0.D+0
        IF( XMLW(2,N).LT.EPSL ) XMLW(2,N) = 0.D+0
        IF( XLO(2,N).LT.EPSL ) XLO(2,N) = 0.D+0
        IF( XLW(2,N).LT.EPSL ) XLW(2,N) = 0.D+0
  110 CONTINUE
!
!---  Check initial dissolved VOC mass fraction  ---
!
      DO 130 N = 1,NFLD
        XMLOX = (PSO(2,N)+1.D+0)/HCOW
        IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.2 ) XMLOX = 1.50D+0*XMLOX
        IF( XMLO(2,N).GT.XMLOX .OR. XMLO(2,N).LT.0.D+0 ) THEN
          CHMSG = 'Dissolved-Oil Mole Fraction: Node: '
          INDX = 16
          IMSG = N
          RLMSG = XMLO(2,N)
          CALL WRMSGS( INDX )
        ENDIF
  130 CONTINUE
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
!
!---    Fluid-pair interfacial tension and scaling factor   ---
!
        CALL SFT4( XLO(2,N),IZN )
        IF( PN(2,N).LE.-1.D+20 ) THEN
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &        (PL(2,N)*BNL+ENPR)/BNL )
          PN(2,N) = PN(2,N) + 1.D+4
        ENDIF
        PX = MAX( PL(2,N),PG(2,N),PN(2,N) ) + PATM
        ENPR = SCHR(2,IZN)*RHORL*GRAV
        HGL = MAX( ZERO,BGL*(PG(2,N)-PL(2,N)) )
        HNL = MAX( ZERO,BNL*(PN(2,N)-PL(2,N)) )
        HGN = MAX( ZERO,BGN*(PG(2,N)-PN(2,N)) )
!
!---  Aqueous and aqueous-gas conditions  ---
!
        IF( HGN.GE.HNL .OR. HNL.LE.ENPR ) THEN
          IF( HGL.GT.ENPR ) THEN
            IF( SCHR(11,IZN).GT.EPSL .AND. SNT(2,N).GT.EPSL ) THEN
              NPHAZ(2,N) = 4
              PVO(2,N) = PSO(2,N)
              PVW(2,N) = PSW(2,N)
            ELSEIF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
              NPHAZ(2,N) = 2
              PVO(2,N) = XMLO(2,N)*HCOW
              PVW(2,N) = PSW(2,N)
            ELSE
              NPHAZ(2,N) = 6
              PVO(2,N) = PSO(2,N)
              PVW(2,N) = PSW(2,N)
            ENDIF
            PN(2,N) = (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN)
          ELSE
            IF( SCHR(11,IZN).GT.EPSL .AND. SNT(2,N).GT.EPSL ) THEN
              NPHAZ(2,N) = 3
              PVO(2,N) = PSO(2,N)
              PVW(2,N) = PSW(2,N)
            ELSEIF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
              NPHAZ(2,N) = 1
              PVO(2,N) = XMLO(2,N)*HCOW
              PVW(2,N) = PSW(2,N)
            ELSE
              NPHAZ(2,N) = 5
              PVO(2,N) = PSO(2,N)
              PVW(2,N) = PSW(2,N)
            ENDIF
            PN(2,N) = (PL(2,N)*BNL + ENPR)/BNL
          ENDIF
!
!---  Aqueous-NAPL and aqueous-NAPL-gas conditions  ---
!
        ELSE
          IF( HGN.GT.ENPR ) THEN
            NPHAZ(2,N) = 6
            PVO(2,N) = PSO(2,N)
            PVW(2,N) = PSW(2,N)
          ELSE
            PVO(2,N) = PSO(2,N)
            PVW(2,N) = PSW(2,N)
            NPHAZ(2,N) = 5
          ENDIF
        ENDIF
!
!---  Call the saturation function to initialize the
!     apparent total-liquid and apparent aqueous saturations  ---
!
        INDX = 2
        ASNT(N) = SNT(2,N)/(1.D+0-SCHR(4,IZN))
        ASNR(N) = SNR(2,N)/(1.D+0-SCHR(4,IZN))
        TRPNLX = SCHR(9,IZN)/(SCHR(9,IZN)+TRPNL(2,N))
        ASNTMX = SCHR(11,IZN)*TRPNLX
        ASNRMX = SCHR(10,IZN)
        CALL KSP4( PG(2,N),PN(2,N),PL(2,N),SLX,SNX,SGX,RKLX,RKNX,
     &    AST(N),ASL(N),ASLMINX,ASTMAXX,ASNT(N),ASNR(N),
     &    ASNTMX,ASNRMX,INDX,IZN,N )
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
      CALL TRPNL4
!
!---  Compute the initial effective aqueous and total-liquid
!     minimum saturations, and check for excessive trapped-NAPL and
!     trapped-gas saturations  ---
!
      DO 220 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 220
        IZN = IZ(N)
!
!---    Fluid-pair interfacial tension and scaling factor   ---
!
        CALL SFT4( XLO(2,N),IZN )
!
!---  Trapped NAPL check  ---
!
        TRPNLX = SCHR(9,IZN)/(SCHR(9,IZN)+TRPNL(2,N))
        ASNTMX = SCHR(11,IZN)*TRPNLX
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
        IF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.4 ) SN(2,N) = SNT(2,N)
!
!---    Load old time step values of
!       apparent minimum aqueous saturation, and
!       apparent maximum total-liquid saturation  ---
!
        ASLMIN(1,N) = ASLMIN(2,N)
        ASTMAX(1,N) = ASTMAX(2,N)
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
#ifdef future_coding
!
!---  Initialize well pressures  ---
!
      DO 260 NS = 1,NSR
        IF( ISRT(NS).GE.20 .AND. ISRT(NS).LE.29 ) THEN
          I1X = ISRDM(1,NS)
          J1X = ISRDM(3,NS)
          K1X = ISRDM(5,NS)
          N = ND(I1X,J1X,K1X)
          DO 250 M = 1,ISVC+2
            PLWB(M,NS) = PL(2,N) + 5.D-1*RHOL(2,N)*GRAVZ*DZGF(NS)
  250     CONTINUE
        ENDIF
  260 CONTINUE
#endif
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CISC4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Compute initial solute concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, May 1994.
!     Last Modified by MD White, Battelle, June 30, 1994.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/CISC4'
      DO 140 NSL = 1,NSOLU
        DO 110 N = 1,NFLD
          IF( IXP(N).EQ.0 ) GOTO 110
          IZN = IZ(N)
          SVLX = PORD(2,N)*SL(2,N)+SMALL
          SVNX = PORD(2,N)*SN(2,N)+SMALL
!
!-------  Linear Kd sorption  ---
!
          IF( IPCSL(IZN,NSL).EQ.1 ) THEN
            PCSLX = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
!
!-------  Linear, Freundlich, or Langmuir sorption  ---
!
          ELSE
            PCSLX = PCSL(1,IZN,NSL)
          ENDIF
!
!-------  Equilibrium initial conditions:
!         declared aqueous solute concentration  ---
!
          IF( ICT(N,NSL).EQ.2 ) THEN
!
!---------  Aqueous-NAPL adsorption linear isotherm and
!           Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
            IF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*(SVLX + PCSLX)
              ELSEIF( PCSLX.LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*(SVLX + SVNX*PCLN(1,NSL))
              ELSE
                C(N,NSL) = C(N,NSL)*(SVLX + SVNX*PCLN(1,NSL) +
     &            PCSLX)
              ENDIF
!
!---------  Aqueous-NAPL adsorption Freundlich isotherm and
!           Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*(SVLX + PCSLX)
              ELSEIF( PCSLX.LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*SVLX + C(N,NSL)**PCLN(2,NSL)*
     &            SVNX*PCLN(1,NSL)
              ELSE
                C(N,NSL) = C(N,NSL)*SVLX + C(N,NSL)**PCLN(2,NSL)*
     &            SVNX*PCLN(1,NSL) + C(N,NSL)*PCSLX
              ENDIF
!
!---------  Aqueous-NAPL adsorption Langmuir isotherm and
!           Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*(SVLX + PCSLX)
              ELSEIF( PCSLX.LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*(SVLX + (SVNX*PCLN(1,NSL))/
     &            (1.D+0 + C(N,NSL)*PCLN(2,NSL)))
              ELSE
                C(N,NSL) = C(N,NSL)*(SVLX + (SVNX*PCLN(1,NSL))/
     &            (1.D+0 + C(N,NSL)*PCLN(2,NSL)) + PCSLX)
              ENDIF
!
!---------  Aqueous-NAPL adsorption linear isotherm and
!           Solid-Aqueous adsorption Freundlich isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*SVLX + C(N,NSL)**PCSL(2,IZN,NSL)*
     &            PCSL(1,IZN,NSL)
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*SVLX + C(N,NSL)*SVNX*PCLN(1,NSL)
              ELSE
                C(N,NSL) = C(N,NSL)*SVLX + C(N,NSL)*SVNX*PCLN(1,NSL) +
     &            C(N,NSL)**PCSL(2,IZN,NSL)*PCSL(1,IZN,NSL)
              ENDIF
!
!---------  Aqueous-NAPL adsorption Freundlich isotherm and
!           Solid-Aqueous adsorption Freundlich isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*SVLX + C(N,NSL)**PCSL(2,IZN,NSL)*
     &            PCSL(1,IZN,NSL)
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*SVLX + C(N,NSL)**PCLN(2,NSL)*
     &            SVNX*PCLN(1,NSL)
              ELSE
                C(N,NSL) = C(N,NSL)*SVLX + C(N,NSL)**PCLN(2,NSL)*
     &            SVNX*PCLN(1,NSL) + C(N,NSL)**PCSL(2,IZN,NSL)*
     &            PCSL(1,IZN,NSL)
              ENDIF
!
!---------  Aqueous-NAPL adsorption Langmuir isotherm and
!           Solid-Aqueous adsorption Fruendlich isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*SVLX + C(N,NSL)**PCSL(2,IZN,NSL)*
     &            PCSL(1,IZN,NSL)
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*SVLX + (C(N,NSL)*SVNX*PCLN(1,NSL))/
     &            (1.D+0 + C(N,NSL)*PCLN(2,NSL))
              ELSE
                C(N,NSL) = C(N,NSL)*SVLX + (C(N,NSL)*SVNX*PCLN(1,NSL))/
     &            (1.D+0 + C(N,NSL)*PCLN(2,NSL)) +
     &            C(N,NSL)**PCSL(2,IZN,NSL)*PCSL(1,IZN,NSL)
              ENDIF
!
!---------  Aqueous-NAPL adsorption linear isotherm and
!           Solid-Aqueous adsorption Langmuir isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*(SVLX + PCSL(1,IZN,NSL)/
     &            (1.D+0 + C(N,NSL)*PCSL(2,IZN,NSL)))
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*(SVLX + SVNX*PCLN(1,NSL))
              ELSE
                C(N,NSL) = C(N,NSL)*(SVLX + SVNX*PCLN(1,NSL) +
     &            PCSL(1,IZN,NSL)/(1.D+0 + C(N,NSL)*PCSL(2,IZN,NSL)))
              ENDIF
!
!---------  Aqueous-NAPL adsorption Freundlich isotherm and
!           Solid-Aqueous adsorption Langmuir isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*SVLX + (C(N,NSL)*PCSL(1,IZN,NSL))/
     &            (1.D+0 + C(N,NSL)*PCSL(2,IZN,NSL))
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*SVLX + C(N,NSL)**PCLN(2,NSL)
     &            *SVNX*PCLN(1,NSL)
              ELSE
                C(N,NSL) = C(N,NSL)*SVLX + C(N,NSL)**PCLN(2,NSL)
     &            *SVNX*PCLN(1,NSL) + (C(N,NSL)*PCSL(1,IZN,NSL))/
     &            (1.D+0 + C(N,NSL)*PCSL(2,IZN,NSL))
              ENDIF
!
!---------  Aqueous-NAPL adsorption Langmuir isotherm and
!           Solid-Aqueous adsorption Langmuir isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*(SVLX + PCSL(1,IZN,NSL)/
     &            (1.D+0 + C(N,NSL)*PCSL(2,IZN,NSL)))
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                C(N,NSL) = C(N,NSL)*(SVLX + (SVNX*PCLN(1,NSL))/
     &            (1.D+0 + C(N,NSL)*PCLN(2,NSL)))
              ELSE
                C(N,NSL) = C(N,NSL)*(SVLX + (SVNX*PCLN(1,NSL))/
     &            (1.D+0 + C(N,NSL)*PCLN(2,NSL)) + PCSL(1,IZN,NSL)/
     &            (1.D+0 + C(N,NSL)*PCSL(2,IZN,NSL)))
              ENDIF
            ENDIF
!
!-------  Equilibrium initial conditions:
!         declared NAPL solute concentration  ---
!
          ELSEIF( ICT(N,NSL).EQ.3 ) THEN
!
!---------  Aqueous-NAPL adsorption linear isotherm and
!           Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
            IF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                INDX = 17
                N_DB = N
                CHMSG = 'Infinite Solute Concentration: '//
     &            'Declared NAPL Solute Concentration: @ Node'
                CALL WRMSGS(INDX)
              ELSEIF( PCSLX.LE.EPSL ) THEN
                CLX = C(N,NSL)/PCLN(1,NSL)
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX
              ELSE
                CLX = C(N,NSL)/PCLN(1,NSL)
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX + CLX*PCSLX
              ENDIF
!
!---------  Aqueous-NAPL adsorption Freundlich isotherm and
!           Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                INDX = 17
                N_DB = N
                CHMSG = 'Infinite Solute Concentration: '//
     &            'Declared NAPL Solute Concentration: @ Node'
                CALL WRMSGS(INDX)
              ELSEIF( PCSLX.LE.EPSL ) THEN
                CLX = (C(N,NSL)/PCLN(1,NSL))**(1/PCLN(2,NSL))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX
              ELSE
                CLX = (C(N,NSL)/PCLN(1,NSL))**(1/PCLN(2,NSL))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX + CLX*PCSLX
              ENDIF
!
!---------  Aqueous-NAPL adsorption Langmuir isotherm and
!           Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                INDX = 17
                N_DB = N
                CHMSG = 'Infinite Solute Concentration: '//
     &            'Declared NAPL Solute Concentration: @ Node'
                CALL WRMSGS(INDX)
              ELSEIF( PCSLX.LE.EPSL ) THEN
                CLX = -(C(N,NSL)/(-PCLN(1,NSL) + C(N,NSL)*PCLN(2,NSL)))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX
              ELSE
                CLX = -(C(N,NSL)/(-PCLN(1,NSL) + C(N,NSL)*PCLN(2,NSL)))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX + CLX*PCSLX
              ENDIF
!
!---------  Aqueous-NAPL adsorption linear isotherm and
!           Solid-Aqueous adsorption Freundlich isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                INDX = 17
                N_DB = N
                CHMSG = 'Infinite Solute Concentration: '//
     &            'Declared NAPL Solute Concentration: @ Node'
                CALL WRMSGS(INDX)
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                CLX = C(N,NSL)/PCLN(1,NSL)
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX
              ELSE
                CLX = C(N,NSL)/PCLN(1,NSL)
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX +
     &            CLX**PCSL(2,IZN,NSL)*PCSL(1,IZN,NSL)
              ENDIF
!
!---------  Aqueous-NAPL adsorption Freundlich isotherm and
!           Solid-Aqueous adsorption Freundlich isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                INDX = 17
                N_DB = N
                CHMSG = 'Infinite Solute Concentration: '//
     &            'Declared NAPL Solute Concentration: @ Node'
                CALL WRMSGS(INDX)
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                CLX = (C(N,NSL)/PCLN(1,NSL))**(1/PCLN(2,NSL))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX
              ELSE
                CLX = (C(N,NSL)/PCLN(1,NSL))**(1/PCLN(2,NSL))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX +
     &            CLX**PCSL(2,IZN,NSL)*PCSL(1,IZN,NSL)
              ENDIF
!
!---------  Aqueous-NAPL adsorption Langmuir isotherm and
!           Solid-Aqueous adsorption Fruendlich isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                INDX = 17
                N_DB = N
                CHMSG = 'Infinite Solute Concentration: '//
     &            'Declared NAPL Solute Concentration: @ Node'
                CALL WRMSGS(INDX)
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                CLX = -(C(N,NSL)/(-PCLN(1,NSL) + C(N,NSL)*PCLN(2,NSL)))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX
              ELSE
                CLX = -(C(N,NSL)/(-PCLN(1,NSL) + C(N,NSL)*PCLN(2,NSL)))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX +
     &            CLX**PCSL(2,IZN,NSL)*PCSL(1,IZN,NSL)
              ENDIF
!
!---------  Aqueous-NAPL adsorption linear isotherm and
!           Solid-Aqueous adsorption Langmuir isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                INDX = 17
                N_DB = N
                CHMSG = 'Infinite Solute Concentration: '//
     &            'Declared NAPL Solute Concentration: @ Node'
                CALL WRMSGS(INDX)
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                CLX = C(N,NSL)/PCLN(1,NSL)
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX
              ELSE
                CLX = C(N,NSL)/PCLN(1,NSL)
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX +
     &            (CLX*PCSL(1,IZN,NSL))/(1.D+0 + CLX*PCSL(2,IZN,NSL))
              ENDIF
!
!---------  Aqueous-NAPL adsorption Freundlich isotherm and
!           Solid-Aqueous adsorption Langmuir isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                INDX = 17
                N_DB = N
                CHMSG = 'Infinite Solute Concentration: '//
     &            'Declared NAPL Solute Concentration: @ Node'
                CALL WRMSGS(INDX)
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                CLX = (C(N,NSL)/PCLN(1,NSL))**(1/PCLN(2,NSL))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX
              ELSE
                CLX = (C(N,NSL)/PCLN(1,NSL))**(1/PCLN(2,NSL))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX +
     &            (CLX*PCSL(1,IZN,NSL))/(1.D+0 + CLX*PCSL(2,IZN,NSL))
              ENDIF
!
!---------  Aqueous-NAPL adsorption Langmuir isotherm and
!           Solid-Aqueous adsorption Langmuir isotherm  ---
!
            ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
              IF( SVNX*PCLN(1,NSL).LE.EPSL ) THEN
                INDX = 17
                N_DB = N
                CHMSG = 'Infinite Solute Concentration: '//
     &            'Declared NAPL Solute Concentration: @ Node'
                CALL WRMSGS(INDX)
              ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
                CLX = -(C(N,NSL)/(-PCLN(1,NSL) + C(N,NSL)*PCLN(2,NSL)))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX
              ELSE
                CLX = -(C(N,NSL)/(-PCLN(1,NSL) + C(N,NSL)*PCLN(2,NSL)))
                C(N,NSL) = CLX*SVLX + C(N,NSL)*SVNX +
     &            (CLX*PCSL(1,IZN,NSL))/(1.D+0 + CLX*PCSL(2,IZN,NSL))
              ENDIF
            ENDIF
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CISC4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Compute primary variable increments
!
!     Phase Conditions
!     Number  Description                              Primary Unknowns
!       1     Aqueous w/ Dissolved Oil                 PL,XMLO
!       2     Aqueous-Gas w/ Dissolved Oil             PL,XMLO
!       3     Aqueous w/ Trapped-Residual NAPL         PL,SN
!       4     Aqueous-Gas w/ Trapped-Residual NAPL     PL,SN
!       5     Aqueous-NAPL                             PL,PN
!       6     Aqueous-NAPL-Gas                         PL,PN
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White and M Oostrom, Battelle, May 9, 1994.
!     Last Modified by MD White, PNNL, 5 September 2002.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!      DATA TOL1,TOL2 /1.D-6,1.D-4/
      DATA TOL1,TOL2 /0.D+0,0.D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM4'
!
!---  Phase options, compute phase condition   ---
!
      DO 100 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 100
        IZN = IZ(N)
!
!---    Fluid-pair interfacial tension and scaling factor   ---
!
        XMLOX = MIN( XMLO(2,N),PSO(2,N)/HCOW )
        XMLWX = MAX( 1.D+0-XMLOX,0.D+0 )
        XLOX = XMLOX*WTMO/(XMLOX*WTMO + XMLWX*WTMW)
        CALL SFT4( XLOX,IZN )
        ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---    Capillary pressures   ---
!
        HGL = MAX( ZERO,BGL*(PG(2,N)-PL(2,N)) )
        HNL = MAX( ZERO,BNL*(PN(2,N)-PL(2,N)) )
        HGN = MAX( ZERO,BGN*(PG(2,N)-PN(2,N)) )
!
!---    Aqueous or aqueous-gas w/ dissolved-oil conditions  ---
!
        IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.2 ) THEN
          IF( HGL.GT.ENPR ) THEN
            PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &        (PL(2,N)*BNL+ENPR)/BNL )
            IF( XMLO(2,N).GT.VSLC(1)*(PSO(2,N)/HCOW) ) THEN
              NPHAZ(2,N) = 6
            ELSEIF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
              NPHAZ(2,N) = 2
            ELSE
              NPHAZ(2,N) = 6
            ENDIF
          ELSE
            PN(2,N) = (PL(2,N)*BNL + ENPR)/BNL
            IF( XMLO(2,N).GT.VSLC(1)*(PSO(2,N)/HCOW) ) THEN
              NPHAZ(2,N) = 5
            ELSEIF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
              NPHAZ(2,N) = 1
            ELSE
              NPHAZ(2,N) = 5
            ENDIF
          ENDIF
!
!---    Aqueous or aqueous-gas w/ trapped-oil conditions  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.4 ) THEN
          TRPNLX = SCHR(9,IZN)/(SCHR(9,IZN)+TRPNL(2,N))
          ASNTMX = SCHR(11,IZN)*TRPNLX
          IF( ASNTMX.GE.EPSL ) THEN
            INDX = 4
            ASLMINX = ASLMIN(1,N)
            ASTMAXX = ASTMAX(1,N)
            CALL KSP4( PG(2,N),PN(2,N),PL(2,N),SLX,SN(2,N),SGX,
     &        RKLX,RKNX,ASTX,ASLX,ASLMINX,ASTMAXX,ASNTX,
     &        ASNRX,ASNTMX,ASNRMX,INDX,IZN,N )
          ENDIF
          ASNX = SN(2,N)/(1.D+0-SCHR(4,IZN))
          IF( ASTX.LT.1.D+0 ) THEN
            IF( (ASNX-ASNTX)/(ASNTMX+EPSL).GT.1.D-1 ) THEN
              NPHAZ(2,N) = 6
            ELSEIF( ASNX.LT.EPSL ) THEN
              IF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
                NPHAZ(2,N) = 2
              ELSE
                NPHAZ(2,N) = 4
              ENDIF
            ELSE
              NPHAZ(2,N) = 4
            ENDIF
          ELSE
            IF( (ASNX-ASNTX)/(ASNTMX+EPSL).GT.1.D-1 ) THEN
              NPHAZ(2,N) = 5
            ELSEIF( ASNX.LT.EPSL ) THEN
              IF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
                NPHAZ(2,N) = 1
              ELSE
                NPHAZ(2,N) = 3
              ENDIF
            ELSE
              NPHAZ(2,N) = 3
            ENDIF
          ENDIF
!
!---    Aqueous-NAPL or aqueous-NAPL-gas conditions  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 .OR. NPHAZ(2,N).EQ.6 ) THEN
          TRPNLX = SCHR(9,IZN)/(SCHR(9,IZN)+TRPNL(2,N))
          ASNTMX = SCHR(11,IZN)*TRPNLX
          IF( ASNTMX.GE.EPSL ) THEN
           INDX = 1
           ASLMINX = ASLMIN(1,N)
           ASTMAXX = ASTMAX(1,N)
           CALL KSP4( PG(2,N),PN(2,N),PL(2,N),SLX,SNX,SGX,
     &       RKLX,RKNX,ASTX,ASLX,ASLMINX,ASTMAXX,ASNTX,
     &       ASNRX,ASNTMX,ASNRMX,INDX,IZN,N )
            ASNX = SNX/(1.D+0-SCHR(4,IZN))
            IF( ASNX.GT.EPSL ) THEN
              RFNX = MIN( (ASNX-ASNTX)/(ASNTX+EPSL),1.D+0 )
            ELSE
              RFNX = 1.D+0
            ENDIF
          ELSE
            RFNX = 1.D+0
          ENDIF
          IF( HGN.GE.(HNL-TOL1) .OR. (HNL-TOL1).LE.ENPR ) THEN
            PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &        (PL(2,N)*BNL+ENPR)/BNL )
            IF( HGL.GT.ENPR ) THEN
              IF( SCHR(11,IZN).GT.EPSL .AND. SNT(2,N).GT.EPSL ) THEN
                NPHAZ(2,N) = 4
              ELSEIF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
                NPHAZ(2,N) = 2
              ELSE
                NPHAZ(2,N) = 6
              ENDIF
            ELSE
              IF( SCHR(11,IZN).GT.EPSL .AND. SNT(2,N).GT.EPSL ) THEN
                NPHAZ(2,N) = 3
              ELSEIF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
                NPHAZ(2,N) = 1
              ELSE
                NPHAZ(2,N) = 5
              ENDIF
            ENDIF
          ELSE
            PNCX = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &        (PL(2,N)*BNL+ENPR)/BNL )
            IF( PN(2,N)-PNCX.LT.TOL2 ) PN(2,N) = PNCX
            IF( HGN.GT.ENPR ) THEN
              IF( RFNX.LT.1.D-2 ) THEN
                NPHAZ(2,N) = 4
              ELSE
                NPHAZ(2,N) = 6
              ENDIF
            ELSE
              IF( RFNX.LT.1.D-2 ) THEN
                NPHAZ(2,N) = 3
              ELSE
                NPHAZ(2,N) = 5
              ENDIF
            ENDIF
          ENDIF
        ENDIF
  100 CONTINUE
!
!---  Load previous time step condition for initial time step  ---
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
!
!---    Fluid-pair interfacial tension and scaling factor   ---
!
        XMLOX = MIN( XMLO(2,N),PSO(2,N)/HCOW )
        XMLWX = MAX( 1.D+0-XMLOX,0.D+0 )
        XLOX = XMLOX*WTMO/(XMLOX*WTMO + XMLWX*WTMW)
        CALL SFT4( XLOX,IZN )
        ENPR = SCHR(2,IZN)*RHORL*GRAV
        XMLOX = PSO(2,N)/HCOW
!
!---  Aqueous or aqueous-gas w/ dissolved-oil conditions  ---
!
        IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.2 ) THEN
          IF( ISCHR(IZN).GE.31 .AND. ISCHR(IZN).LE.34 ) THEN
            HDX = MAX( (PG(2,N)-PL(2,N))/RHORL/GRAV,0.D+0 )
            ASLX = -SCHR(14,IZN)*SCHR(3,IZN)*SCHR(1,IZN)*
     &        ((HDX*SCHR(1,IZN))**(SCHR(3,IZN)-1.D+0))*
     &        ((1.D+0 + (HDX*SCHR(1,IZN))**SCHR(3,IZN))**
     &         (-1.D+0 - SCHR(14,IZN)))
            DPX = 1.D-1*EXP(SCHR(3,IZN))/(EXP(ABS(ASLX)))
          ELSE
            DPX = 1.D-1
          ENDIF
!!
!!---      Node adjacent to surface spill node  ---
!!
!          IF( KD(N).EQ.KSPS(ID(N),JD(N)) ) THEN
!            DNR(IEQW,N) = -1.D-4
!          ELSE
            DNR(IEQW,N) = -MAX( DPX,1.D-6*ABS(PG(2,N)-PL(2,N)) )
!          ENDIF
          DNR(IEQO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
!
!---    Aqueous or aqueous-gas w/ trapped-oil conditions  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.4 ) THEN
          PNCX = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      (PL(2,N)*BNL+ENPR)/BNL )
          DNR(IEQW,N) = -MAX( 1.D-3,1.D-5*ABS(PN(2,N)-PNCX) )
          DNRO = MAX(1.D-9,1.D-4*SN(2,N))
          DNR(IEQO,N) = SIGN( DNRO,5.D-1-SN(2,N) )
!
!---    Aqueous-NAPL conditions  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
          DNR(IEQW,N) = -MAX( 1.D-3,1.D-5*ABS(PN(2,N)-PL(2,N))  )
          DNR(IEQO,N) = MAX( 1.D-3,1.D-5*ABS(PN(2,N)-PL(2,N))  )
          IF( SN(2,N)/(1.D+0-SCHR(4,IZN)).GT.5.D-1 ) THEN
            DNR(IEQW,N) = -DNR(IEQW,N)
            DNR(IEQO,N) = -DNR(IEQO,N)
          ENDIF
!
!---    Aqueous-NAPL-gas conditions  ---
!
        ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
          DNR(IEQW,N) = -MAX( 1.D-3,1.D-5*ABS(PN(2,N)-PL(2,N)) )
          DNR(IEQO,N) = MAX( 1.D-3,1.D-5*ABS(PG(2,N)-PN(2,N)) )
        ENDIF
!
!---    Increment the primary variables  ---
!
        DO 200 M = 2,ISVC+2
          PL(M,N) = PL(2,N)
          PN(M,N) = PN(2,N)
          PG(M,N) = PG(2,N)
          XMLO(M,N) = XMLO(2,N)
          SN(M,N) = SN(2,N)
!
!---      Aqueous or aqueous-gas w/ dissolved-oil conditions  ---
!
          IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.2 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQO+2 ) THEN
              XMLO(M,N) = XMLO(M,N) + DNR(IEQO,N)
              IF( XMLO(M,N).LT.EPSL ) XMLO(M,N) = 0.D+0
            ENDIF
!
!---        Fluid-pair interfacial tension and scaling factor   ---
!
            XMLOX = MIN( XMLO(M,N),PSO(2,N)/HCOW )
            XMLWX = MAX( 1.D+0-XMLOX,0.D+0 )
            XLOX = XMLOX*WTMO/(XMLOX*WTMO + XMLWX*WTMW)
            CALL SFT4( XLOX,IZN )
            ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---        Critical NAPL pressure   ---
!
            PN(M,N) = MAX( (BNL*PL(M,N)+BGN*PG(M,N))/(BNL+BGN),
     &        (PL(M,N)*BNL+ENPR)/BNL )
!
!---      Aqueous or aqueous-gas w/ trapped-oil conditions  ---
!
          ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.4 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQO+2 ) THEN
              SN(M,N) = SN(M,N) + DNR(IEQO,N)
            ENDIF
!
!---        Fluid-pair interfacial tension and scaling factor   ---
!
            XMLOX = PSO(2,N)/HCOW
            XMLWX = MAX( 1.D+0-XMLOX,0.D+0 )
            XLOX = XMLOX*WTMO/(XMLOX*WTMO + XMLWX*WTMW)
            CALL SFT4( XLOX,IZN )
            ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---        Critical NAPL pressure   ---
!
            PN(M,N) = MAX( (BNL*PL(M,N)+BGN*PG(M,N))/(BNL+BGN),
     &        (PL(M,N)*BNL+ENPR)/BNL )
!
!---      Aqueous-NAPL or aqueous-NAPL-gas conditions  ---
!
          ELSEIF( NPHAZ(2,N).EQ.5 .OR. NPHAZ(2,N).EQ.6 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQO+2 ) THEN
              PN(M,N) = PN(M,N) + DNR(IEQO,N)
            ENDIF
          ENDIF
  200   CONTINUE
  300 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM4 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KSP4( PGX,PNX,PLX,SLX,SNX,SGX,RKLX,RKNX,
     &                 ASTX,ASLX,ASLMINX,ASTMAXX,ASNTX,ASNRX,
     &                 ASNTMX,ASNRMX,INDX,IZN,NX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
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
!     Last Modified by MD White on October 27, 1997.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/KSP4'
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
!---    Residual NAPL Function  ---
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
!---    No Residual NAPL Function ---
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
!---    Free Corey relative permeability function  ---
!
        ELSEIF( IRPL(IZN).EQ.7 ) THEN
          RKLX(1) = RPLC(1,IZN)*(ASLX**(RPLC(2,IZN)))
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
!---    Free Corey relative permeability function  ---
!
        ELSEIF( IRPL(IZN).EQ.7 ) THEN
          RKNX = RPNC(1,IZN)*(ASNX**(RPNC(2,IZN)))
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
!---    Residual NAPL Function  ---
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
!---    Free Corey relative permeability function  ---
!
        ELSEIF( IRPL(IZN).EQ.7 ) THEN
          RKLX(1) = RPLC(1,IZN)*(ASLX**(RPLC(2,IZN)))
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
!---    Free Corey relative permeability function  ---
!
        ELSEIF( IRPL(IZN).EQ.7 ) THEN
          RKNX = RPNC(1,IZN)*(ASNX**(RPNC(2,IZN)))
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
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
!
!---    Dissolved oil condition, assume no trapped NAPL,
!       no residual NAPL, no free NAPL, and
!       zero NAPL relative permeability  ---
!
        ELSEIF( INDX.EQ.3 ) THEN
          ASNTX = 0.D+0
          ASNFX = 0.D+0
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
          ASLMINX = ASLX
          ASTMAXX = ASTX
!
!---    Entrapped NAPL phase condition, compute an estimate
!       of the NAPL pressure, assume no residual NAPL,
!       and zero NAPL relative permeability  ---
!
        ELSEIF( INDX.EQ.4 ) THEN
          ASNX = SNX/(1.D+0-SLRX)
          CALL CAP4( ASLMINX,PGX,PLX,PNX,ASNTMX,ASNX,INDX,IZN,NX )
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
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
!
!---    Full entrapment  ---
!
        ELSE
          IF( (ASLX-ASLMINX).LT.TOL .OR. ASNTMX.LT.TOL ) THEN
            ASNTX = 0.D+0
          ELSE
            PLND = 1.D+0/ASNTMX - 1.D+0
            ASNTX = (1.D+0-ASLMINX)/(1.D+0 + PLND*(1.D+0-ASLMINX))
            ASNTX = ASNTX-(1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ENDIF
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
          ASNFX = ASTX-ASLX
          IF( ASNFX.LT.TOL ) ASNFX = 0.D+0
        ENDIF
        ESLX = ASLX-ASNTX
        IF( ESLX.LT.TOL ) ESLX = 0.D+0
        SLX = (ESLX)*(1.D+0-SLRX) + SLRX
        SNX = (ASNFX+ASNTX)*(1.D+0-SLRX)
        SGX = ASGFX*(1.D+0-SLRX)
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
!---    Residual NAPL  ---
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
!---    Free Corey relative permeability function  ---
!
        ELSEIF( IRPL(IZN).EQ.7 ) THEN
          RKLX(1) = RPLC(1,IZN)*(ASLX**(RPLC(2,IZN)))
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
!---    Free Corey relative permeability function  ---
!
        ELSEIF( IRPL(IZN).EQ.7 ) THEN
          RKNX = RPNC(1,IZN)*(ASNX**(RPNC(2,IZN)))
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
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
!
!---    Dissolved oil condition, assume no trapped NAPL,
!       no residual NAPL, no free NAPL, and
!       zero NAPL relative permeability  ---
!
        ELSEIF( INDX.EQ.3 ) THEN
          ASNTX = 0.D+0
          ASNFX = 0.D+0
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
          ASLMINX = ASLX
!
!---    Entrapped NAPL phase condition, compute an estimate
!       of the NAPL pressure, assume no residual NAPL,
!       and zero NAPL relative permeability  ---
!
        ELSEIF( INDX.EQ.4 ) THEN
          ASNX = SNX/(1.D+0-SLRX)
          CALL CAP4( ASLMINX,PGX,PLX,PNX,ASNTMX,ASNX,INDX,IZN,NX )
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
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
!
!---    Full entrapment  ---
!
        ELSE
          IF( (ASLX-ASLMINX).LT.TOL .OR. ASNTMX.LT.TOL ) THEN
            ASNTX = 0.D+0
          ELSE
            PLND = 1.D+0/ASNTMX - 1.D+0
            ASNTX = (1.D+0-ASLMINX)/(1.D+0 + PLND*(1.D+0-ASLMINX))
            ASNTX = ASNTX-(1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ENDIF
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
          ASNFX = ASTX-ASLX
          IF( ASNFX.LT.TOL ) ASNFX = 0.D+0
        ENDIF
        ESLX = ASLX-ASNTX
        SLX = (ESLX)*(1.D+0-SLRX) + SLRX
        SNX = (ASNFX+ASNTX)*(1.D+0-SLRX)
        SGX = ASGFX*(1.D+0-SLRX)
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
!---    Free Corey relative permeability function  ---
!
        ELSEIF( IRPL(IZN).EQ.7 ) THEN
          RKLX(1) = RPLC(1,IZN)*(ASLX**(RPLC(2,IZN)))
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
!---    Free Corey relative permeability function  ---
!
        ELSEIF( IRPL(IZN).EQ.7 ) THEN
          RKNX = RPNC(1,IZN)*(ASNX**(RPNC(2,IZN)))
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
      RKLX(1) = MAX( RKLX(1),1.D-6 )
      RKLX(2) = MAX( RKLX(2),1.D-6 )
      RKLX(3) = MAX( RKLX(3),1.D-6 )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KSP4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Water-Oil Mode
!
!     Load the current time step values into the old time step
!     variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SPILL
      USE SOURC
      USE SOLTN
      USE REACT
      USE HYST
      USE GRID
      USE FDVP
      USE FDVN
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
      SUB_LOG(ISUB_LOG) = '/LDO4'
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
        SNT(1,N) = SNT(2,N)
        SNR(1,N) = SNR(2,N)
        SG(1,N) = SG(2,N)
        PL(1,N) = PL(2,N)
        PN(1,N) = PN(2,N)
        PVW(1,N) = PVW(2,N)
        PVO(1,N) = PVO(2,N)
        XLW(1,N) = XLW(2,N)
        XLO(1,N) = XLO(2,N)
        XSO(1,N) = XSO(2,N)
        XMLW(1,N) = XMLW(2,N)
        XMLO(1,N) = XMLO(2,N)
        RHOL(1,N) = RHOL(2,N)
        RHOML(1,N) = RHOML(2,N)
        RHON(1,N) = RHON(2,N)
        VISL(1,N) = VISL(2,N)
        VISN(1,N) = VISN(2,N)
        RKL(1,1,N) = RKL(1,2,N)
        RKL(2,1,N) = RKL(2,2,N)
        RKL(3,1,N) = RKL(3,2,N)
        RKN(1,N) = RKN(2,N)
        TORL(1,N) = TORL(2,N)
        TORN(1,N) = TORN(2,N)
        NPHAZ(1,N) = NPHAZ(2,N)
        TRPNL(1,N) = TRPNL(2,N)
        IF( SN(2,N).GT.EPSL ) ASLMIN(1,N) = MIN( ASL(N),ASLMIN(2,N) )
        ASTMAX(1,N) = MAX( AST(N),ASTMAX(2,N) )
        ASLMIN(2,N) = ASLMIN(1,N)
        ASTMAX(2,N) = ASTMAX(1,N)
        DO 90 NSL = 1,NSOLU
          CO(N,NSL) = C(N,NSL)
   90   CONTINUE
#ifdef ecke
        DO 92 NEQ = 1,NEQC+NEQK
          NSL = NEQ + NSOLU 
          CO(N,NSL) = C(N,NSL)
   92   CONTINUE
        DO 94 NSP = 1,NSPR
          SP_CO(N,NSP) = SP_C(N,NSP)
   94   CONTINUE
#endif
  100 CONTINUE
#ifdef awtap
      DO 110 NS = 1,NSR
        QLW(1,NS) = QLW(1,NS) + QLW(3,NS)*DT
        QNW(1,NS) = QNW(1,NS) + QNW(3,NS)*DT
        QTW(1,NS) = QTW(1,NS) + QTW(3,NS)*DT
  110 CONTINUE
!
!---  Initialization for well type sources  ---
!
      DO 120 NS = 1,NSR
        PLWB(1,NS) = PLWB(2,NS)
        PGW(1,NS) = PGW(2,NS)
  120 CONTINUE
#endif
!
!---  Assign old time step values
!     for surface spill variables  ---
!
      IF( ISLC(49).EQ.1 ) THEN
        DO 130 N = 1,IJFLD
          HNSP(1,N) = HNSP(2,N)
          HLSP(1,N) = HLSP(2,N)
          RHONSP(1,N) = RHONSP(2,N)
          RHOLSP(1,N) = RHOLSP(2,N)
          VISNSP(1,N) = VISNSP(2,N)
          VISLSP(1,N) = VISLSP(2,N)
  130   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO4 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP4( INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Control saturation, relative permeability, porosity, and
!     tortuosity calculations
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/PROP4'
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 200
        IZN = IZ(N)
        N_DB = N
#ifdef ecke
        POR0(1,N) = POR0(1,N)
        POR0(2,N) = POR0(2,N)
#endif
        NC = 0
   10   CONTINUE
        DO 100 M = ISVC+2,2,-1
!
!---      Fluid-pair interfacial tension and scaling factor   ---
!
          XMLOX = MIN( XMLO(M,N),PSO(2,N)/HCOW )
          XMLWX = MAX( 1.D+0-XMLOX,0.D+0 )
          XLOX = XMLOX*WTMO/(XMLOX*WTMO + XMLWX*WTMW)
          CALL SFT4( XLOX,IZN )
!
!---      Maximum effective trapped-NAPL saturation  ---
!
          TRPNLX = SCHR(9,IZN)/(SCHR(9,IZN)+TRPNL(2,N))
          ASNTMX = SCHR(11,IZN)*TRPNLX
!
!---      Maximum effective residual-NAPL saturation  ---
!
          ASNRMX = SCHR(10,IZN)
!
!---      For phase conditions with dissolved oil compute saturations
!         assuming no fluid entrapment  ---
!
          INDX = 1
          IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.2 ) INDX = 3
          IF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.4 ) INDX = 4
          ASLMINX = ASLMIN(1,N)
          ASTMAXX = ASTMAX(1,N)
          CALL KSP4( PG(M,N),PN(M,N),PL(M,N),SL(M,N),SN(M,N),SG(M,N),
     &      RKL(1,M,N),RKN(M,N),ASTX,ASLX,ASLMINX,ASTMAXX,
     &      ASNTX,ASNRX,ASNTMX,ASNRMX,INDX,IZN,N )
          IF( INDX.LT.0 ) THEN
            ISUB_LOG = ISUB_LOG-1
            RETURN
          ENDIF
          IF( M.EQ.2 ) THEN
            AST(N) = ASTX
            ASL(N) = ASLX
            ASL(N) = BGL
            ASNT(N) = ASNTX
            ASNR(N) = ASNRX
            ASLMIN(2,N) = ASLMINX
            ASTMAX(2,N) = ASTMAXX
          ENDIF
          SNT(M,N) = ASNTX*(1.D+0-SCHR(4,IZN))
          SNR(M,N) = ASNRX*(1.D+0-SCHR(4,IZN))
!
!---  Compute porosity and tortuosity  ---
!
          PX = MAX( PG(M,N),PL(M,N),PN(M,N) )+PATM
          CALL PORSTY( N,PX,PCMP(N),PORD(M,N),PORT(M,N))
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL(M,N),SG(M,N),
     &      SN(M,N),PORD(M,N),TORL(M,N),TORGX,TORN(M,N) )
  100   CONTINUE
        IF( NPHAZ(2,N).EQ.5 .OR. NPHAZ(2,N).EQ.6 ) THEN
          IF( ABS(SN(4,N)-SN(2,N)).LT.1.D+4*EPSL ) THEN
            NC = NC + 1
            DNR(1,N) = DNR(1,N)*2.D+0
            DNR(2,N) = DNR(2,N)*2.D+0
            PL(3,N) = PL(2,N) + DNR(1,N)
            PN(4,N) = PN(2,N) + DNR(2,N)
            IF( NC.GT.32 ) THEN
              CHMSG = 'No Change in NAPL Saturation: Node: '
              INDX = 15
              IMSG = N
              RLMSG = SN(2,N)
              CALL WRMSGS( INDX )
            END IF
            GOTO 10
          ENDIF
        ENDIF
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDBC4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Read input file for boundary condition information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     Last Modified by MD White, PNNL, November 2, 1999.
!     Last Modified by MD White, PNNL, November 8, 1999.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
#ifdef ecke
      CHARACTER*64 SDUM
#endif
      CHARACTER*64 ADUM,BDUM(LUK+LSOLU),FDUM
      CHARACTER*64 UNTS
      CHARACTER*32 CHTYP(17)
      CHARACTER*512 CHDUM
      REAL*8 VAR(LBTM,LBCV)
      INTEGER ITYP(LUK+LSOLU+1),NCHB(LUK+LSOLU)
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
     &            'Inflow NAPL ','Seepage Face '/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDBC4'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Boundary Conditions Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      NBC = 0
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Boundary Condition Cards'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO 400 NB = 1, NLIN
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
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
        VARB = 'Boundary Type'
        CALL RDCHR(ISTART,ICOMMA,NCHB(IEQW),CHDUM,BDUM(IEQW))
        CALL RDCHR(ISTART,ICOMMA,NCHB(IEQO),CHDUM,BDUM(IEQO))
        IF( IEQC.GT.0 ) THEN
          DO 5 NSL = 1,NSOLU
            BDUM(NSL+LUK) = 'zero flux'
            IDFLT = 1
            VARB = 'Solute Boundary Type'
            CALL RDCHR(ISTART,ICOMMA,NCHB(NSL+LUK),CHDUM,BDUM(NSL+LUK))
    5     CONTINUE
        ENDIF
#ifdef ecke
!
!---      Aqueous species boundary condition types, 
!         allowing for returns in input lines  ---
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
!---      NAPL species boundary condition types, 
!         allowing for returns in input lines  ---
!
          CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
          IF( INDX.EQ.0 ) THEN
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
            ISTART = 1
          ENDIF
          BDUM(NSOLU+LUK+2) = 'zero flux'
          IDFLT = 1
          VARB = 'NAPL Reactive Species Boundary Condition Type: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(NSOLU+LUK+2))
!
!---      Number of reactive species  ---
!
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
          VARB = 'Number of Reactive Species'
          CALL RDINT(ISTART,ICOMMA,CHDUM,IBCSPX(1))
          DO 16 NSPX = 2,IBCSPX(1)+1
            IBCSPX(NSPX) = 0
   16     CONTINUE
!
!---      Loop over number of reactive species  ---
!
          DO 20 NSPX = 1,IBCSPX(1)
!
!---        Allow for returns in input lines  ---
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
!---        Aqueous species  ---
!
            DO 17 M = 1,NSPL
              IF( SPNML(M).EQ.SDUM ) THEN
                IBCSPX(NSPX+1) = M
                GOTO 19
              ENDIF
   17       CONTINUE
!
!---        NAPL species  ---
!
            DO 18 M = 1,NSPN
              IF( SPNMN(M).EQ.SDUM ) THEN
                MX = M + NSPL + NSPS + NSPG
                GOTO 19
              ENDIF
   18       CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Aqueous or NAPL Species Name: '
     &         // SDUM(1:NCH)
            CALL WRMSGS( INDX )
   19       CONTINUE
   20     CONTINUE
!        ENDIF
#endif
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
        ELSE
          INDX = 4
          NCH = NCHB(IEQW)
          CHMSG = 'Unrecognized Boundary Condition: '//BDUM(IEQW)(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Declare dissolved-oil boundary condition type ---
!
        IEQDO = LUK+NSOLU+1
        ITYP(IEQDO) = 1
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
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'inflow napl').NE.0 ) THEN
               ITYP(NSL+LUK) = 16
            ELSEIF(INDEX(BDUM(NSL+LUK)(1:),
     &        'volumetric conc').NE.0 ) THEN
               ITYP(NSL+LUK) = 1
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'aqueous conc').NE.0 )THEN
               ITYP(NSL+LUK) = 8
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
#ifdef ecke
!
!---    Read aqueous reactive species boundary condition type(s) ---
!
        IF( ISLC(40).EQ.1 ) THEN
          IF( INDEX(BDUM(NSOLU+LUK+1)(1:),'outflow').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 7
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'initial co').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 12
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'inflow aqu').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 14
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'inflow napl').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 16
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'inflow').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 13
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),
     &      'volumetric conc').NE.0 ) THEN
             ITYP(NSOLU+LUK+1) = 1
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'aqueous conc').NE.0 )THEN
            ITYP(NSOLU+LUK+1) = 8
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'napl conc').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 10
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'zero flux').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 3
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Reactive Species Boundary Condition: '
     &        //BDUM(NSOLU+LUK+1)
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Read NAPL reactive species boundary condition type(s) ---
!
          IF( INDEX(BDUM(NSOLU+LUK+2)(1:),'outflow').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 7
            IF( ITYP(NSOLU+LUK+1).NE.7 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-NAPL Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'initial co').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 12
            IF( ITYP(NSOLU+LUK+1).NE.1 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.8 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.10 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.12 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-NAPL Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'inflow').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 13
            IF( ITYP(NSOLU+LUK+1).NE.13 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.14 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.16 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-NAPL Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'inflow aqu').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 14
            IF( ITYP(NSOLU+LUK+1).NE.13 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.14 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.16 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-NAPL Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'inflow napl').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 16
            IF( ITYP(NSOLU+LUK+1).NE.13 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.14 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.16 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-NAPL Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),
     &      'volumetric conc').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 1
            IF( ITYP(NSOLU+LUK+1).NE.1 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.8 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.10 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.12 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-NAPL Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'aqueous con').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 8
            IF( ITYP(NSOLU+LUK+1).NE.1 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.8 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.10 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.12 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-NAPL Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'napl conc').NE.0 )THEN
            ITYP(NSOLU+LUK+2) = 10
            IF( ITYP(NSOLU+LUK+1).NE.8 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.10 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.12 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-NAPL Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'zero flux').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 3
            IF( ITYP(NSOLU+LUK+1).NE.3 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-NAPL Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Reactive Species Boundary Condition: '
     &        //BDUM(NSOLU+LUK+2)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
#endif
!
!---  Write boundary condition type(s) ---
!
        WRITE(IWR,'(A)') 'Boundary Condition Type: '
        WRITE(IWR,'(2X,2A)') 'Aqueous: ',CHTYP(ITYP(IEQW))
        WRITE(IWR,'(2X,2A)') 'NAPL:    ',CHTYP(ITYP(IEQO))
        IF( IEQC.GT.0 ) THEN
          DO 30 NSL = 1,NSOLU
            IDB = INDEX( SOLUT(NSL)(1:),'  ') - 1
            WRITE(IWR,'(2X,2A)') CHTYP(ITYP(NSL+LUK)),SOLUT(NSL)(1:IDB)
   30     CONTINUE
        ENDIF
#ifdef ecke
!
!---    Write aqueous and NAPL species boundary condition type(s) ---
!
        IF( ISLC(40).EQ.1 ) THEN
          WRITE(IWR,'(2X,3A)') 'Aqueous: ',
     &      CHTYP(ITYP(NSOLU+LUK+1)),' Reactive Species'
          WRITE(IWR,'(2X,3A)') 'NAPL: ',
     &      CHTYP(ITYP(NSOLU+LUK+2)),' Reactive Species'
        ENDIF
#endif
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
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          VARB = 'Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH)
     &,': ',VAR(NTM,1)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,VAR(NTM,1),INDX)
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
          ELSE
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          ENDIF
          VARB = 'Dissolved Oil Relative Saturation'
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),': '
          ISX = ISTART
          ICX = ICOMMA
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            ISTART = ISX
            ICOMMA = ICX
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            WRITE(IWR,'(1PE11.4)') VAR(NTM,4)
          IF( ITYP(IEQO).EQ.1 ) THEN
            VARB = 'NAPL Pressure'
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
          ELSEIF( ITYP(IEQO).EQ.2 ) THEN
            VARB = 'Volumetric NAPL Flux'
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
          ELSEIF( ITYP(IEQO).EQ.3 ) THEN
            VARB = 'NAPL Pressure'
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
          ELSEIF( ITYP(IEQO).EQ.11 ) THEN
            VARB = 'Base NAPL Pressure'
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
          IF( IEQC.GT.0 ) THEN
            DO 90 NSL = 1,NSOLU
              IDB = INDEX( SOLUT(NSL)(1:),'  ') - 1
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
   90       CONTINUE
          ENDIF
#ifdef ecke
!
!---      Loop over reactive species inputs  ---
!
          DO 92 NSPX = 1,IBCSPX(1)
            NSP = IBCSPX(NSPX+1)
            M = NSOLU+LBCU+NSPX
!
!---        Initial input line  ---
!
            IF( NSPX.EQ.1 ) THEN
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              ISTART = 1
            ENDIF
!
!---        Allow for returns in input lines  ---
!
            CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
            IF( INDX.EQ.0 ) THEN
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              ISTART = 1
            ENDIF
!
!---        Aqueous species  ---
!
            IF( NSP.LE.NSPL ) THEN
              IF( ITYP(NSOLU+LUK+1).EQ.8
     &          .OR. ITYP(NSOLU+LUK+1).EQ.14 ) THEN
                VARB = 'Aqueous Concentration, '
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,M))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                IDB = INDEX( SPNML(NSP)(1:),'  ') - 1
                WRITE(IWR,'(2X,A,2X,3A,1PE11.4,$)') 
     &            SPNML(NSP)(1:IDB),VARB(1:IVR),UNTS(1:NCH),
     &            ': ',VAR(NTM,M)
                INDX = 0
                IUNM = -3
                IUNMOL = 1
                CALL RDUNIT(UNTS,VAR(NTM,M),INDX)
!
!---            Convert aqueous concentration from kmol/m^3 to
!               mol/m^3  ---
!
                VAR(NTM,M) = VAR(NTM,M)*1.D+3
                WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &            VAR(NTM,M),', mol/m^3)'
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
!---        NAPL species  ---
!
            ELSEIF( NSP.GT.(NSPL+NSPS+NSPG) .AND. 
     &        NSP.LE.(NSPL+NSPS+NSPG+NSPN) ) THEN
              IF( ITYP(NSOLU+LUK+2).EQ.8
     &          .OR. ITYP(NSOLU+LUK+2).EQ.10
     &          .OR. ITYP(NSOLU+LUK+2).EQ.16 ) THEN
                VARB = 'NAPL Concentration, '
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,M))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                IDB = INDEX( SPNML(NSP)(1:),'  ') - 1
                WRITE(IWR,'(2X,A,2X,3A,1PE11.4,$)') 
     &            SPNML(NSP)(1:IDB),VARB(1:IVR),UNTS(1:NCH),
     &            ': ',VAR(NTM,M)
                INDX = 0
                IUNM = -3
                IUNMOL = 1
                CALL RDUNIT(UNTS,VAR(NTM,M),INDX)
!
!---            Convert NAPL concentration from kmol/m^3 to
!               mol/m^3  ---
!
                VAR(NTM,M) = VAR(NTM,M)*1.D+3
                WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &            VAR(NTM,M),', mol/m^3)'
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
   92     CONTINUE
#endif
!
!---  Check for nonascending boundary condition times  ---
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
        DO 108 NTM = 1,IBCMX
          DO 102 M = 1,LBCU
            BC(M,NTM,NB) = VAR(NTM,M)
  102     CONTINUE
          DO 104 NSL = 1,NSOLU
            BC(NSL+LBCU,NTM,NB) = VAR(NTM,NSL+LBCU)
  104     CONTINUE
#ifdef ecke
          IF( ISLC(40).EQ.1 ) THEN
            DO 106 NSPX = 1,IBCSPX(1)
              M = NSOLU+LBCU+NSPX
              BC(M,NTM,NB) = VAR(NTM,M)
  106       CONTINUE
          ENDIF
#endif
  108   CONTINUE
!
!---  Assign values to boundary variables  ---
!
        NBCL = 0
        DO 320 K = K1X, K2X
          DO 310 J = J1X, J2X
            DO 300 I = I1X, I2X
              N = ND(I,J,K)
              IX = I
              JX = J
              KX = K
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
              IBCT(IEQW,NBC) = ITYP(IEQW)
              IBCT(IEQO,NBC) = ITYP(IEQO)
              IBCT(IEQDO,NBC) = ITYP(IEQDO)
              IF( IEQC.GT.0 ) THEN
                DO 110 NSL = 1,NSOLU
                  IBCT(NSL+LUK,NBC) = ITYP(NSL+LUK)
  110           CONTINUE
              ENDIF
#ifdef ecke
              IF( ISLC(40).EQ.1 ) THEN
                IBCT(NSOLU+LUK+1,NBC) = ITYP(NSOLU+LUK+1)
                IBCT(NSOLU+LUK+2,NBC) = ITYP(NSOLU+LUK+2)
                DO 120 NSP = 1,LSPBC+1
                  IBCSP(NSP,NBC) = IBCSPX(NSP)
  120           CONTINUE
              ENDIF
#endif
              IBCM(NBC) = IBCMX
              IBCIN(NBC) = NB
!
!---  Assign hydraulic gradient type boundary condition
!     for the base surface  ---
!
              IF( ITYP(IEQW).EQ.11 ) THEN
                IF( NBCL.EQ.1 ) THEN
                  IBCT(IEQW,NBC) = -11
                  NBHG(1,NBC) = NBC
                ELSE
                  NBHG(1,NBC) = NBHG(1,NBC-1)
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDBC4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDOT4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Reads the oil transport card for solid-aqueous 
!     partition coefficients and hydraulic dispersivities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1996.
!     Last Modified by MD White, PNNL, 1 August 1996.
!     Last Modified by MD White, PNNL, 8 January 2004.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/RDOT4'
!
!---  Write card information to output file  ---
!
      CARD = 'Oil Transport Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      ISLC(12) = 1
      IDSPD = 0
!
!---  Check for advection-dispersion-diffusion
!     transport scheme, default to upwind  ---
!
      IADTSX = 0
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
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
          IADTSX = 1
        ENDIF
        ISTART = 1
        VARB = 'Rock Name: '
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
        WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': '
     &,DPLD(IROCK)
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
        WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': '
     &,DPTD(IROCK)
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
        ENDIF
!
!---  Read next rock/soil type or scaling group  ---
!
        IF( N.LT.NROCK ) WRITE(IWR,'(/)')
        GOTO 10
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDOT4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIC4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Read input file for initial conditions information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!     Last Modified by MD White, PNNL, November 2, 1997.
!     Last Modified by MD White, PNNL, November 2, 1999.
!     Last Modified by SK Wurstner, PNNL, December 13, 2007.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SPILL
      USE SOLTN
      USE REACT
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
      CHARACTER*64 ADUM,BDUM,FDUM,UNTS
      CHARACTER*24 CHLB(3)
      CHARACTER*512 CHDUM
      INTEGER IDOM(6)
      REAL*8 VAR(5)
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
      SUB_LOG(ISUB_LOG) = '/RDIC4'
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
        ELSEIF( INDEX( ADUM(1:),'residual napl' ).NE.0 ) THEN
          VARB = 'Initial Residual NAPL Saturation'
        ELSEIF( INDEX( ADUM(1:),'trapped gas' ).NE.0 ) THEN
          VARB = 'Initial Trapped Gas Saturation'
        ELSEIF( INDEX( ADUM(1:),'oil conc' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-Oil Concentration'
        ELSEIF( INDEX( ADUM(1:),'oil mole' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-Oil Mole Fraction'
        ELSEIF( INDEX( ADUM(1:),'oil mass' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-Oil Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'napl' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'spill' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'height' ).NE.0 ) THEN
          VARB = 'Initial NAPL Spill Height'
          IUNM = 1
        ELSEIF( INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'spill' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'height' ).NE.0 ) THEN
          VARB = 'Initial Aqueous Spill Height'
          IUNM = 1
        ELSEIF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          VARB = 'Initial Solute Concentration'
          IUNM = -3
#ifdef ecke
        ELSEIF( INDEX( ADUM(1:),'specie' ).NE.0 ) THEN
          VARB = 'Reactive Species Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          VARB = 'Initial Reactive Species Concentration'
!
!---      Set species units  ---
!
          IUNMOL = 1
          IF( INDEX(ADUM(1:),'napl').NE.0 ) THEN
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
#endif
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
          WRITE(IWR,'(2X,4A,1PE11.4)') ADUM(1:NCHA),
     &      ' (Default Value), ',UNTS(1:NCHU),': ',VAR(1)
          CALL RDCHR(ISTART,ICOMMA,NCHF,CHDUM,FDUM)
          NCHF = INDEX(FDUM,'  ')-1
          OPEN(UNIT=26,FILE=FDUM(1:NCHF),STATUS='OLD',FORM='FORMATTED')
          WRITE(IWR,'(/,2A)') 'Initial Conditions File: ',FDUM(1:NCHF)
          INDX = 0
          CALL RDUNIT( UNTS,VAR(1),INDX )
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
!---  Read aqueous pressure  ---
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
!---    Read NAPL pressure or auto prime NAPL pressures  ---
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
!---    Read aqueous-trapped NAPL saturation  ---
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
!---    Read residual NAPL saturation  ---
!
        ELSEIF( INDEX(ADUM(1:),'residual napl').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            CALL RDINFS( SNR,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( SNR,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( SNR,VAR,ADDER,IDOM,INDX )
          ENDIF
!
!---    Read temperature  ---
!
        ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = 1
            CALL RDINFS( T,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( T,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( T,VAR,ADDER,IDOM,INDX )
          ENDIF
!
!---    Read NAPL spill height  ---
!
        ELSEIF( INDEX( ADUM(1:),'napl' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'spill' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'height' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Option with NAPL Spill Height: ' //
     &        ADUM(1:NCHA)
            CALL WRMSGS( INDX )
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            INDX = 4
            CHMSG = 'Unrecognized Option with NAPL Spill Height: '//
     &        ADUM(1:NCHA)
            CALL WRMSGS( INDX )
          ELSE
            CALL RDINIS( SN,VAR,ADDER,IDOM,INDX )
          ENDIF
          IF( IDOM(5).NE.KFLD .AND. IDOM(6).NE.KFLD ) THEN
            INDX = 4
            CHMSG = 'Specify Domain Top with NAPL Spill Height: '//
     &        ADUM(1:NCHA)
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Load spill height variable  ---
!
          K = KFLD
          DO 202 J=IDOM(3),IDOM(4)
          DO 202 I=IDOM(1),IDOM(2)
            NS = (J-1)*IFLD + I
            N = ND(I,J,K)
            HNSP(2,NS) = SN(2,N)
 202      CONTINUE
!
!---    Read Aqueous spill height  ---
!
        ELSEIF( INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'spill' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'height' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Option with Aqueous ' //
     &        'Spill Height: ' // ADUM(1:NCHA)
            CALL WRMSGS( INDX )
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            INDX = 4
            CHMSG = 'Unrecognized Option with Aqueous ' //
     &        'Spill Height: ' // ADUM(1:NCHA)
            CALL WRMSGS( INDX )
          ELSE
            CALL RDINIS( SN,VAR,ADDER,IDOM,INDX )
          ENDIF
          IF( IDOM(5).NE.KFLD .AND. IDOM(6).NE.KFLD ) THEN
            INDX = 4
            CHMSG = 'Specify Domain Top with Aqueous' // 
     &        ' Spill Height: '// ADUM(1:NCHA)
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Load spill height variable  ---
!
          K = KFLD
          DO 204 J=IDOM(3),IDOM(4)
          DO 204 I=IDOM(1),IDOM(2)
            NS = (J-1)*IFLD + I
            N = ND(I,J,K)
            HLSP(2,NS) = SN(2,N)
 204      CONTINUE
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
            IUNK = 1
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
            IUNK = 1
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
!---    Read solute concentrations  ---
!
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'napl').NE.0 ) THEN
            IVAR = 3
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
            IVAR = 2
          ELSE
            IVAR = 1
          ENDIF
          IF( INDEX( UNTS(1:),'bd' ).NE.0 ) IVAR = -IVAR
          DO 240 NSL = 1,NSOLU
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
              GOTO 250
            ENDIF
  240     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute: '//BDUM(1:NCHB)
          CALL WRMSGS( INDX )
  250     CONTINUE
#ifdef ecke
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
!---      NAPL reactive species  ---
!
          DO 530 NSPX = 1,NSPN
            IDB = INDEX(SPNMN(NSPX)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SPNMN(NSPX)(1:IDB) ) THEN
              NSP = NSPX + NSPL + NSPS + NSPG
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
#endif
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: '//
     &      ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIC4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINPT4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Read input file cards.
!     Direct control to card reader subroutines.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/RDINTP4'
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
      INDX = 18
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
      INDX = 18
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
      INDX = 18
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
        CALL RDSP4
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
      INDX = 18
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
        GOTO 1300
      ELSE
        GOTO 1000
      ENDIF
 1010 CONTINUE
      INDX = 18
      CHMSG = 'Missing NAPL Relative Permeability Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for oil properties card --
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
      INDX = 18
      CHMSG = 'Missing Oil Properties Card'
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
        CALL RDTF4
        REWIND(IRD)
        GOTO 2100
      ELSE
        GOTO 2000
      ENDIF
 2010 CONTINUE
      IF( IEQC.EQ.0 ) THEN
        REWIND(IRD)
      ELSE
        INDX = 18
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
        CALL RDTP4
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
        CHMSG = 'Missing Solute/Rock Interaction Card'
        CALL WRMSGS( INDX )
      ENDIF
 2200 CONTINUE
#ifdef ecke
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
        GOTO 2259
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
!---  Search input file for NAPL species card  ---
!
 2259 READ(IRD,'(A)', END=2260) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2259
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'napl').NE.0 .AND.
     &  INDEX(CHDUM(2:),'specie').NE.0) THEN
        CALL RDNPSP
        REWIND(IRD)
        GOTO 2300
      ELSE
        GOTO 2259
      ENDIF
 2260 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing NAPL Species Card'
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
!---  Search input file for equilibrium reactions card  ---
!
 2409 READ(IRD,'(A)', END=2410) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2409
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'equil').NE.0 .AND.
     &  INDEX(CHDUM(2:),'react').NE.0) THEN
        CALL RDEQRC
        REWIND(IRD)
        GOTO 2500
      ELSE
        GOTO 2400
      ENDIF
 2410 CONTINUE
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
!---  Search input file for kinetic equations card  ---
!
 2800 CONTINUE
 2809 READ(IRD,'(A)', END=2810) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2809
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'kinetic').NE.0 .AND.
     &  INDEX(CHDUM(2:),'equat').NE.0) THEN
        CALL RDKNEQ
        REWIND(IRD)
        GOTO 2900
      ELSE
        GOTO 2800
      ENDIF
 2810 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Kinetic Equations Card'
        CALL WRMSGS( INDX )
      ENDIF
      REWIND(IRD)
!
!---  Search input file for conservation equations card  ---
!
 2900 CONTINUE
 2909 READ(IRD,'(A)', END=2910) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2909
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'conservation').NE.0 .AND.
     &  INDEX(CHDUM(2:),'equat').NE.0) THEN
        CALL RDCNEQ
        REWIND(IRD)
        GOTO 3000
      ELSE
        GOTO 2900
      ENDIF
 2910 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Conservation Equations Card'
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
#endif
!
!---  Search input file for oil transport card --
!
 3200 CONTINUE
 3209 READ(IRD,'(A)', END=3210) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3209
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    ( INDEX(CHDUM(2:),'oil transport').NE.0 .OR.
     &    INDEX(CHDUM(2:),'voc transport').NE.0 )) THEN
        CALL RDOT4
        REWIND(IRD)
        GOTO 4000
      ELSE
        GOTO 3200
      ENDIF
 3210 CONTINUE
      INDX = 1
      CHMSG = 'Missing Oil Transport Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for initial conditions card --
!
 4000 CONTINUE
 4009 READ(IRD,'(A)', END=4010) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4009
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'initial').NE.0 ) THEN
        CALL RDIC4
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
        INDX = 18
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
        CALL RDBC4
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
        CALL RDSR4
        REWIND(IRD)
        IF( LWELL.EQ.1 ) THEN
          GOTO 4250
        ELSE
          GOTO 4300
        ENDIF
      ELSE
        GOTO 4200
      ENDIF
 4210 CONTINUE
      INDX = 1
      CHMSG = 'Missing Source Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for well card ---
!
 4250 CONTINUE
 4259 READ(IRD,'(A)', END= 4260) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4259
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'well').NE.0 ) THEN
        CALL RDWELL4
        REWIND(IRD)
        GOTO 4300
      ELSE
        GOTO 4250
      ENDIF
 4260 CONTINUE
      INDX = 1
      CHMSG = 'Missing Well Card'
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
        CALL RDOU4
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
        CALL RDSF4
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
!---  End of RDINPT4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDOU4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Read input file for output information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE OUTPU
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
      CHARACTER*64 ADUM,UNTS,SOLNM
      CHARACTER*512 CHDUM
      CHARACTER*6 FORM
      LOGICAL ISKIP
      EXTERNAL ICOUNT
!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM
      DATA FORM / '(I6,$)' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDOU4'
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
      WRITE(IWR,'(/,A,I6)') 'Reference Node No. and Indices: ',NREF
      DO 100 N = 1,NREF
        ISTART = 1
        CALL RDINPL( CHDUM )
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
          CHMSG = 'Unrecognized Solute: '//SOLNM(1:NCS)
          CALL WRMSGS( INDX )
  120     CONTINUE
        ENDIF
#ifdef ecke
        IF( INDEX( ADUM(1:),'species' ).NE.0 ) THEN
          IF( ISLC(40).EQ.0 ) THEN
            NVC = NVC -1
            GOTO 200
          ENDIF
          VARB = 'Reference Node Variable: Reactive Species Name: '
          CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SPNM)
!
!---      Conservation- or kinetic-component species  ---
!
          IF( INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
            DO 130 NSL = NSOLU+1,NSOLU+NEQC+NEQK
              IF( SPNM.EQ.SOLUT(NSL) ) GOTO 150
  130       CONTINUE
          ENDIF
!
!---      Aqueous species  ---
!
          DO 132 M = 1,NSPL
            NSP = M
            IF( SPNM.EQ.SPNML(M) ) GOTO 150
  132     CONTINUE
!
!---      Solid species  ---
!
          DO 134 M = 1,NSPS
            NSP = M+NSPL
            IF( SPNM.EQ.SPNMS(M) ) GOTO 150
  134     CONTINUE
!
!---      NAPL species  ---
!
          DO 136 M = 1,NSPN
            NSP = M+NSPL+NSPS+NSPG
            IF( SPNM.EQ.SPNMN(M) ) GOTO 150
  136     CONTINUE
!
!---      Unrecognized species name  ---
!
          INDX = 4
          CHMSG = 'Unrecognized Reference-Node Reactive Species Name: '
     &      // SPNM
          CALL WRMSGS( INDX )
          NVC = NVC -1
          GOTO 200
  150     CONTINUE
        ENDIF
#endif
        ISKIP = .FALSE.
        IF( INDEX(ADUM(1:),'coupled-well').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'total-liquid depth').NE.0 ) THEN
            IREF(NV) = 177
          ELSEIF( INDEX(ADUM(1:),'aqueous depth').NE.0 ) THEN
            IREF(NV) = 178
          ELSEIF( INDEX(ADUM(1:),'napl depth').NE.0 ) THEN
            IREF(NV) = 179
          ELSEIF( INDEX(ADUM(1:),'pressure').NE.0 ) THEN
            IREF(NV) = 180
          ELSEIF( INDEX(ADUM(1:),'aqueous sat').NE.0 ) THEN
            IREF(NV) = 181
          ELSEIF( INDEX(ADUM(1:),'napl sat').NE.0 ) THEN
            IREF(NV) = 182
            CHREF(182) = 'SNW'
          ELSEIF( INDEX(ADUM(1:),'oil water mass').NE.0 ) THEN
            IREF(NV) = 183
            CHREF(183) = 'XLOW'
          ELSEIF( INDEX(ADUM(1:),'total-liquid pumping rat').NE.0 ) THEN
            IREF(NV) = 184
            CHREF(184) = 'QT_W'
            IUNM = 3
            IUNS = -1
            ISKIP = .TRUE.
          ELSEIF( INDEX(ADUM(1:),'aqueous pumping rat').NE.0 ) THEN
            IREF(NV) = 185
            CHREF(185) = 'QL_W'
            IUNM = 3
            IUNS = -1
            ISKIP = .TRUE.
          ELSEIF( INDEX(ADUM(1:),'napl pumping rat').NE.0 ) THEN
            IREF(NV) = 186
            CHREF(186) = 'QN_W'
            IUNM = 3
            IUNS = -1
            ISKIP = .TRUE.
          ELSEIF( INDEX(ADUM(1:),'total-liquid pumping int').NE.0 ) THEN
            IREF(NV) = 187
            CHREF(187) = 'QTWI'
            IUNM = 3
            ISKIP = .TRUE.
          ELSEIF( INDEX(ADUM(1:),'aqueous pumping int').NE.0 ) THEN
            IREF(NV) = 188
            CHREF(188) = 'QLWI'
            IUNM = 3
            ISKIP = .TRUE.
          ELSEIF( INDEX(ADUM(1:),'napl pumping int').NE.0 ) THEN
            IREF(NV) = 189
            CHREF(189) = 'QNWI'
            IUNM = 3
            ISKIP = .TRUE.
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Reference Node Variable: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous pressure').NE.0 ) THEN
          IREF(NV) = 1
        ELSEIF( INDEX(ADUM(1:),'napl pressure').NE.0 ) THEN
          IREF(NV) = 3
        ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 ) THEN
          IREF(NV) = 4
        ELSEIF( INDEX(ADUM(1:),'phase condition').NE.0 ) THEN
          IREF(NV) = 5
        ELSEIF( INDEX(ADUM(1:),'aqueous gauge pressure').NE.0 ) THEN
          IREF(NV) = 6
        ELSEIF( INDEX(ADUM(1:),'napl gauge pressure').NE.0 ) THEN
          IREF(NV) = 8
        ELSEIF( INDEX(ADUM(1:),'apparent aqueous sat').NE.0 ) THEN
          IREF(NV) = 9
        ELSEIF( INDEX(ADUM(1:),'apparent total sat').NE.0 ) THEN
          IREF(NV) = 10
        ELSEIF( INDEX(ADUM(1:),'residual napl sat').NE.0 ) THEN
          IREF(NV) = 210
        ELSEIF( INDEX(ADUM(1:),'mobile napl sat').NE.0 ) THEN
          IREF(NV) = 211
        ELSEIF( INDEX(ADUM(1:),'free napl sat').NE.0 ) THEN
          IREF(NV) = 212
        ELSEIF( INDEX(ADUM(1:),'effective trapped napl').NE.0 ) THEN
          IREF(NV) = 18
        ELSEIF( INDEX(ADUM(1:),'trapped napl sat').NE.0 ) THEN
          IREF(NV) = 106
        ELSEIF( INDEX(ADUM(1:),'aqueous sat').NE.0 ) THEN
          IREF(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'napl sat').NE.0 ) THEN
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
        ELSEIF( INDEX(ADUM(1:),'napl hydraulic head').NE.0 ) THEN
          IREF(NV) = 29
        ELSEIF( INDEX(ADUM(1:),'rock/soil type').NE.0 ) THEN
          IREF(NV) = 30
        ELSEIF( INDEX(ADUM(1:),'aqueous relative perm').NE.0 ) THEN
          IREF(NV) = 31
        ELSEIF( INDEX(ADUM(1:),'napl relative perm').NE.0 ) THEN
          IREF(NV) = 33
        ELSEIF( INDEX(ADUM(1:),'aqueous density').NE.0 ) THEN
          IREF(NV) = 34
        ELSEIF( INDEX(ADUM(1:),'napl density').NE.0 ) THEN
          IREF(NV) = 36
        ELSEIF( INDEX(ADUM(1:),'total water mass').NE.0 ) THEN
          IREF(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'total oil mass').NE.0 ) THEN
          IREF(NV) = 39
        ELSEIF( INDEX(ADUM(1:),'water mass source int').NE.0 ) THEN
          IREF(NV) = 40
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
        ELSEIF( INDEX(ADUM(1:),'xnc aqueous vol').NE.0 ) THEN
          IREF(NV) = 87
        ELSEIF( INDEX(ADUM(1:),'ync aqueous vol').NE.0 ) THEN
          IREF(NV) = 88
        ELSEIF( INDEX(ADUM(1:),'znc aqueous vol').NE.0 ) THEN
          IREF(NV) = 89
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
        ELSEIF( INDEX(ADUM(1:),'oil mass source rate').NE.0 ) THEN
          IREF(NV) = 142
        ELSEIF( INDEX(ADUM(1:),'well depth').NE.0 ) THEN
          IREF(NV) = 144
        ELSEIF( INDEX(ADUM(1:),'aqueous well flow rate').NE.0 ) THEN
          IREF(NV) = 145
        ELSEIF( INDEX(ADUM(1:),'aqueous well flow integral').NE.0 ) THEN
          IREF(NV) = 146
        ELSEIF( INDEX(ADUM(1:),'napl well flow rate').NE.0 ) THEN
          IREF(NV) = 154
        ELSEIF( INDEX(ADUM(1:),'napl well flow integral').NE.0 ) THEN
          IREF(NV) = 155
        ELSEIF( INDEX(ADUM(1:),'total well flow rate').NE.0 ) THEN
          IREF(NV) = 156
        ELSEIF( INDEX(ADUM(1:),'total well flow integral').NE.0 ) THEN
          IREF(NV) = 157
        ELSEIF( INDEX(ADUM(1:),'integrated water mass').NE.0 ) THEN
          IREF(NV) = 191
        ELSEIF( INDEX(ADUM(1:),'integrated oil mass').NE.0 ) THEN
          IREF(NV) = 193
        ELSEIF( INDEX(ADUM(1:),'integrated aqueous water').NE.0 ) THEN
          IREF(NV) = 194
        ELSEIF( INDEX(ADUM(1:),'integrated aqueous oil').NE.0 ) THEN
          IREF(NV) = 196
        ELSEIF( INDEX(ADUM(1:),'sorbed oil mass frac').NE.0 ) THEN
          IREF(NV) = 241
        ELSEIF( INDEX(ADUM(1:),'sorbed oil mass').NE.0 ) THEN
          IREF(NV) = 240
        ELSEIF( INDEX(ADUM(1:),'sorbed oil volumetric conc').NE.0 ) THEN
          IREF(NV) = 242
        ELSEIF( INDEX(ADUM(1:),'aqueous spill').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'head').NE.0 .OR.
     &    INDEX(ADUM(1:),'height').NE.0 ) ) THEN
          IREF(NV) = 287
        ELSEIF( INDEX(ADUM(1:),'napl spill').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'head').NE.0 .OR.
     &    INDEX(ADUM(1:),'height').NE.0 ) ) THEN
          IREF(NV) = 288
#ifdef ecke
        ELSEIF( (INDEX(ADUM(1:),'solute volumetric conc').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species volumetric conc').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 1
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 1
            CHREF(INDX) = 'SP'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute volumetric conc').NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 1
#endif
#ifdef ecke
        ELSEIF( (INDEX(ADUM(1:),'solute aqueous conc').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species aqueous conc').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 2
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 2
            CHREF(INDX) = 'SPL'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute aqueous conc').NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 2
#endif
#ifdef ecke
        ELSEIF( (INDEX(ADUM(1:),'solute napl conc').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species napl conc').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 4
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 4
            CHREF(INDX) = 'SPN'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute napl conc').NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 4
#endif
        ELSEIF( INDEX(ADUM(1:),'solute aqueous mol').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 5
        ELSEIF( INDEX(ADUM(1:),'solute napl mol').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 7
        ELSEIF( INDEX(ADUM(1:),'x solute flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 8
        ELSEIF( INDEX(ADUM(1:),'y solute flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 9
        ELSEIF( INDEX(ADUM(1:),'z solute flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 10
#ifdef ecke
        ELSEIF( (INDEX(ADUM(1:),'solute source').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species source').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 11
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 11
            CHREF(INDX) = 'SPSR'
            UNREF(INDX) = 'mol/s'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute source').NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 11
#endif
#ifdef ecke
        ELSEIF( (INDEX(ADUM(1:),'solute integrated mass').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species integrated mass').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 23
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 23
            CHREF(INDX) = 'SPIM'
            UNREF(INDX) = 'mol'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute integrated mass').NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 23
#endif
#ifdef ecke
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
        ELSEIF( INDEX(ADUM(1:),'species napl conc').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 4
          IREF(NV) = INDX
          CHREF(INDX) = 'SPN'
          UNPLOT(INDX) = 'mol/m^3'
          UNREF(INDX) = 'mol/m^3'
        ELSEIF( INDEX(ADUM(1:),'species source').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 11
          IREF(NV) = INDX
        ELSEIF( INDEX(ADUM(1:),'species integrated mass').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 23
          IREF(NV) = INDX
          CHREF(INDX) = 'SPIM'
          UNREF(INDX) = 'mol'
#endif
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
#ifdef ecke
        ELSEIF( INDEX( ADUM(1:),'species' ).NE.0 ) THEN
          WRITE( IWR,'(2X,3A,2X,2A,I2,A)' ) ADUM(1:NCH),', ',
     &      UNREF(IREF(NV))(1:NCU),SPNM(1:NCS),' Species(',NSP,')'
#endif
        ELSE
          WRITE( IWR,'(2X,3A)' ) ADUM(1:NCH),', ',UNREF(IREF(NV))(1:NCU)
        ENDIF
        IF( .NOT. ISKIP ) CALL RDOUUN( IREF(NV) )
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
          CHMSG = 'Unrecognized Solute: '//SOLNM(1:NCS)
          CALL WRMSGS( INDX )
  320     CONTINUE
        ENDIF
#ifdef ecke
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
!---      NAPL species  ---
!
          DO 336 M = 1,NSPN
            NSP = M+NSPL+NSPS+NSPG
            IF( SPNM.EQ.SPNMN(M) ) GOTO 350
  336     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Plot File Reactive Species Name: '
     &      // SPNM
          CALL WRMSGS( INDX )
          NVC = NVC -1
          GOTO 400
  350     CONTINUE
        ENDIF
#endif
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
        ELSEIF( INDEX(ADUM(1:),'napl pressure').NE.0 ) THEN
          IPLOT(NV) = 3
        ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 ) THEN
          IPLOT(NV) = 4
        ELSEIF( INDEX(ADUM(1:),'phase condition').NE.0 ) THEN
          IPLOT(NV) = 5
        ELSEIF( INDEX(ADUM(1:),'aqueous gauge pressure').NE.0 ) THEN
          IPLOT(NV) = 6
        ELSEIF( INDEX(ADUM(1:),'napl gauge pressure').NE.0 ) THEN
          IPLOT(NV) = 8
        ELSEIF( INDEX(ADUM(1:),'apparent aqueous sat').NE.0 ) THEN
          IPLOT(NV) = 9
        ELSEIF( INDEX(ADUM(1:),'apparent total sat').NE.0 ) THEN
          IPLOT(NV) = 10
        ELSEIF( INDEX(ADUM(1:),'effective trapped napl').NE.0 ) THEN
          IPLOT(NV) = 18
        ELSEIF( INDEX(ADUM(1:),'trapped napl sat').NE.0 ) THEN
          IPLOT(NV) = 106
        ELSEIF( INDEX(ADUM(1:),'residual napl sat').NE.0 ) THEN
          IPLOT(NV) = 210
        ELSEIF( INDEX(ADUM(1:),'mobile napl sat').NE.0 ) THEN
          IPLOT(NV) = 211
        ELSEIF( INDEX(ADUM(1:),'free napl sat').NE.0 ) THEN
          IPLOT(NV) = 212
        ELSEIF( INDEX(ADUM(1:),'aqueous sat').NE.0 ) THEN
          IPLOT(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'napl sat').NE.0 ) THEN
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
        ELSEIF( INDEX(ADUM(1:),'napl hydraulic head').NE.0 ) THEN
          IPLOT(NV) = 29
        ELSEIF( INDEX(ADUM(1:),'rock/soil type').NE.0 ) THEN
          IPLOT(NV) = 30
        ELSEIF( INDEX(ADUM(1:),'aqueous relative perm').NE.0 ) THEN
          IPLOT(NV) = 31
        ELSEIF( INDEX(ADUM(1:),'napl relative perm').NE.0 ) THEN
          IPLOT(NV) = 33
        ELSEIF( INDEX(ADUM(1:),'aqueous density').NE.0 ) THEN
          IPLOT(NV) = 34
        ELSEIF( INDEX(ADUM(1:),'napl density').NE.0 ) THEN
          IPLOT(NV) = 36
        ELSEIF( INDEX(ADUM(1:),'total water mass').NE.0 ) THEN
          IPLOT(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'total oil mass').NE.0 ) THEN
          IPLOT(NV) = 39
        ELSEIF( INDEX(ADUM(1:),'water mass source int').NE.0 ) THEN
          IPLOT(NV) = 40
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
        ELSEIF( INDEX(ADUM(1:),'xnc aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 87
        ELSEIF( INDEX(ADUM(1:),'ync aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 88
        ELSEIF( INDEX(ADUM(1:),'znc aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 89
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
        ELSEIF( INDEX(ADUM(1:),'oil mass source rate').NE.0 ) THEN
          IPLOT(NV) = 142
        ELSEIF( INDEX(ADUM(1:),'sorbed oil mass frac').NE.0 ) THEN
          IPLOT(NV) = 241
        ELSEIF( INDEX(ADUM(1:),'sorbed oil mass').NE.0 ) THEN
          IPLOT(NV) = 240
        ELSEIF( INDEX(ADUM(1:),'sorbed oil volumetric conc').NE.0 ) THEN
          IPLOT(NV) = 242
        ELSEIF( INDEX(ADUM(1:),'aqueous spill').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'head').NE.0 .OR.
     &    INDEX(ADUM(1:),'height').NE.0 ) ) THEN
          IPLOT(NV) = 287
        ELSEIF( INDEX(ADUM(1:),'napl spill').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'head').NE.0 .OR.
     &    INDEX(ADUM(1:),'height').NE.0 ) ) THEN
          IPLOT(NV) = 288
#ifdef ecke
        ELSEIF( (INDEX(ADUM(1:),'solute volumetric conc').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species volumetric conc').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 1
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 1
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute volumetric conc').NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 1
#endif
#ifdef ecke
        ELSEIF( (INDEX(ADUM(1:),'solute aqueous conc').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species aqueous conc').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 2
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 2
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute aqueous conc').NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 2
#endif
#ifdef ecke
        ELSEIF( (INDEX(ADUM(1:),'solute napl conc').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species napl conc').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 4
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 4
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute napl conc').NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 4
#endif
        ELSEIF( INDEX(ADUM(1:),'solute aqueous mol').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 5
        ELSEIF( INDEX(ADUM(1:),'solute napl mol').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 7
        ELSEIF( INDEX(ADUM(1:),'x solute flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 8
        ELSEIF( INDEX(ADUM(1:),'y solute flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 9
        ELSEIF( INDEX(ADUM(1:),'z solute flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 10
#ifdef ecke
        ELSEIF( (INDEX(ADUM(1:),'solute source').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species source').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 11
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 11
            UNPLOT(INDX) = 'mol/s'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute source').NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 11
#endif
#ifdef ecke
        ELSEIF( INDEX(ADUM(1:),'species volumetric conc').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 1
          IPLOT(NV) = INDX
          UNPLOT(INDX) = 'mol/m^3'
        ELSEIF( INDEX(ADUM(1:),'species aqueous conc').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 2
          IPLOT(NV) = INDX
          UNPLOT(INDX) = 'mol/m^3'
        ELSEIF( INDEX(ADUM(1:),'species napl conc').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 4
          IPLOT(NV) = INDX
          UNPLOT(INDX) = 'mol/m^3'
        ELSEIF( INDEX(ADUM(1:),'species source').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 11
          IPLOT(NV) = INDX
          UNPLOT(INDX) = 'mol'
#endif
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDOU4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSF4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Read input file surface flux information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     Last Modified by MD White, Battelle, November 3, 1999.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/RDSF4'
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
            WRITE(IWR,'(A)') 'Aqueous Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            ISFT(NS) = 2
            WRITE(IWR,'(A)') 'Aqueous Volumetric Flux Surface'
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
        ELSEIF( INDEX(ADUM(1:),'dissolved oil').NE.0) THEN
          ISFT(NS) = 9
          WRITE(IWR,'(A)') 'Dissolved Oil Flux Surface'
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
        IF( ISFT(NS).EQ.2 ) THEN
          IUNM = 3
          IUNS = -1
        ELSEIF( ISFT(NS).EQ.4 ) THEN
          IUNM = 3
          IUNS = -1
        ELSEIF( ISFT(NS).EQ.5 ) THEN
          IUNKG = 1
          IUNS = -1
        ELSEIF( ISFT(NS).EQ.7 ) THEN
          IUNKG = 1
          IUNS = -1
        ELSEIF( ISFT(NS).EQ.9 ) THEN
          IUNKG = 1
          IUNS = -1
        ELSEIF( ISFT(NS).EQ.10 ) THEN
          IUNKG = 1
          IUNS = -1
        ELSEIF( ISFT(NS).GT.100 ) THEN
          IUNS = -1
        ENDIF
        VAR = 0.D+0
        INDX = 0
        CALL RDUNIT(UNSF(1,NS),VAR,INDX)
        IDFLT = 1
        VARB = 'Surface Flux Integral Variable Unit'
        CALL RDCHR(ISTART,ICOMMA,NCU,CHDUM,UNSF(2,NS))
        IF( ISFT(NS).EQ.2 ) THEN
          IUNM = 3
        ELSEIF( ISFT(NS).EQ.4 ) THEN
          IUNM = 3
        ELSEIF( ISFT(NS).EQ.5 ) THEN
          IUNKG = 1
        ELSEIF( ISFT(NS).EQ.7 ) THEN
          IUNKG = 1
        ELSEIF( ISFT(NS).EQ.9 ) THEN
          IUNKG = 1
        ELSEIF( ISFT(NS).EQ.10 ) THEN
          IUNKG = 1
        ENDIF
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
        WRITE(IWR,'(2A)') VARB(1:IVR),': '
        WRITE (IWR,'(2X,2(A,I6))') 'I = ',ISFC(1,NS),' to ',ISFC(2,NS)
        WRITE (IWR,'(2X,2(A,I6))') 'J = ',ISFC(3,NS),' to ',ISFC(4,NS)
        WRITE (IWR,'(2X,2(A,I6))') 'K = ',ISFC(5,NS),' to ',ISFC(6,NS)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSF4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSP4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Read input file for rock/soil saturation function information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     Last Modified by MD White, PNNL, 5 October 2001.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/RDSP4'
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
!          CA_NLX = CA_NL
!          CA_GNX = CA_GN
!          IF( M.EQ.2 ) CA_NLX = CA_NL + 1.D-6
!          IF( M.EQ.3 ) CA_GNX = CA_GN + 1.D-6
!          RCOSX(M) = SCHR(5,1) - SCHR(6,1)*COS(CA_NLX) -
!     &      SCHR(7,1)*COS(CA_GNX)
!          RSINX(M) = SCHR(6,1)*SIN(CA_NLX) - SCHR(7,1)*SIN(CA_GNX)
!   30   CONTINUE
!        DRCOSX(1) = 1.D+6*(RCOSX(2)-RCOSX(1))
!        DRCOSX(2) = 1.D+6*(RCOSX(3)-RCOSX(1))
!        DRSINX(1) = 1.D+6*(RSINX(2)-RSINX(1))
!        DRSINX(2) = 1.D+6*(RSINX(3)-RSINX(1))
!        DCA_GNX = (DRSINX(1)*RCOSX(1) - DRCOSX(1)*RSINX(1))/
!     &    (DRSINX(2)*DRCOSX(1)-DRSINX(1)*DRCOSX(2))
!        DCA_NLX = (DRCOSX(2)*RSINX(1) - DRSINX(2)*RCOSX(1))/
!     &    (DRSINX(2)*DRCOSX(1)-DRSINX(1)*DRCOSX(2))
!        CA_GN = CA_GN + DCA_GNX
!        CA_NL = CA_NL + DCA_NLX
!        NC = NC + 1
!        IF( NC.GT.32 ) THEN
!          INDX = 3
!          CHMSG = 'Unconverged Interfacial Contact Angle:'
!          CALL WRMSGS( INDX )
!        ENDIF
!        IF( ABS(DCA_GNX).GT.1.D-9 .OR. ABS(DCA_NLX).GT.1.D-9 ) GOTO 20
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
  100 CONTINUE
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
        ISNR = 0
        IF( INDEX(ADUM(1:),'entrap').NE.0 .AND.
     &    INDEX(ADUM(1:),'resid').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            ISCHRX = 34
            ISNR = 5
          ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &      INDEX(ADUM(1:),'corey').NE.0 ) THEN
            ISCHRX = 38
            ISNR = 5
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
            ISNR = 5
          ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &      INDEX(ADUM(1:),'corey').NE.0 ) THEN
            ISCHRX = 37
            ISNR = 5
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
!
!---    Saturation function indices  ---
!
        IF( IJK.GT.0 ) THEN
          DO 330 N = 1,NFLD
            ISCHR(IZ(N)) = ISCHRX
  330     CONTINUE
        ELSE
          ISCHR(IROCK) = ISCHRX
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
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &        ', ',UNTS(1:NCH),': ',SCHR(1,IROCK)
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,SCHR(1,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(1,IROCK),', 1/m)'
          ENDIF
          VARB = 'van Genuchten (n)'
          IF( IJK.GT.0 ) THEN
            INDX = 3
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(3,IROCK))
            WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHR(3,IROCK)
          ENDIF
          VARB = 'van Genuchten (minimum saturation)'
          IF( IJK.GT.0 ) THEN
            INDX = 4
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(4,IROCK))
            WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHR(4,IROCK)
          ENDIF
          VARB = 'van Genuchten (Reference Interfacial Tension)'
          IF( IJK.GT.0 ) THEN
            INDX = 8
            LNDX = LSCHR
            UNTS = 'n/m'
            IUNKG = 1
            IUNS = -2
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(8,IROCK))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &        ', ',UNTS(1:NCH),': ',SCHR(8,IROCK)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHR(8,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(8,IROCK),', N/m)'
          ENDIF
          IF( IJK.GT.0 ) THEN
            DO 340 N = 1,NFLD
              SCHR(2,IZ(N)) = 0.D+0
              SCHR(9,IZ(N)) = 1.D+20
  340       CONTINUE
          ELSE
            SCHR(2,IROCK) = 0.D+0
            SCHR(9,IROCK) = 1.D+20
          ENDIF
          VARB = 'van Genuchten (m)'
          IF( IJK.GT.0 ) THEN
            INDX = 14
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(14,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(14,IROCK)
          ENDIF
          IF( ISCHRX.EQ.32 ) THEN
            VARB = 'van Genuchten (Effective Maximum Trapped NAPL ' //
     &        'Saturation)'
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
            VARB = 'van Genuchten (Critical Trapping Number)'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              INDX = 9
              LNDX = LSCHR
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(9,IROCK))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ', 
     &          SCHR(9,IROCK)
            ENDIF
            IF( IJK.GT.0 ) THEN
              DO 350 N = 1,NFLD
                IF( SCHR(9,IZ(N))/EPSL.LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative or Zero Critical Trapping Number'
                  CALL WRMSGS( INDX )
                ENDIF
  350         CONTINUE
            ELSE
              IF( SCHR(9,IROCK)/EPSL.LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative or Zero Critical Trapping Number'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( ISCHRX.EQ.33 ) THEN
            VARB = 'van Genuchten (Effective Maximum Residual NAPL ' //
     &        'Saturation)'
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
          ELSEIF( ISCHRX.EQ.34 ) THEN
            VARB = 'van Genuchten (Effective Maximum Trapped NAPL ' //
     &        'Saturation)'
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
            VARB = 'van Genuchten (Effective Maximum Residual NAPL ' //
     &        'Saturation)'
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
            VARB = 'van Genuchten (Critical Trapping Number)'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              INDX = 9
              LNDX = LSCHR
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(9,IROCK))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ', 
     &          SCHR(9,IROCK)
            ENDIF
            IF( IJK.GT.0 ) THEN
              DO 360 N = 1,NFLD
                IF( SCHR(9,IZ(N))/EPSL.LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative or Zero Critical Trapping Number'
                  CALL WRMSGS( INDX )
                ENDIF
  360         CONTINUE
            ELSE
              IF( SCHR(9,IROCK)/EPSL.LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative or Zero Critical Trapping Number'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
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
          IF( IJK.GT.0 ) THEN
            INDX = 2
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
          VARB = 'Brooks-Corey (lambda parameter)'
          IF( IJK.GT.0 ) THEN
            INDX = 3
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(3,IROCK))
            WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHR(3,IROCK)
          ENDIF
          VARB = 'Brooks-Corey (residual saturation)'
          IF( IJK.GT.0 ) THEN
            INDX = 4
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(4,IROCK))
            WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHR(4,IROCK)
          ENDIF
          VARB = 'Brooks-Corey (Reference Interfacial Tension)'
          IF( IJK.GT.0 ) THEN
            INDX = 8
            LNDX = LSCHR
            UNTS = 'n/m'
            IUNKG = 1
            IUNS = -2
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(8,IROCK))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHR(8,IROCK)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHR(8,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(8,IROCK),', N/m)'
          ENDIF
          IF( IJK.GT.0 ) THEN
            DO 370 N = 1,NFLD
              SCHR(9,IZ(N)) = 1.D+20
  370       CONTINUE
          ELSE
            SCHR(9,IROCK) = 1.D+20
          ENDIF
          IF( ISCHRX.EQ.36 ) THEN
            VARB = 'Brooks-Corey (Effective Maximum Trapped NAPL ' //
     &        'Saturation)'
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
            VARB = 'Brooks-Corey (Critical Trapping Number)'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              INDX = 9
              LNDX = LSCHR
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(9,IROCK))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ', 
     &          SCHR(9,IROCK)
            ENDIF
            IF( IJK.GT.0 ) THEN
              DO 380 N = 1,NFLD
                IF( SCHR(9,IZ(N))/EPSL.LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative or Zero Critical Trapping Number'
                  CALL WRMSGS( INDX )
                ENDIF
  380         CONTINUE
            ELSE
              IF( SCHR(9,IROCK)/EPSL.LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative or Zero Critical Trapping Number'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( ISCHRX.EQ.37 ) THEN
            VARB = 'Brooks-Corey (Effective Maximum Resdiual NAPL ' //
     &        'Saturation)'
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
          ELSEIF( ISCHRX.EQ.38 ) THEN
            VARB = 'Brooks-Corey (Effective Maximum Trapped NAPL ' //
     &        'Saturation)'
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
            VARB = 'Brooks-Corey (Effective Maximum Residual NAPL ' //
     &        'Saturation)'
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
            VARB = 'Brooks-Corey (Critical Trapping Number)'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              INDX = 9
              LNDX = LSCHR
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(9,IROCK))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ', 
     &          SCHR(9,IROCK)
            ENDIF
            IF( IJK.GT.0 ) THEN
              DO 390 N = 1,NFLD
                IF( SCHR(9,IZ(N))/EPSL.LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative or Zero Critical Trapping Number'
                  CALL WRMSGS( INDX )
                ENDIF
  390         CONTINUE
            ELSE
              IF( SCHR(9,IROCK)/EPSL.LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative or Zero Critical Trapping Number'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ENDIF
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
!---  End of RDSP4 group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSR4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Water-Oil Mode
!
!     Read input file for source information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 VAR(LSTM,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSR4'
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
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'napl').NE.0 .AND. 
     &      INDEX(ADUM(1:),'spill').NE.0 .AND.
     &      INDEX(ADUM(1:),'volu').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'NAPL Spill Volumetric Source'
          ISRTX = 40
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &          INDEX(ADUM(1:),'spill').NE.0 .AND.
     &          INDEX(ADUM(1:),'volu').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Spill Volumetric Source w/ Mass 
     &Frac.'
          ISRTX = 41
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &          INDEX(ADUM(1:),'volu').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Volumetric Source w/ Mass Frac.'
          ISRTX = 3
        ELSEIF( INDEX(ADUM(1:),'napl').NE.0 .AND.
     &          INDEX(ADUM(1:),'volu').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'NAPL Volumetric Source'
          ISRTX = 6
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &          INDEX(ADUM(1:),'liquid').NE.0 .AND.
     &          INDEX(ADUM(1:),'volu').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Total Liquid Volumetric Source'
          ISRTX = 4
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &          INDEX(ADUM(1:),'mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Mass Source w/ Mass Frac.'
          ISRTX = 7
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &          INDEX(ADUM(1:),'liquid').NE.0 .AND.
     &          INDEX(ADUM(1:),'mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Total Liquid Mass Source'
          ISRTX = 8
        ELSEIF( INDEX(ADUM(1:),'napl').NE.0 .AND.
     &          INDEX(ADUM(1:),'mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'NAPL Mass Source'
          ISRTX = 10
#ifdef awtap
        ELSEIF( INDEX(ADUM(1:),'well').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'flow').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Flow Well'
            ISRTX = 20
          ELSEIF( INDEX(ADUM(1:),'head').NE.0 .OR.
     &      INDEX(ADUM(1:),'elevation').NE.0 .OR.
     &      INDEX(ADUM(1:),'depth').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Elevation Head Well'
            ISRTX = 21
          ENDIF
#endif
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
!---  Read source domain indices  ---
!
        IF( ISRTX.EQ.40 .OR. ISRTX.EQ.41 ) THEN
          VARB = 'Source Domain Index'
          ISX = ISTART
          CALL RDINT(ISTART,ICOMMA,CHDUM,I1X)
          CALL RDINT(ISTART,ICOMMA,CHDUM,I2X)
          CALL RDINT(ISTART,ICOMMA,CHDUM,J1X)
          CALL RDINT(ISTART,ICOMMA,CHDUM,J2X)
          ICX = ISTART
          WRITE(IWR,'(/,2X,A)') 'Source Domain:'
          WRITE(IWR,'(4X,A,I6,A,I6)') 'I = ',I1X,' to ',I2X
          WRITE(IWR,'(4X,A,I6,A,I6)') 'J = ',J1X,' to ',J2X
          K1X=1
          K2X=1
        ELSE
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
        ENDIF
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
     &     K2X.GT.KFLD .OR. K2X.LT.K1X ) THEN
          INDX = 4
          CHMSG = 'Invalid Source Domain: ' // CHDUM(ISX:ICX)
          CALL WRMSGS( INDX )
        ENDIF
#ifdef awtap
!
!---  Check for z-direction well sources  ---
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
                ELSEIF( K2X.GE.ISRDM(5,M) .AND. K2X.LE.ISRDM(6,M) ) THEN
                  INDX = 4
                  CHMSG = 'Overlapping Well-Source Domains'
                  CALL WRMSGS( INDX )
                ENDIF
              ENDIF
            END IF
   45     CONTINUE
        ENDIF
#endif
!
!---  Check for sources applied to inactive nodes  ---
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
          WRITE(IWR,'(2X,4A,1PE11.4)')VARB(1:IVR),', ',UNTS(1:NCH),': ',
     &      VAR(NTM,1)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,VAR(NTM,1),INDX)
          IF( ISRTX.EQ.3 ) THEN
            VARB = 'Source Aqueous Volumetric Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Mass Fraction Dissolved Oil'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
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
          ELSEIF( ISRTX.EQ.4 ) THEN
            VARB = 'Source Total Liquid Volumetric Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Volume Fraction Oil'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
          ELSEIF( ISRTX.EQ.7 ) THEN
            VARB = 'Source Aqueous Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Mass Fraction Dissolved Oil'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
          ELSEIF( ISRTX.EQ.8 ) THEN
            VARB = 'Source Total Liquid Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Mass Fraction Oil'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
          ELSEIF( ISRTX.EQ.10 ) THEN
            VARB = 'Source NAPL Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
          ELSEIF( ISRTX.EQ.40 ) THEN
            VARB = 'Source NAPL Spill Volumetric Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
          ELSEIF( ISRTX.EQ.41 ) THEN
            VARB = 'Source Aqueous Spill Volumetric Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
#ifdef future_coding
          ELSEIF( ISRTX.EQ.20 ) THEN
            VARB = 'Source Volumetric Flow Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,4),', m^3/s)'
            VARB = 'Source Well Bore-Hole Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,3),', m)'
            IF( VAR(NTM,4).LT.ZERO ) THEN
              VAR(NTM,2) = 1.D+0
            ELSE
              VAR(NTM,2) = 100.D+0
            ENDIF
            IDFLT = 1
            VARB = 'Source Well Liquid Level Constraint'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,2)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,2),', m)'
            VARB = 'Source Dissolved-Oil Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
#endif
#ifdef awtap
          ELSEIF( ISRTX.EQ.21 ) THEN
            VARB = 'Source Liquid Height in Well'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,4),', m)'
            VARB = 'Source Well Bore-Hole Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,3),', m)'
            VARB = 'Source Dissolved-Oil Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
#endif
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
        DO 130 NTM = 1,ISRM(NS)
          DO 120 M = 1,6
            SRC(M,NTM,NSR) = VAR(NTM,M)
  120     CONTINUE
  130   CONTINUE
  140 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSR4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTF4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Reads solute/fluid interaction card for diffusion and partition
!     coefficients, and internodal diffusion term averaging scheme for
!     two phase (aqueous-NAPL) convective-dispersive transport equation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
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
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICDSX
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDTF4'
!
!---  Write card information to output file  ---
!
      CARD = 'Solute/Fluid Interaction Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
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
        CALL LCASE(CHDUM)
        ISTART = 1
        ADUM(1:) = ' '
        VARB = 'Solute Name'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        DO 100 NS = 1,NSOLU
          IF( SOLUT(NS).EQ.ADUM ) THEN
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
        WRITE(IWR,'(/,3A)') VARB(1:IVR),': ',ADUM
!
!---  Aqueous-phase molecular diffusion coefficient  ---
!
        VARB = 'Aqueous-Phase Molecular Diffusion Coefficient'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SMDL(NSL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),': '
     &,SMDL(NSL)
        INDX = 0
        IUNM = 2
        IUNS = -1
        CALL RDUNIT(UNTS,SMDL(NSL),INDX)
!
!---  NAPL-phase molecular diffusion coefficient  ---
!
        VARB = 'NAPL Molecular Diffusion Coefficient'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SMDN(NSL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),': '
     &,SMDN(NSL)
        INDX = 0
        IUNM = 2
        IUNS = -1
        CALL RDUNIT(UNTS,SMDN(NSL),INDX)
!
!---  Aqueous-NAPL partition coefficient  ---
!
!        IDFLT = 1
!        VARB = 'Aqueous-NAPL Partition Coefficient'
!        CALL RDDPR(ISTART,ICOMMA,CHDUM,PCLN(1,NSL))
!        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
!        WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),': '
!     &,PCLN(1,NSL)
!        INDX = 0
!        CALL RDUNIT(UNTS,PCLN(1,NSL),INDX)
!        PCLN(1,NSL) = MAX( PCLN(1,NSL),1.D-20 )
!        PCLN(1,NSL) = MIN( PCLN(1,NSL),1.D+20 )
!
!---  Read aqueous-NAPL adsorption function  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Aqueous-NAPL Adsorption Function'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'linear').NE.0 ) THEN
          IPCLN(NSL) = 0
          WRITE(IWR,'(A)') 'Linear Isotherm'
          VARB = 'k parameter'
          IDFLT = 1
          PCLN(1,NSL) = 1.D+20
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCLN(1,NSL))
          WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCLN(1,NSL)
          PCLN(1,NSL) = MAX( PCLN(1,NSL),1.D-20 )
        ELSEIF( INDEX(ADUM(1:),'freundlich').NE.0 ) THEN
          IPCLN(NSL) = 2
          WRITE(IWR,'(A)') 'Freundlich Isotherm'
          VARB = 'k parameter'
          IDFLT = 1
          PCLN(1,NSL) = 1.D+20
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCLN(1,NSL))
          WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCLN(1,NSL)
          PCLN(1,NSL) = MAX( PCLN(1,NSL),1.D-20 )
          VARB = 'n Parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCLN(2,NSL))
          WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCLN(2,NSL)
        ELSEIF( INDEX(ADUM(1:),'langmuir').NE.0 ) THEN
          IPCLN(NSL) = 3
          WRITE(IWR,'(A)') 'Langmuir Isotherm'
          VARB = 'a parameter'
          IDFLT = 1
          PCLN(1,NSL) = 1.D+20
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCLN(1,NSL))
          WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCLN(1,NSL)
          PCLN(1,NSL) = MAX( PCLN(1,NSL),1.D-20 )
          VARB = 'b parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCLN(2,NSL))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',PCLN(2,NSL)
          INDX = 0
          IUNM = 3
          CALL RDUNIT(UNTS,PCLN(2,NSL),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PCLN(2,NSL),', m^3)'
        ENDIF
#ifdef radioactive
!
!---  Half-life  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
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
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
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
        CALL LCASE(CHDUM)
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDTF4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTP4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Reads the solute/porous media interaction card for the
!     dispersivities, species half-lives, and partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/RDTP4'
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
        VARB = 'Longitudinal Dispersivity'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DISPL(IROCK))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),', ',
     &DISPL(IROCK)
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
        WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),', ',
     &DISPT(IROCK)
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
          CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
!
!---  Search known solutes for matching name  ---
!
          DO 300 NSS = 1,NSOLU
            IF( ADUM.EQ.SOLUT(NSS)) THEN
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
!          IDFLT = 1
!          VARB = 'Solid-Aqueous Partition Coefficient'
!          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
!          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
!          WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),': ',
!     &      PCSL(1,IROCK,NSL)
!          INDX = 0
!          IUNM = 3
!          IUNKG = -1
!          CALL RDUNIT(UNTS,PCSL(1,IROCK,NSL),INDX)
!          PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-12 )
!
!---  Solid-aqueous partition coefficient  ---
!
        VARB = 'Solid-aqueous Adsorption Function'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'linear kd').NE.0 ) THEN
          IPCSL(IROCK,NSL) = 1
          WRITE(IWR,'(A)') 'Linear Isotherm w/ Kd'
          VARB = 'Kd parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),': ',
     &    PCSL(1,IROCK,NSL)
          INDX = 0
          IUNM = 3
          IUNKG = -1
          CALL RDUNIT(UNTS,PCSL(1,IROCK,NSL),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PCSL(1,IROCK,NSL),', m^3/kg)'
          PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
        ELSEIF( INDEX(ADUM(1:),'linear').NE.0 ) THEN
          IPCSL(IROCK,NSL) = 0
          WRITE(IWR,'(A)') 'Linear Isotherm'
          VARB = 'k parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
          WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCSL(1,IROCK,NSL)
          PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
        ELSEIF( INDEX(ADUM(1:),'freundlich').NE.0 ) THEN
          IPCSL(IROCK,NSL) = 2
          WRITE(IWR,'(A)') 'Freundlich Isotherm'
          VARB = 'k parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
          WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCSL(1,IROCK,NSL)
          PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
          VARB = 'n Parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(2,IROCK,NSL))
          WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCSL(2,IROCK,NSL)
        ELSEIF( INDEX(ADUM(1:),'langmuir').NE.0 ) THEN
          IPCSL(IROCK,NSL) = 3
          WRITE(IWR,'(A)') 'Langmuir Isotherm'
          VARB = 'a parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
          WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCSL(1,IROCK,NSL)
          PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
          VARB = 'b parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(2,IROCK,NSL))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',PCSL(2,IROCK,NSL)
          INDX = 0
          IUNM = 3
          CALL RDUNIT(UNTS,PCSL(2,IROCK,NSL),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PCSL(2,IROCK,NSL),', m^3)'
        ENDIF
  500   CONTINUE
!
!---  Read next rock/soil type or scaling group  ---
!
      IF( N.LT.NROCK ) WRITE(IWR,'(/)')
      GOTO 10
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDTP4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Compute the maximum relative residuals
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
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
      SUB_LOG(ISUB_LOG) = '/RSDL4'
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
        MPL = IM(IEQW,NMD)
        MPN = IM(IEQO,NMD)
        IF( ISKP(IZN).EQ.1 ) GOTO 200
!
!---    Aqueous or aqueous-gas w/ dissolved-oil system
!       prior to iteration  ---
!
        IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.2 ) THEN
          RSDX = ABS(BLU(MPL))/(ABS(PL(2,N))+PATM)
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          RSDX = ABS(BLU(MPN))/MAX( PSO(2,N)/HCOW,1.D-12 )
          RSDX = 0.D+0
          IF( RSDX.GT.RSD(IEQO) ) THEN
            RSD(IEQO) = RSDX
            NSD(IEQO) = N
          ENDIF
!
!---    Aqueous or aqueous-gas w/ trapped-oil system
!       prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.4 ) THEN
          RSDX = ABS(BLU(MPL))/(ABS(PL(2,N))+PATM)
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          RSDX = ABS(BLU(MPN))
          IF( RSDX.GT.RSD(IEQO) ) THEN
            RSD(IEQO) = RSDX
            NSD(IEQO) = N
          ENDIF
!
!---    Aqueous-oil or aqueous-oil-gas system
!       prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 .OR. NPHAZ(2,N).EQ.6 ) THEN
          ACPL = PORD(2,N)*RHOL(2,N)*SL(2,N)*DTI*VOL(N)+SMALL
          RSDX = ABS(RSDL(IEQW,N)/ACPL)
          RSDX = MIN( RSDX, ABS(BLU(MPL))/(ABS(PL(2,N))+PATM) )
          IF( 1.D+5-(PN(2,N)-PL(2,N)).LE.EPSL ) RSDX = 0.D+0
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          ACPN = PORD(2,N)*RHON(2,N)*SN(2,N)*DTI*VOL(N)+SMALL
          RSDX = ABS(RSDL(IEQO,N)/ACPN)
          RSDX = MIN( RSDX, ABS(BLU(MPN))/(ABS(PN(2,N))+PATM) )
          IF( RSDX.GT.RSD(IEQO) ) THEN
            RSD(IEQO) = RSDX
            NSD(IEQO) = N
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
        WRITE(ISC,'(A)') 'Execution Note: Newton Convergence Failure'
        WRITE(IWR,'(A)') 'Execution Note: Newton Convergence Failure'
        WRITE(ISC,'(A,1PE11.4,$)') '  Water Residual: ',RSD(IEQW)
        WRITE(ISC,'(A,I6,$)') '  Node: ',NSD(IEQW)
        WRITE(ISC,'(A,I1,$)') '  Previous Phase: ',NPHAZ(1,NSD(IEQW))
        WRITE(ISC,'(A,I1)') '  Current Phase: ',NPHAZ(2,NSD(IEQW))
        IF( NSD(IEQO).GT.0 ) THEN
          WRITE(ISC,'(A,1PE11.4,$)') '  Oil Residual: ',RSD(IEQO)
          WRITE(ISC,'(A,I6,$)') ' Node: ',NSD(IEQO)
          WRITE(ISC,'(A,I1,$)') '  Previous Phase: ',NPHAZ(1,NSD(IEQO))
          WRITE(ISC,'(A,I1)') '  Current Phase: ',NPHAZ(2,NSD(IEQO))
        ENDIF
        WRITE(IWR,'(A,1PE11.4,$)') '  Water Residual: ',RSD(IEQW)
        WRITE(IWR,'(A,I6,$)') '  Node: ',NSD(IEQW)
        WRITE(IWR,'(A,I1,$)') '  Previous Phase: ',NPHAZ(1,NSD(IEQW))
        WRITE(IWR,'(A,I1)') '  Current Phase: ',NPHAZ(2,NSD(IEQW))
        IF( NSD(IEQO).GT.0 ) THEN
          WRITE(IWR,'(A,1PE11.4,$)') '  Oil Residual: ',RSD(IEQO)
          WRITE(IWR,'(A,I6,$)') '  Node: ',NSD(IEQO)
          WRITE(IWR,'(A,I1,$)') '  Previous Phase: ',NPHAZ(1,NSD(IEQO))
          WRITE(IWR,'(A,I1)') '  Current Phase: ',NPHAZ(2,NSD(IEQO))
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
       WRITE(ISC,'(A,1PE11.4,1X,$)') '  Time step reduced from ',VARX
       WRITE(ISC,'(2A,1PE11.4,1X,A)') UNTM(1:NCH),' to ',VAR,UNTM(1:NCH)
       WRITE(IWR,'(A,1PE11.4,1X,$)') '  Time step reduced from ',VARX
       WRITE(IWR,'(2A,1PE11.4,1X,A)') UNTM(1:NCH),' to ',VAR,UNTM(1:NCH)
          DO 400 N = 1,NFLD
            PL(2,N) = PL(1,N)
            PN(2,N) = PN(1,N)
            XMLO(2,N) = XMLO(1,N)
            SN(2,N) = SN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  400     CONTINUE
#ifdef awtap
          DO 402 N = 1,NSR
            PLWB(2,N) = PLWB(1,N)
            PGW(2,N) = PGW(1,N)
  402     CONTINUE
#endif
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
          WRITE(ISC,'(A)') 'Execution Note:  Reduction Limit Exceeded'
          WRITE(IWR,'(A)') 'Execution Note:  Reduction Limit Exceeded'
          ICNV = 4
        ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDL4 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SBND4( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Modify the Jacobian matrix for the solute transport equation
!     to incorporate boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1994.
!     Last Modified by M. D. White, Battelle, PNL, October 16, 1995.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE REACT
      USE PORMED
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FDVP
      USE FDVN
      USE CONST
      USE BCVP
      USE BCVN
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
      REAL*8 BCX(LSPBC+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SBND4'
!
!---  Loop over number of specified boundary conditions  ---
!
      NBCT = MIN( NSL+LUK,NSOLU+LUK+1 )
      DO 200 NB = 1,NBC
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 200
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) GOTO 200
#ifdef ecke
!
!---    Solute transport  ---
!
        IF( NSL.LE.NSOLU ) THEN
#endif
          IBCTX = IBCT(NSL+LUK,NB)
#ifdef ecke
!
!---    Reactive species transport  ---
!
        ELSE
          IBCTX = IBCT(NSOLU+LUK+1,NB)
        ENDIF
#endif
!
!---    Zero flux boundary condition  ---
!
        IF( IBCTX.EQ.3 ) GOTO 200
!
!---    Single boundary condition time  ---
!
        IF( IBCM(NB).EQ.1 ) THEN
#ifdef ecke
!
!---      Solute transport  ---
!
          IF( NSL.LE.NSOLU ) THEN
#endif
            BCX(1) = BC(NSL+LBCU,1,MB)
            IF( IBCTX.EQ.12 ) BCX(1) = CBO(NB,NSL)
#ifdef ecke
!
!---      Reactive species transport  ---
!
          ELSE
            BCX(1) = 0.D+0
            DO 10 NSPX = 1,IBCSP(1,NB)
              NSP = IBCSP(NSPX+1,NB)
              MX = NSOLU+LBCU+NSPX
              BCX(NSPX+1) = BC(MX,1,MB)
!
!---          Aqueous species ---
!
              IF( NSP.LE.NSPL ) THEN
                IF( IBCT(NSOLU+LUK+1,NB).EQ.12 ) 
     &            BCX(NSPX+1) = SP_CBO(NB,NSP)
!
!---          NAPL species ---
!
              ELSE
                IF( IBCT(NSOLU+LUK+2,NB).EQ.12 ) 
     &            BCX(NSPX+1) = SP_CBO(NB,NSP)
              ENDIF
   10       CONTINUE
          ENDIF
#endif
!
!---    Multiple boundary condition times  ---
!
        ELSE
          DO 100 M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
             TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
             DTBC = MIN( BC(1,M,MB)-TMZ,DT )
             TFBC = (TMZ-5.D-1*DTBC-BC(1,M-1,MB))/TDBC
#ifdef ecke
!
!---          Solute transport  ---
!
             IF( NSL.LE.NSOLU ) THEN
#endif
               BCX(1) = BC(NSL+LBCU,M-1,MB) +
     &           TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
               IF( IBCT(NSL+LUK,NB).EQ.12 ) BCX(1) = CBO(NB,NSL)
#ifdef ecke
!
!---          Reactive species transport  ---
!
              ELSE
                BCX(1) = 0.D+0
                DO 20 NSPX = 1,IBCSP(1,NB)
                  NSP = IBCSP(NSPX+1,NB)
                  MX = NSOLU+LBCU+NSPX
                  BCX(NSPX+1) = BC(MX,M-1,MB) +
     &              TFBC*(BC(MX,M,MB)-BC(MX,M-1,MB))
!
!---              Aqueous species ---
!
                  IF( NSP.LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.12 ) 
     &                BCX(NSPX+1) = SP_CBO(NB,NSP)
!
!---              NAPL species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.12 ) 
     &                BCX(NSPX+1) = SP_CBO(NB,NSP)
                  ENDIF
   20           CONTINUE
              ENDIF
#endif
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
!-----  Diffusion coefficients at node adjacent to boundary  ---
!
        SVLP = PORD(2,N)*SL(2,N)
        SVNP = PORD(2,N)*SN(2,N)
        TCOR = (T(2,N)+TABS)/TSPRF
        SMDLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
        DLP = TORL(2,N)*SVLP*SMDLP
        SMDNP = SMDN(NSL)*TCOR*(VISRN/VISN(2,N))
        DNP = TORN(2,N)*SVNP*SMDNP
!
!-----  Phase fraction factors at node adjacent to boundary  ---
!
        FCLP = YL(N,NSL)/(SVLP+SMALL)
        FCNP = YN(N,NSL)/(SVNP+SMALL)
!
!-----  Phase fraction factors at boundary  ---
!
        XVLB = PORDB(2,NB)*SLB(2,NB)
        XVNB = PORDB(2,NB)*SNB(2,NB)
#ifdef ecke
        IF( NSL.LE.NSOLU ) THEN
#endif
!
!-----    Linear Kd sorption  ---
!
          IF( IPCSL(IZN,NSL).EQ.1 ) THEN
            PCSLX = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
!
!-----    Linear, Freundlich, or Langmuir sorption  ---
!
          ELSE
            PCSLX = PCSL(1,IZN,NSL)
          ENDIF
!
!-----    Aqueous-NAPL adsorption linear isotherm and
!         Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
          CTX = 1.D+0
          IF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
            CALL LN_LN( CTX,CNX,CLX,XVLB,XVNB,N,IZN,NSL )
!
!-----    Aqueous-NAPL adsorption Freundlich isotherm and
!         Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
          ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
            CALL FR_LN( CTX,CNX,CLX,XVLB,XVNB,N,IZN,NSL )
!
!-----    Aqueous-NAPL adsorption Langmuir isotherm and
!         Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
          ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
            CALL LG_LN( CTX,CNX,CLX,XVLB,XVNB,N,IZN,NSL )
!
!-----    Aqueous-NAPL adsorption linear isotherm and
!         Solid-Aqueous adsorption Freundlich isotherm  ---
!
          ELSEIF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
            CALL LN_FR( CTX,CNX,CLX,XVLB,XVNB,N,IZN,NSL )
!
!-----    Aqueous-NAPL adsorption Freundlich isotherm and
!         Solid-Aqueous adsorption Freundlich isotherm  ---
!
          ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
            CALL FR_FR( CTX,CNX,CLX,XVLB,XVNB,N,IZN,NSL )
!
!-----    Aqueous-NAPL adsorption Langmuir isotherm and
!         Solid-Aqueous adsorption Fruendlich isotherm  ---
!
          ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
            CALL LG_FR( CTX,CNX,CLX,XVLB,XVNB,N,IZN,NSL )
!
!-----    Aqueous-NAPL adsorption linear isotherm and
!         Solid-Aqueous adsorption Langmuir isotherm  ---
!
          ELSEIF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
            CALL LN_LG( CTX,CNX,CLX,XVLB,XVNB,N,IZN,NSL )
!
!-----    Aqueous-NAPL adsorption Freundlich isotherm and
!         Solid-Aqueous adsorption Langmuir isotherm  ---
!
          ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
            CALL FR_LG( CTX,CNX,CLX,XVLB,XVNB,N,IZN,NSL )
!
!-----    Aqueous-NAPL adsorption Langmuir isotherm and
!         Solid-Aqueous adsorption Langmuir isotherm  ---
!
          ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
            CALL LG_LG( CTX,CNX,CLX,XVLB,XVNB,N,IZN,NSL )
          ENDIF
          YLB(NB,NSL) = XVLB*CLX/CTX
          YNB(NB,NSL) = XVNB*CNX/CTX
          FCL = YLB(NB,NSL)/(XVLB+SMALL)
          FCN = YNB(NB,NSL)/(XVNB+SMALL)
!
!-----    Convert boundary phase concentrations to
!         volumetric concentrations  ---
!
          IF( IBCT(NSL+LUK,NB).EQ.8 .OR. 
     &      IBCT(NSL+LUK,NB).EQ.14 ) THEN
            BCX(1) = BCX(1)/FCL
          ELSEIF( IBCT(NSL+LUK,NB).EQ.10 .OR. 
     &      IBCT(NSL+LUK,NB).EQ.16 ) THEN
            BCX(1) = BCX(1)/FCN
          ENDIF
          CB(NB,NSL) = BCX(1)
#ifdef ecke
        ELSE
!
!---      Convert species concentrations to total-component
!         concentrations  ---
!
          IF( NSL.LE.NSOLU+NEQC ) THEN
            NEQ = NSL-NSOLU
            YSPLX = 0.D+0
            YSPNX = 0.D+0
            DO 130 NSP = 1,IEQ_C(1,NEQ)
              DO 120 NSPX = 1,IBCSP(1,NB)
                IF( IBCSP(NSPX+1,NB).EQ.IEQ_C(NSP+1,NEQ) ) THEN
!
!---              Aqueous species ---
!
                  IF( IEQ_C(NSP+1,NEQ).LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.14 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+1,NB).EQ.10 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.16 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVNB
                    ENDIF
                    YSPLX = YSPLX + EQ_C(NSP,NEQ)*BCX(NSPX+1)
!
!---              NAPL species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.14 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+2,NB).EQ.10 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.16 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVNB
                    ENDIF
                    YSPNX = YSPNX + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                  ENDIF                    
                  BCX(1) = BCX(1) + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
  120         CONTINUE
  130       CONTINUE
!
!---        Linked aqueous oil   ---
!
            IF( ISPLK(12).EQ.NSL ) BCX(1) = 1.D+3*XLOB(2,NB)*
     &          RHOLB(2,N)*SLB(2,NB)*PORDB(2,NB)/WTMO
!
!---        Convert species concentrations to total-kinetic
!           concentrations  ---
!
          ELSEIF( NSL.LE.NSOLU+NEQC+NEQK ) THEN
            NEQ = NSL-NSOLU-NEQC
            YSPLX = 0.D+0
            DO 150 NSP = 1,IEQ_K(1,NEQ)
              DO 140 NSPX = 1,IBCSP(1,NB)
                IF( IBCSP(NSPX+1,NB).EQ.IEQ_K(NSP+1,NEQ) ) THEN
!
!---              Aqueous species ---
!
                  IF( IEQ_K(NSP+1,NEQ).LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.14 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+1,NB).EQ.10 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.16 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVNB
                    ENDIF
                    YSPLX = YSPLX + EQ_K(NSP,NEQ)*BCX(NSPX+1)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.14 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+2,NB).EQ.10 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.16 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVNB
                    ENDIF
                  ENDIF                    
                  BCX(1) = BCX(1) + EQ_K(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
  140         CONTINUE
  150       CONTINUE
!
!---        Linked aqueous oil   ---
!
            IF( ISPLK(12).EQ.NSL ) BCX(1) = 1.D+3*XLOB(2,NB)*
     &          RHOLB(2,N)*SLB(2,NB)*PORDB(2,NB)/WTMO
          ENDIF
          IF( ABS(BCX(1))/EPSL.LT.EPSL ) THEN
            YSPLX = 0.D+0
          ELSE
            YSPLX = YSPLX/BCX(1)
          ENDIF
!
!---      Phase-volumetric concentration ratios  ---
!
          YLBX = MAX( MIN( 1.D+0,YSPLX ),0.D+0 )
          YNBX = 1.D+0-YLBX
          FCL = 0.D+0
          IF( XVLB/EPSL.GT.EPSL ) FCL = YLBX/XVLB
          FCN = 0.D+0
          IF( XVNB/EPSL.GT.EPSL ) FCN = YNBX/XVNB
!
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = YLBX
          YNB(NB,NSL) = YNBX
          CB(NB,NSL) = BCX(1)
        ENDIF
#endif
!
!-----  Bottom boundary  ---
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
          ELSE
            DPLB = 0.D+0
            DPNB = 0.D+0
          ENDIF
          FLB = AFZ(NPZ)*WL(1,NPZ)
          FNB = AFZ(NPZ)*WN(1,NPZ)
          CRLB = ABS( WL(1,NPZ) )*DT/(DZGF(N)*XVLB+SMALL)
          CRNB = ABS( WN(1,NPZ) )*DT/(DZGF(N)*XVNB+SMALL)
!
!-------  Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.10 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*XVLB*SMDLB
            INDX = 16
            DLZ = DIFMN(DLB,DLP,DZGF(N),DZGF(N),WL(1,NPZ),INDX)
            DLZ = AFZ(NPZ)*(DLZ+DPLB)/(5.D-1*DZGF(N))
            PCOR = (PNB(2,NB)+PATM)/PATM
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*XVNB*SMDNB
            INDX = 16
            DNZ = DIFMN(DNB,DNP,DZGF(N),DZGF(N),WN(1,NPZ),INDX)
            DNZ = AFZ(NPZ)*(DNZ+DPNB)/(5.D-1*DZGF(N))
!
!---------  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              IF( FLB.GE.ZERO ) THEN
                WCLZ = BCX(1)*FCL*FLB
              ELSEIF( FLB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCLT = YL(NBT,NSL)/(SL(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBT,NSL)*FCLT)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                WCLZ = BCX(1)*FLB*THETA*FCL
     &            + C(N,NSL)*FLB*(1.D+0-THETA)*FCLP
              ELSEIF( FLB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCLZ = C(N,NSL)*FLB*FCLP
              ENDIF
              IF( FNB.GE.ZERO ) THEN
                WCNZ = BCX(1)*FCN*FNB
              ELSEIF( FNB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCNT = YN(NBT,NSL)/(SN(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBT,NSL)*FCNT)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRNB,ISLC(1) )
                WCNZ = BCX(1)*FNB*THETA*FCN
     &            + C(N,NSL)*FNB*(1.D+0-THETA)*FCNP
              ELSEIF( FNB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCNZ = C(N,NSL)*FNB*FCNP
              ENDIF
              AB = DLZ*FCL + DNZ*FCN
              AP = DNZ*FCNP + DLZ*FCLP
              WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ+WCNZ)/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCNZ
!
!-----------  TVD Transport for interior surface 
!             adjacent to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVLX = SL(2,NBT)*PORD(2,NBT)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(SVLX+SMALL)
                FCLT = YL(NBT,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNT = AFZ(NQZ)*WN(1,NQZ)
              IF( FNT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVNX = SN(2,NBT)*PORD(2,NBT)
                CRNT = ABS( WN(1,NQZ) )*DT/DZGP(NQZ)/(SVNX+SMALL)
                FCNT = YN(NBT,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANB = MAX( FNB,ZERO ) +
     &          DNZ*MAX((ONE-(TENTH*ABS(FNB)/(DNZ+SMALL)))**5,ZERO)
              AP = (ALB-FLB)*FCLP + (ANB-FNB)*FCNP
              AB = ALB*FCL + ANB*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AB*BCX(1)
!
!-------  Outflow ---
!
          ELSEIF( IBCTX.EQ.7 ) THEN
            FLB = MIN( FLB,0.D+0 )
            FNB = MIN( FNB,0.D+0 )
!
!---------  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              WCLZ = 0.D+0
              IF( FLB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCLT = YL(NBT,NSL)/(SL(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBT,NSL)*FCLT)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                WCLZ = BCX(1)*FLB*THETA*FCL
     &            + C(N,NSL)*FLB*(1.D+0-THETA)*FCLP
              ELSEIF( FLB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCLZ = C(N,NSL)*FLB*FCLP
              ENDIF
              WCNZ = 0.D+0
              IF( FNB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCNT = YN(NBT,NSL)/(SN(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBT,NSL)*FCNT)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRNB,ISLC(1) )
                WCNZ = BCX(1)*FNB*THETA*FCN
     &            + C(N,NSL)*FNB*(1.D+0-THETA)*FCNP
              ELSEIF( FNB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCNZ = C(N,NSL)*FNB*FCNP
              ENDIF
              AB = 0.D+0
              AP = 0.D+0
              WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ+WCNZ)/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCNZ
!
!-----------  TVD Transport for interior surface 
!             adjacent to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVLX = SL(2,NBT)*PORD(2,NBT)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(SVLX+SMALL)
                FCLT = YL(NBT,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNT = AFZ(NQZ)*WN(1,NQZ)
              IF( FNT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVNX = SN(2,NBT)*PORD(2,NBT)
                CRNT = ABS( WN(1,NQZ) )*DT/DZGP(NQZ)/(SVNX+SMALL)
                FCNT = YN(NBT,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANB = MAX( FNB,ZERO )
              AP = (ALB-FLB)*FCLP + (ANB-FNB)*FCNP
              AB = ALB*FCL + ANB*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AB*BCX(1)
!
!-------- Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR.
     &      IBCTX.EQ.14 .OR. IBCTX.EQ.16 ) THEN
            FLB = MAX( FLB,0.D+0 )
            FNB = MAX( FNB,0.D+0 )
!
!---------  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              WCLZ = 0.D+0
              IF( FLB.GE.ZERO ) WCLZ = BCX(1)*FCL*FLB
              WCNZ = 0.D+0
              IF( FNB.GE.ZERO ) WCNZ = BCX(1)*FCN*FNB
              AB = 0.D+0
              AP = 0.D+0
              WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ+WCNZ)/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCNZ
!
!-----------  TVD Transport for interior surface 
!             adjacent to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVLX = SL(2,NBT)*PORD(2,NBT)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(SVLX+SMALL)
                FCLT = YL(NBT,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNT = AFZ(NQZ)*WN(1,NQZ)
              IF( FNT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVNX = SN(2,NBT)*PORD(2,NBT)
                CRNT = ABS( WN(1,NQZ) )*DT/DZGP(NQZ)/(SVNX+SMALL)
                FCNT = YN(NBT,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANB = MAX( FNB,ZERO )
              AP = (ALB-FLB)*FCLP + (ANB-FNB)*FCNP
              AB = ALB*FCL + ANB*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AB*BCX(1)
          ENDIF
!
!-----  South boundary  ---
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
          ELSE
            DPLS = 0.D+0
            DPNS = 0.D+0
          ENDIF
          FLS = AFY(NPY)*VL(1,NPY)
          FNS = AFY(NPY)*VN(1,NPY)
          CRLS = ABS( VL(1,NPY) )*DT/(RP(I)*DYGF(N)*XVLB+SMALL)
          CRNS = ABS( VN(1,NPY) )*DT/(RP(I)*DYGF(N)*XVNB+SMALL)
!
!---  Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.10 .OR.
     &      IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*XVLB*SMDLB
            INDX = 16
            DLY = DIFMN(DLB,DLP,DYGF(N),DYGF(N),VL(1,NPY),INDX)
            DLY = AFY(NPY)*(DLY+DPLS)/RP(I)/(5.D-1*DYGF(N))
            PCOR = (PNB(2,NB)+PATM)/PATM
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*XVNB*SMDNB
            INDX = 16
            DNY = DIFMN(DNB,DNP,DYGF(N),DYGF(N),VN(1,NPY),INDX)
            DNY = AFY(NPY)*(DNY+DPNS)/RP(I)/(5.D-1*DYGF(N))
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              IF( FLS.GE.ZERO ) THEN
                VCLY = BCX(1)*FCL*FLS
              ELSEIF( FLS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCLN = YL(NBN,NSL)/(SL(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBN,NSL)*FCLN)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                VCLY = BCX(1)*FLS*THETA*FCL
     &            + C(N,NSL)*FLS*(1.D+0-THETA)*FCLP
              ELSEIF( FLS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCLY = C(N,NSL)*FLS*FCLP
              ENDIF
              IF( FNS.GE.ZERO ) THEN
                VCNY = BCX(1)*FCN*FNS
              ELSEIF( FNS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCNN = YN(NBN,NSL)/(SN(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBN,NSL)*FCNN)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRNS,ISLC(1) )
                VCNY = BCX(1)*FNS*THETA*FCN
     &            + C(N,NSL)*FNS*(1.D+0-THETA)*FCNP
              ELSEIF( FNS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCNY = C(N,NSL)*FNS*FCNP
              ENDIF
              AS = DLY*FCL + DNY*FCN
              AP = DNY*FCNP + DLY*FCLP
              VC(NPY,NSL) = VC(NPY,NSL) + (VCLY+VCNY)/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCNY
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNN = AFY(NQY)*VN(1,NQY)
              IF( FNN.GE.ZERO ) THEN
                NBN = N+IFLD
                SVNX = SN(2,NBN)*PORD(2,NBN)
                CRNN = ABS( VN(1,NQY) )*DT/DYGP(NQY)/(SVNX*RP(I)+SMALL)
                FCNN = YN(NBN,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANS = MAX( FNS,ZERO ) +
     &          DNY*MAX((ONE-(TENTH*ABS(FNS)/(DNY+SMALL)))**5,ZERO)
              AP = (ALS-FLS)*FCLP + (ANS-FNS)*FCNP
              AS = ALS*FCL + ANS*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AS*BCX(1)
!
!---  Outflow ---
!
          ELSEIF( IBCTX.EQ.7 ) THEN
            FLS = MIN( FLS,0.D+0 )
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
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                VCLY = BCX(1)*FLS*THETA*FCL
     &            + C(N,NSL)*FLS*(1.D+0-THETA)*FCLP
              ELSEIF( FLS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCLY = C(N,NSL)*FLS*FCLP
              ENDIF
              VCNY = 0.D+0
              IF( FNS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCNN = YN(NBN,NSL)/(SN(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBN,NSL)*FCNN)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRNS,ISLC(1) )
                VCNY = BCX(1)*FNS*THETA*FCN
     &            + C(N,NSL)*FNS*(1.D+0-THETA)*FCNP
              ELSEIF( FNS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCNY = C(N,NSL)*FNS*FCNP
              ENDIF
              AS = 0.D+0
              AP = 0.D+0
              VC(NPY,NSL) = VC(NPY,NSL) + (VCLY+VCNY)/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCNY
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNN = AFY(NQY)*VN(1,NQY)
              IF( FNN.GE.ZERO ) THEN
                NBN = N+IFLD
                SVNX = SN(2,NBN)*PORD(2,NBN)
                CRNN = ABS( VN(1,NQY) )*DT/DYGP(NQY)/(SVNX*RP(I)+SMALL)
                FCNN = YN(NBN,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANS = MAX( FNS,ZERO )
              AP = (ALS-FLS)*FCLP + (ANS-FNS)*FCNP
              AS = ALS*FCL + ANS*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AS*BCX(1)
!
!--- Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR.
     &      IBCTX.EQ.14 .OR. IBCTX.EQ.16 ) THEN
            FLS = MAX( FLS,0.D+0 )
            FNS = MAX( FNS,0.D+0 )
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              VCLY = 0.D+0
              IF( FLS.GE.ZERO ) VCLY = BCX(1)*FCL*FLS
              VCNY = 0.D+0
              IF( FNS.GE.ZERO ) VCNY = BCX(1)*FCN*FNS
              AS = 0.D+0
              AP = 0.D+0
              VC(NPY,NSL) = VC(NPY,NSL) + (VCLY+VCNY)/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCNY
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNN = AFY(NQY)*VN(1,NQY)
              IF( FNN.GE.ZERO ) THEN
                NBN = N+IFLD
                SVNX = SN(2,NBN)*PORD(2,NBN)
                CRNN = ABS( VN(1,NQY) )*DT/DYGP(NQY)/(SVNX*RP(I)+SMALL)
                FCNN = YN(NBN,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANS = MAX( FNS,ZERO )
              AP = (ALS-FLS)*FCLP + (ANS-FNS)*FCNP
              AS = ALS*FCL + ANS*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AS*BCX(1)
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
          ELSE
            DPLW = 0.D+0
            DPNW = 0.D+0
          ENDIF
          FLW = AFX(NPX)*UL(1,NPX)
          FNW = AFX(NPX)*UN(1,NPX)
          CRLW = ABS( UL(1,NPX) )*DT/(DXGF(N)*XVLB+SMALL)
          CRNW = ABS( UN(1,NPX) )*DT/(DXGF(N)*XVNB+SMALL)
!
!---  Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.10 .OR.
     &      IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*XVLB*SMDLB
            INDX = 16
            DLX = DIFMN(DLB,DLP,DXGF(N),DXGF(N),UL(1,NPX),INDX)
            DLX = AFX(NPX)*(DLX+DPLW)/(5.D-1*DXGF(N))
            AL = MAX( UL(1,NPX),ZERO ) +
     &        DLX*MAX((ONE-(TENTH*ABS(UL(1,NPX))/(DLX+SMALL)))**5,ZERO)
            PCOR = (PNB(2,NB)+PATM)/PATM
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*XVNB*SMDNB
            INDX = 16
            DNX = DIFMN(DNB,DNP,DXGF(N),DXGF(N),UN(1,NPX),INDX)
            DNX = AFX(NPX)*(DNX+DPNW)/(5.D-1*DXGF(N))
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              IF( FLW.GE.ZERO ) THEN
                UCLX = BCX(1)*FCL*FLW
              ELSEIF( FLW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCLE = YL(NBE,NSL)/(SL(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBE,NSL)*FCLE)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA)*FCLP
     &            + BCX(1)*FLW*THETA*FCL
              ELSEIF( FLW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCLX = C(N,NSL)*FLW*FCLP
              ENDIF
              IF( FNW.GE.ZERO ) THEN
                UCNX = BCX(1)*FCN*FNW
              ELSEIF( FNW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCNE = YN(NBE,NSL)/(SN(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBE,NSL)*FCNE)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRNW,ISLC(1) )
                UCNX = C(N,NSL)*FNW*(1.D+0-THETA)*FCNP
     &            + BCX(1)*FNW*THETA*FCN
              ELSEIF( FNW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCNX = C(N,NSL)*FNW*FCNP
              ENDIF
              AW = DLX*FCL + DNX*FCN
              AP = DLX*FCLP + DNX*FCNP
              UC(NPX,NSL) = UC(NPX,NSL) + (UCLX+UCNX)/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCNX
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNE = AFX(NQX)*UN(1,NQX)
              IF( FNE.GE.ZERO ) THEN
                NBE = N+1
                SVNX = SN(2,NBE)*PORD(2,NBE)
                CRNE = ABS( UN(1,NQX) )*DT/DXGP(NQX)/(SVNX+SMALL)
                FCNE = YN(NBE,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANW = MAX(FNW,ZERO)
     &          + DNX*MAX((ONE-(TENTH*ABS(FNW)/(DNX+SMALL)))**5,ZERO)
              AP = (ALW-FLW)*FCLP + (ANW-FNW)*FCNP
              AW = ALW*FCL + ANW*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AW*BCX(1)
!
!---  Outflow ---
!
          ELSEIF( IBCTX.EQ.7 ) THEN
            FLW = MIN( FLW,0.D+0 )
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
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA)*FCLP
     &            + BCX(1)*FLW*THETA*FCL
              ELSEIF( FLW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCLX = C(N,NSL)*FLW*FCLP
              ENDIF
              UCNX = 0.D+0
              IF( FNW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCNE = YN(NBE,NSL)/(SN(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBE,NSL)*FCNE)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRNW,ISLC(1) )
                UCNX = C(N,NSL)*FNW*(1.D+0-THETA)*FCNP
     &            + BCX(1)*FNW*THETA*FCN
              ELSEIF( FNW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCNX = C(N,NSL)*FNW*FCNP
              ENDIF
              AW = 0.D+0
              AP = 0.D+0
              UC(NPX,NSL) = UC(NPX,NSL) + (UCLX+UCNX)/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCNX
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNE = AFX(NQX)*UN(1,NQX)
              IF( FNE.GE.ZERO ) THEN
                NBE = N+1
                SVNX = SN(2,NBE)*PORD(2,NBE)
                CRNE = ABS( UN(1,NQX) )*DT/DXGP(NQX)/(SVNX+SMALL)
                FCNE = YN(NBE,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANW = MAX(FNW,ZERO)
              AP = (ALW-FLW)*FCLP + (ANW-FNW)*FCNP
              AW = ALW*FCL + ANW*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AW*BCX(1)
!
!--- Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR.
     &      IBCTX.EQ.14 .OR. IBCTX.EQ.16 ) THEN
            FLW = MAX( FLW,0.D+0 )
            FNW = MAX( FNW,0.D+0 )
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              UCLX = 0.D+0
              IF( FLW.GE.ZERO ) UCLX = BCX(1)*FCL*FLW
              UCNX = 0.D+0
              IF( FNW.GE.ZERO ) UCNX = BCX(1)*FCN*FNW
              AW = 0.D+0
              AP = 0.D+0
              UC(NPX,NSL) = UC(NPX,NSL) + (UCLX+UCNX)/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCNX
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNE = AFX(NQX)*UN(1,NQX)
              IF( FNE.GE.ZERO ) THEN
                NBE = N+1
                SVNX = SN(2,NBE)*PORD(2,NBE)
                CRNE = ABS( UN(1,NQX) )*DT/DXGP(NQX)/(SVNX+SMALL)
                FCNE = YN(NBE,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANW = MAX(FNW,ZERO)
              AP = (ALW-FLW)*FCLP + (ANW-FNW)*FCNP
              AW = ALW*FCL + ANW*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AW*BCX(1)
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
          ELSE
            DPLE = 0.D+0
            DPNE = 0.D+0
          ENDIF
          FLE = AFX(NQX)*UL(1,NQX)
          FNE = AFX(NQX)*UN(1,NQX)
          CRLE = ABS( UL(1,NQX) )*DT/(DXGF(N)*XVLB+SMALL)
          CRNE = ABS( UN(1,NQX) )*DT/(DXGF(N)*XVNB+SMALL)
!
!---  Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.10 .OR.
     &      IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*XVLB*SMDLB
            INDX = 16
            DLX = DIFMN(DLP,DLB,DXGF(N),DXGF(N),UL(1,NQX),INDX)
            DLX = AFX(NQX)*(DLX+DPLE)/(5.D-1*DXGF(N))
            AL = MAX( -UL(1,NQX),ZERO ) +
     &        DLX*MAX((ONE-(TENTH*ABS(UL(1,NQX))/(DLX+SMALL)))**5,ZERO)
            PCOR = (PNB(2,NB)+PATM)/PATM
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*XVNB*SMDNB
            INDX = 16
            DNX = DIFMN(DNP,DNB,DXGF(N),DXGF(N),UN(1,NQX),INDX)
            DNX = AFX(NQX)*(DNX+DPNE)/(5.D-1*DXGF(N))
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              IF( FLE.LT.ZERO ) THEN
                UCLX = BCX(1)*FCL*FLE
              ELSEIF( FLE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCLW = YL(NBW,NSL)/(SL(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBW,NSL)*FCLW)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                UCLX =  C(N,NSL)*FLE*(1.D+0-THETA)*FCLP
     &          +  BCX(1)*FLE*THETA*FCL
              ELSEIF( FLE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCLX =  C(N,NSL)*FLE*FCLP
              ENDIF
              IF( FNE.LT.ZERO ) THEN
                UCNX = BCX(1)*FCN*FNE
              ELSEIF( FNE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCNW = YN(NBW,NSL)/(SN(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBW,NSL)*FCNW)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRNE,ISLC(1) )
                UCNX =  C(N,NSL)*FNE*(1.D+0-THETA)*FCNP
     &          +  BCX(1)*FNE*THETA*FCN
              ELSEIF( FNE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCNX =  C(N,NSL)*FNE*FCNP
              ENDIF
              AE = DLX*FCL + DNX*FCN
              AP = DLX*FCLP + DNX*FCNP
              UC(NQX,NSL) = UC(NQX,NSL) + (UCLX+UCNX)/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCNX
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNW = AFX(NPX)*UN(1,NPX)
              IF( FNW.LT.ZERO ) THEN
                NBW = N-1
                SVNX = SN(2,NBW)*PORD(2,NBW)
                CRNW = ABS( UN(1,NPX) )*DT/DXGP(NPX)/(SVNX+SMALL)
                FCNW = YN(NBW,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANE = MAX( -FNE,ZERO ) +
     &          DNX*MAX((ONE-(TENTH*ABS(FNE)/(DNX+SMALL)))**5,ZERO)
              AP = (ALE+FLE)*FCLP + (ANE+FNE)*FCNP
              AE = ALE*FCL + ANE*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AE*BCX(1)
!
!---  Outflow ---
!
          ELSEIF( IBCTX.EQ.7 ) THEN
            FLE = MAX( FLE,0.D+0 )
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
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                UCLX =  C(N,NSL)*FLE*(1.D+0-THETA)*FCLP
     &          +  BCX(1)*FLE*THETA*FCL
              ELSEIF( FLE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCLX =  C(N,NSL)*FLE*FCLP
              ENDIF
              UCNX = 0.D+0
              IF( FNE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCNW = YN(NBW,NSL)/(SN(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBW,NSL)*FCNW)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRNE,ISLC(1) )
                UCNX =  C(N,NSL)*FNE*(1.D+0-THETA)*FCNP
     &          +  BCX(1)*FNE*THETA*FCN
              ELSEIF( FNE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCNX =  C(N,NSL)*FNE*FCNP
              ENDIF
              AE = 0.D+0
              AP = 0.D+0
              UC(NQX,NSL) = UC(NQX,NSL) + (UCLX+UCNX)/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCNX
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNW = AFX(NPX)*UN(1,NPX)
              IF( FNW.LT.ZERO ) THEN
                NBW = N-1
                SVNX = SN(2,NBW)*PORD(2,NBW)
                CRNW = ABS( UN(1,NPX) )*DT/DXGP(NPX)/(SVNX+SMALL)
                FCNW = YN(NBW,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANE = MAX( -FNE,ZERO )
              AP = (ALE+FLE)*FCLP + (ANE+FNE)*FCNP
              AE = ALE*FCL + ANE*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AE*BCX(1)
!
!--- Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR.
     &      IBCTX.EQ.14 .OR. IBCTX.EQ.16 ) THEN
            FLE = MIN( FLE,0.D+0 )
            FNE = MIN( FNE,0.D+0 )
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              UCLX = 0.D+0
              IF( FLE.LT.ZERO ) UCLX = BCX(1)*FCL*FLE
              UCNX = 0.D+0
              IF( FNE.LT.ZERO ) UCNX = BCX(1)*FCN*FNE
              AE = 0.D+0
              AP = 0.D+0
              UC(NQX,NSL) = UC(NQX,NSL) + (UCLX+UCNX)/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCNX
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNW = AFX(NPX)*UN(1,NPX)
              IF( FNW.LT.ZERO ) THEN
                NBW = N-1
                SVNX = SN(2,NBW)*PORD(2,NBW)
                CRNW = ABS( UN(1,NPX) )*DT/DXGP(NPX)/(SVNX+SMALL)
                FCNW = YN(NBW,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANE = MAX( -FNE,ZERO )
              AP = (ALE+FLE)*FCLP + (ANE+FNE)*FCNP
              AE = ALE*FCL + ANE*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AE*BCX(1)
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
          ELSE
            DPLN = 0.D+0
            DPNN = 0.D+0
          ENDIF
          FLN = AFY(NQY)*VL(1,NQY)
          FNN = AFY(NQY)*VN(1,NQY)
          CRLN = ABS( VL(1,NQY) )*DT/(RP(I)*DYGF(N)*XVLB+SMALL)
          CRNN = ABS( VN(1,NQY) )*DT/(RP(I)*DYGF(N)*XVNB+SMALL)
!
!---  Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.10 .OR.
     &      IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*XVLB*SMDLB
            INDX = 16
            DLY = DIFMN(DLP,DLB,DYGF(N),DYGF(N),VL(1,NQY),INDX)
            DLY = AFY(NQY)*(DLY+DPLN)/RP(I)/(5.D-1*DYGF(N))
            AL = MAX( -VL(1,NQY),ZERO ) +
     &        DLY*MAX((ONE-(TENTH*ABS(VL(1,NQY))/(DLY+SMALL)))**5,ZERO)
            PCOR = (PNB(2,NB)+PATM)/PATM
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*XVNB*SMDNB
            INDX = 16
            DNY = DIFMN(DNP,DNB,DYGF(N),DYGF(N),VN(1,NQY),INDX)
            DNY = AFY(NQY)*(DNY+DPNN)/RP(I)/(5.D-1*DYGF(N))
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              IF( FLN.LT.ZERO ) THEN
                VCLY = BCX(1)*FCL*FLN
              ELSEIF( FLN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCLS = YL(NBS,NSL)/(SL(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBS,NSL)*FCLS)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                VCLY =  BCX(1)*FLN*THETA*FCL
     &            + C(N,NSL)*FLN*(1.D+0-THETA)*FCLP
              ELSEIF( FLN.GE.ZERO .AND. J.GT.1 ) THEN
                VCLY =  C(N,NSL)*FLN*FCLP
              ENDIF
              IF( FNN.LT.ZERO ) THEN
                VCNY = BCX(1)*FCN*FNN
              ELSEIF( FNN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCNS = YN(NBS,NSL)/(SN(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBS,NSL)*FCNS)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRNN,ISLC(1) )
                VCNY =  BCX(1)*FNN*THETA*FCN
     &            + C(N,NSL)*FNN*(1.D+0-THETA)*FCNP
              ELSEIF( FNN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCNY =  C(N,NSL)*FNN*FCNP
              ENDIF
              AN = DLY*FCL + DNY*FCN
              AP = DLY*FCLP + DNY*FCNP
              VC(NQY,NSL) = VC(NQY,NSL) + (VCLY+VCNY)/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCNY
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNS = AFY(NPY)*VN(1,NPY)
              IF( FNS.LT.ZERO ) THEN
                NBS = N-IFLD
                SVNX = SN(2,NBS)*PORD(2,NBS)
                CRNS = ABS( VN(1,NPY) )*DT/DYGP(NPY)/(SVNX*RP(I)+SMALL)
                FCNS = YN(NBS,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANN = MAX( -FNN,ZERO ) +
     &          DNY*MAX((ONE-(TENTH*ABS(FNN)/(DNY+SMALL)))**5,ZERO)
              AP = (ALN+FLN)*FCLP + (ANN+FNN)*FCNP
              AN = ALN*FCL + ANN*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AN*BCX(1)
!
!---  Outflow ---
!
          ELSEIF( IBCTX.EQ.7 ) THEN
            FLN = MAX( FLN,0.D+0 )
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
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                VCLY =  BCX(1)*FLN*THETA*FCL
     &            + C(N,NSL)*FLN*(1.D+0-THETA)*FCLP
              ELSEIF( FLN.GE.ZERO .AND. J.GT.1 ) THEN
                VCLY =  C(N,NSL)*FLN*FCLP
              ENDIF
              VCNY = 0.D+0
              IF( FNN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCNS = YN(NBS,NSL)/(SN(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBS,NSL)*FCNS)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRNN,ISLC(1) )
                VCNY =  BCX(1)*FNN*THETA*FCN
     &            + C(N,NSL)*FNN*(1.D+0-THETA)*FCNP
              ELSEIF( FNN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCNY =  C(N,NSL)*FNN*FCNP
              ENDIF
              AN = 0.D+0
              AP = 0.D+0
              VC(NQY,NSL) = VC(NQY,NSL) + (VCLY+VCNY)/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCNY
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNS = AFY(NPY)*VN(1,NPY)
              IF( FNS.LT.ZERO ) THEN
                NBS = N-IFLD
                SVNX = SN(2,NBS)*PORD(2,NBS)
                CRNS = ABS( VN(1,NPY) )*DT/DYGP(NPY)/(SVNX*RP(I)+SMALL)
                FCNS = YN(NBS,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANN = MAX( -FNN,ZERO )
              AP = (ALN+FLN)*FCLP + (ANN+FNN)*FCNP
              AN = ALN*FCL + ANN*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AN*BCX(1)
!
!--- Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR.
     &      IBCTX.EQ.14 .OR. IBCTX.EQ.16 ) THEN
            FLN = MIN( FLN,0.D+0 )
            FNN = MIN( FNN,0.D+0 )
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              VCLY = 0.D+0
              IF( FLN.LT.ZERO ) VCLY = BCX(1)*FCL*FLN
              VCNY = 0.D+0
              IF( FNN.LT.ZERO ) VCNY = BCX(1)*FCN*FNN
              AN = 0.D+0
              AP = 0.D+0
              VC(NQY,NSL) = VC(NQY,NSL) + (VCLY+VCNY)/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCNY
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNS = AFY(NPY)*VN(1,NPY)
              IF( FNS.LT.ZERO ) THEN
                NBS = N-IFLD
                SVNX = SN(2,NBS)*PORD(2,NBS)
                CRNS = ABS( VN(1,NPY) )*DT/DYGP(NPY)/(SVNX*RP(I)+SMALL)
                FCNS = YN(NBS,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANN = MAX( -FNN,ZERO )
              AP = (ALN+FLN)*FCLP + (ANN+FNN)*FCNP
              AN = ALN*FCL + ANN*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AN*BCX(1)
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
          ELSE
            DPLT = 0.D+0
            DPNT = 0.D+0
          ENDIF
          FLT = AFZ(NQZ)*WL(1,NQZ)
          FNT = AFZ(NQZ)*WN(1,NQZ)
          CRLT = ABS( WL(1,NQZ) )*DT/(DZGF(N)*XVLB+SMALL)
          CRNT = ABS( WN(1,NQZ) )*DT/(DZGF(N)*XVNB+SMALL)
!
!---  Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.10 .OR.
     &      IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*XVLB*SMDLB
            INDX = 16
            DLZ = DIFMN(DLP,DLB,DZGF(N),DZGF(N),WL(1,NQZ),INDX)
            DLZ = AFZ(NQZ)*(DLZ+DPLT)/(5.D-1*DZGF(N))
            AL = MAX( -WL(1,NQZ),ZERO ) +
     &        DLZ*MAX((ONE-(TENTH*ABS(WL(1,NQZ))/(DLZ+SMALL)))**5,ZERO)
            PCOR = (PNB(2,NB)+PATM)/PATM
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*XVNB*SMDNB
            INDX = 16
            DNZ = DIFMN(DNP,DNB,DZGF(N),DZGF(N),WN(1,NQZ),INDX)
            DNZ = AFZ(NQZ)*(DNZ+DPNT)/(5.D-1*DZGF(N))
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              IF( FLT.LT.ZERO ) THEN
                WCLZ = BCX(1)*FCL*FLT
              ELSEIF( FLT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCLB = YL(NBB,NSL)/(SL(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBB,NSL)*FCLB)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                WCLZ =  C(N,NSL)*FLT*(1.D+0-THETA)*FCLP
     &            + BCX(1)*FLT*THETA*FCL
              ELSEIF( FLT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCLZ =  C(N,NSL)*FLT*FCLP
              ENDIF
              IF( FNT.LT.ZERO ) THEN
                WCNZ = BCX(1)*FCN*FNT
              ELSEIF( FNT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCNB = YN(NBB,NSL)/(SN(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBB,NSL)*FCNB)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRNT,ISLC(1) )
                WCNZ =  C(N,NSL)*FNT*(1.D+0-THETA)*FCNP
     &            + BCX(1)*FNT*THETA*FCN
              ELSEIF( FNT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCNZ =  C(N,NSL)*FNT*FCNP
              ENDIF
              AT = DLZ*FCL + DNZ*FCN
              AP = DLZ*FCLP + DNZ*FCNP
              WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ+WCNZ)/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCNZ
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNB = AFZ(NPZ)*WN(1,NPZ)
              IF( FNB.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVNX = SN(2,NBB)*PORD(2,NBB)
                CRNB = ABS( WN(1,NPZ) )*DT/DZGP(NPZ)/(SVNX+SMALL)
                FCNB = YN(NBB,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANT = MAX( -FNT,ZERO ) +
     &          DNZ*MAX((ONE-(TENTH*ABS(FNT)/(DNZ+SMALL)))**5,ZERO)
              AP = (ALT+FLT)*FCLP + (ANT+FNT)*FCNP
              AT = ALT*FCL + ANT*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AT*BCX(1)
!
!---  Outflow ---
!
          ELSEIF( IBCTX.EQ.7 ) THEN
            FLT = MAX( FLT,0.D+0 )
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
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                WCLZ =  C(N,NSL)*FLT*(1.D+0-THETA)*FCLP
     &            + BCX(1)*FLT*THETA*FCL
              ELSEIF( FLT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCLZ =  C(N,NSL)*FLT*FCLP
              ENDIF
              WCNZ = 0.D+0
              IF( FNT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCNB = YN(NBB,NSL)/(SN(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBB,NSL)*FCNB)
     &            /(BCX(1)*FCN-C(N,NSL)*FCNP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRNT,ISLC(1) )
                WCNZ =  C(N,NSL)*FNT*(1.D+0-THETA)*FCNP
     &            + BCX(1)*FNT*THETA*FCN
              ELSEIF( FNT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCNZ =  C(N,NSL)*FNT*FCNP
              ENDIF
              AT = 0.D+0
              AP = 0.D+0
              WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ+WCNZ)/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCNZ
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNB = AFZ(NPZ)*WN(1,NPZ)
              IF( FNB.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVNX = SN(2,NBB)*PORD(2,NBB)
                CRNB = ABS( WN(1,NPZ) )*DT/DZGP(NPZ)/(SVNX+SMALL)
                FCNB = YN(NBB,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANT = MAX( -FNT,ZERO )
              AP = (ALT+FLT)*FCLP + (ANT+FNT)*FCNP
              AT = ALT*FCL + ANT*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AT*BCX(1)
!
!--- Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR.
     &      IBCTX.EQ.14 .OR. IBCTX.EQ.16 ) THEN
            FLT = MIN( FLT,0.D+0 )
            FNT = MIN( FNT,0.D+0 )
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              WCLZ = 0.D+0
              IF( FLT.LT.ZERO ) WCLZ = BCX(1)*FCL*FLT
              WCNZ = 0.D+0
              IF( FNT.LT.ZERO ) WCNZ = BCX(1)*FCN*FNT
              AT = 0.D+0
              AP = 0.D+0
              WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ+WCNZ)/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCNZ
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
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
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
              FNB = AFZ(NPZ)*WN(1,NPZ)
              IF( FNB.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVNX = SN(2,NBB)*PORD(2,NBB)
                CRNB = ABS( WN(1,NPZ) )*DT/DZGP(NPZ)/(SVNX+SMALL)
                FCNB = YN(NBB,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX(1)*FCN)
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
              ANT = MAX( -FNT,ZERO )
              AP = (ALT+FLT)*FCLP + (ANT+FNT)*FCNP
              AT = ALT*FCL + ANT*FCN
            ENDIF
            BLU(MP) = BLU(MP) + AT*BCX(1)
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
!---  End of SBND4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFT4( XLOX,IZN )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Water-Oil Mode
!
!     Gas-aqueous interfacial tension and contact angle as a function
!     of dissolved-oil concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 April 2005.
!     Last Modified by MD White, PNNL, 6 April 2005.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/SFT4'
!
!---  Fluid-pair interfacial tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
!
!---  Gas-aqueous interfacial tension as a function of dissolved-oil
!     concentration  ---
!
      IF( ABS(SFCSF(1)).GT.EPSL ) THEN
        RHOLOX = XLOX*RHORL
        SIG_GL = SFCSF(1) + SFCSF(2)*EXP( -SFCSF(3)*RHOLOX )
        IF( ABS(SFCSF(4)).GT.EPSL ) THEN
          SIG_GL = SIG_GL + SFCSF(4)*EXP( -SFCSF(5)*RHOLOX )
        ENDIF
!
!---    NAPL-aqueous interfacial tension as a function of gas-NAPL
!       contact angle, gas-NAPL surface tension, and gas-aqueous
!       surface tension  ---
!
!        SIG_NL = SQRT( -2.D+0*CA_GN*SIG_GL*SIG_GN + 
!     &    (SIG_GL**2) + (SIG_GN**2) )
!
!---    NAPL-aqueous contact angle as a function of NAPL-aqueous,
!       gas-NAPL, and gas-aqueous interfacial tensions  ---
!
        CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &    /(2.D+0*SIG_GL*SIG_NL)
        CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &    /(2.D+0*SIG_GL*SIG_GN)
      ENDIF
!
!---  Fluid-pair scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFT4 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
#ifdef awtap
      SUBROUTINE SORC4( HCFL,HCFN )
#else
      SUBROUTINE SORC4
#endif
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Water-Oil Mode
!
!     Compute source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE REACT
      USE PORMED
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FDVP
      USE FDVN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
#ifdef awtap
      REAL*8 HCFL(LFZ,LSR),HCFN(LFZ,LSR)
#endif
      REAL*8 SRX(8)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORC4'
#ifdef ecke
!
!---  Load oil sources for linked aqueous oil   ---
!
      IF( ISPLK(6).NE.0 ) THEN
!
!---    Load oil sources associated with chemical
!       reactions  ---
!
        DO 20 N = 1,NFLD
          SRCOX = SRCO(1,N)*DTI
          DO 10 M = 2,ISVC+2
            SRCO(M,N) = SRCOX
   10     CONTINUE
   20   CONTINUE
!
!---    Zero source terms  ---
!
        IF( NSR.GT.0 ) THEN
          DO 40 N = 1,NFLD
            DO 30 M = 2,ISVC+2
              SRCW(M,N) = 0.D+0
   30       CONTINUE
   40     CONTINUE
        ENDIF
      ELSE
#endif
!
!---    Zero source terms  ---
!
        IF( NSR.GT.0 ) THEN
          DO 60 N = 1,NFLD
            DO 50 M = 2,ISVC+2
              SRCW(M,N) = 0.D+0
              SRCO(M,N) = 0.D+0
   50       CONTINUE
   60     CONTINUE
        ENDIF
#ifdef ecke
      ENDIF
#endif
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
#ifdef awtap
        IF( ISRT(NS).LT.20 .OR. ISRT(NS).GT.29 ) THEN
#endif
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
!
!---            Aqueous Volumetric Rate w/ Component Mass Fractions ---
!
                IF( ISRT(NS).EQ.3 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    SRCW(M,N) = SRCW(M,N) +
     &                SRX(4)*RHOL(M,N)*(1.D+0-SRX(6))
                    SRCO(M,N) = SRCO(M,N) + SRX(4)*RHOL(M,N)*SRX(6)
                  ELSE
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOL(M,N)*XLW(M,N)
                    SRCO(M,N) = SRCO(M,N) + SRX(4)*RHOL(M,N)*XLO(M,N)
                  ENDIF
!
!---            Total Liquid Volumetric Rate 
!               w/ Component Volume Fractions ---
!
                ELSEIF( ISRT(NS).EQ.4 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    SRCW(M,N) = SRCW(M,N) +
     &                SRX(4)*RHOL(M,N)*(1.D+0-SRX(6))
                    SRCO(M,N) = SRCO(M,N) + SRX(4)*RHON(M,N)*SRX(6)
                  ELSE
                    RKLX = (RKL(1,M,N)*RKL(2,M,N)*RKL(3,M,N))**THIRD
                    RLX = RKLX/VISL(M,N)
                    RNX = RKN(M,N)/VISN(M,N)
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOL(M,N)*XLW(M,N)*
     &                (RLX/(RLX+RNX))
                    SRCO(M,N) = SRCO(M,N) + SRX(4)*(RHOL(M,N)*XLO(M,N)*
     &                (RLX/(RLX+RNX)) + RHON(M,N)*(RNX/(RLX+RNX)))
                  ENDIF
!
!---            NAPL Volumetric Rate  ---
!
                ELSEIF( ISRT(NS).EQ.6 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    SRCO(M,N) = SRCO(M,N) + SRX(4)*RHON(M,N)
                  ELSE
                    SRCO(M,N) = SRCO(M,N) + SRX(4)*RHON(M,N)
                  ENDIF
!
!---            Aqueous Mass Rate w/ Component Mass Fractions ---
!
                ELSEIF( ISRT(NS).EQ.7 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*(1.D+0-SRX(6))
                    SRCO(M,N) = SRCO(M,N) + SRX(4)*SRX(6)
                  ELSE
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*XLW(M,N)
                    SRCO(M,N) = SRCO(M,N) + SRX(4)*XLO(M,N)
                  ENDIF
!
!---            Total Liquid Mass Rate 
!               w/ Component Mass Fractions ---
!
                ELSEIF( ISRT(NS).EQ.8 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*(1.D+0-SRX(6))
                    SRCO(M,N) = SRCO(M,N) + SRX(4)*SRX(6)
                  ELSE
                    RKLX = (RKL(1,M,N)*RKL(2,M,N)*RKL(3,M,N))**THIRD
                    RLX = RKLX/VISL(M,N)
                    RNX = RKN(M,N)/VISN(M,N)
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*XLW(M,N)*
     &                (RLX/(RLX+RNX))
                    SRCO(M,N) = SRCO(M,N) + SRX(4)*(XLO(M,N)*
     &                (RLX/(RLX+RNX)) + (RNX/(RLX+RNX)))
                  ENDIF
!
!---            NAPL Mass Rate ---
!
                ELSEIF( ISRT(NS).EQ.10 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    SRCO(M,N) = SRCO(M,N) + SRX(4)
                  ELSE
                    SRCO(M,N) = SRCO(M,N) + SRX(4)
                  ENDIF
                ENDIF
  400         CONTINUE
  500   CONTINUE
#ifdef awtap
!
!---  Constant liquid level well  ---
!
        ELSEIF( ISRT(NS).EQ.21 ) THEN
          I = ISRDM(1,NS)
          J = ISRDM(3,NS)
          K1X = ISRDM(5,NS)
          K2X = ISRDM(6,NS)
          RDW = SRX(3)
          N = ND(I,J,K1X)
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. (I.EQ.1) ) THEN
            RDE = XGR(I+1,J,K1X)
          ELSE
            RDE = SQRT( AFZ(NSZ(N))/GPI )
          ENDIF
          NC = 0
!
!---  Compute the well bottom pressure on the first iteration  ---
!
!         IF( NITER.GT.1 ) GOTO 509
!
!---  Determine well gas pressure  ---
!
          PGW(2,NS) = 0.D+0
          IFLAG = 0
          DO 504 K = K1X,KFLD
            N = ND(I,J,K)
            IF( NPHAZ(1,N).EQ.2 .OR. NPHAZ(1,N).EQ.4 .OR.
     &        NPHAZ(1,N).EQ.6 ) THEN
              PGW(2,NS) = PGW(2,NS) + PG(2,N)*DZGF(N)
              DZW = DZW + DZGF(N)
              IFLAG = 1
            ENDIF
  504     CONTINUE
          IF( IFLAG.EQ.1 ) PGW(2,NS) = PGW(2,NS)/DZW
!
!---  Loop over the screened interval up to the fluid height
!     in the well to determine the average aqueous-NAPL
!     relative mobility   ---
!
          RMBX = 0.D+0
          NC = 0
          HZ = 0.D+0
          DO 506 K = K1X,K2X
            N = ND(I,J,K)
            IZN = IZ(N)
            HZ = HZ + 5.D-1*DZGF(N)
            IF( HZ.GT.SRX(4) ) GOTO 506
            HZ = HZ + 5.D-1*DZGF(N)
            NC = NC+1
            RKLMX = (RKL(1,2,N)*RKL(2,2,N)*RKL(3,2,N))**(1./3.)
            RMBX = RMBX + RKLMX*VISN(2,N)/(RKN(2,N)*VISL(2,N)+EPSL)
  506     CONTINUE
          REALX = REAL(NC)
          RMBX = RMBX/REALX
!
!---  Determine the pressure at the bottom of the well based on
!     the density of the fluid in the well; where the fluid density
!     is dependent on average relative mobility of the
!     fluids in the well  ---
!
          NC = 128
          REALX = REAL(NC)
          DZW = SRX(4)/REALX
          PLWB(2,NS) = PGW(2,NS)
          DO 508 K = 1,NC
            PWX = PLWB(2,NS)+PATM
            CALL WATLQD( T(2,N),PWX,RHOWL )
            CALL VOCLQD( T(2,N),PWX,PVO(2,N),RHOWN )
            RHOWT = (RMBX*RHOWL + RHOWN)/(1.D+0 + RMBX)
            PLWB(2,NS) = PLWB(2,NS) + DZW*GRAVZ*RHOWT
  508     CONTINUE
!
!---  Loop over the screened interval of the well to determine
!     the pumping rate for the aqueous and NAPL  ---
!
!          PLWB(2,NS) = SRX(4)*RHORL*GRAVZ
!          PGW(2,NS) = 0.D+0
  509     CONTINUE
          DO 580 M = 2,ISVC+2
          PLW = PLWB(2,NS)
          QWSL = 0.D+0
          QWSN = 0.D+0
          DRD2 = (RDE**2-RDW**2)
          DO 510 K = K1X,K2X
            N = ND(I,J,K)
            IZN = IZ(N)
!
!---  Compute the average radial permeability  ---
!
            IF( PERM(1,IZN)/EPSL.GT.EPSL )  THEN
              IF( PERM(2,IZN)/EPSL.GT.EPSL ) THEN
                PERMX = SQRT( PERM(1,IZN)*PERM(2,IZN) )
              ELSE
                PERMX = PERM(1,IZN)
              ENDIF
            ELSE
              PERMX = PERM(2,IZN)
            ENDIF
!
!---  Compute the local well fluid density and pressure  ---
!
            PWX = MAX( PLW+PATM,PGW(2,NS)+PATM )
            CALL WATLQD( T(2,N),PWX,RHOWL )
            CALL VOCLQD( T(2,N),PWX,PVO(2,N),RHOWN )
            RHOWT = (RMBX*RHOWL + RHOWN)/(1.D+0 + RMBX)
            BUOY = RHOWT*GRAVZ*DZGF(N)
            PLWP = MAX( PLW-5.D-1*BUOY,PGW(2,NS) )
!
!---  Consider well bore seepage  ---
!
            AFL = MIN( (PL(M,N)-PGW(2,NS))/BUOY,1.D+0 )
            AFN = MIN( (PN(M,N)-PGW(2,NS))/BUOY,1.D+0 )
!
!---  Compute the aqueous volumetric flux  ---
!
            IF( AFL/EPSL.GT.EPSL ) THEN
              RKLMX = (RKL(1,M,N)*RKL(2,M,N)*RKL(3,M,N))**(1./3.)
              HCFL(K,NS) = 2.D+0*GPI*PERMX*DRD2*DZGF(N)*RKLMX/
     &          (VISL(M,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
     &          *(PLWP-PL(M,N))*AFL
            ELSE
              HCFL(K,NS) = 0.D+0
            ENDIF
!
!---  Compute the NAPL volumetric flux  ---
!
            IF( AFN/EPSL.GT.EPSL ) THEN
              HCFN(K,NS) = 2.D+0*GPI*PERMX*DRD2*DZGF(N)*RKN(M,N)/
     &          (VISN(M,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
     &          *(PLWP-PN(M,N))*AFN
            ELSE
              HCFN(K,NS) = 0.D+0
            ENDIF
!
!---  Sum the volumetric fluxes  ---
!
            IF( M.EQ.2 ) THEN
              QWSL = QWSL + HCFL(K,NS)
              QWSN = QWSN + HCFN(K,NS)
            ENDIF
!
!---  Compute the well pressure at the next elevation  ---
!
            PWX = PLWP+PATM
            CALL WATLQD( T(2,N),PWX,RHOWL )
            CALL VOCLQD( T(2,N),PWX,PVO(2,N),RHOWN )
            RHOWT = (RMBX*RHOWL + RHOWN)/(1.D+0 + RMBX)
            BUOY = RHOWT*GRAVZ*DZGF(N)
            PLW = PLW-BUOY
  510     CONTINUE
!
!---  Assign well-bottom pressure across indices  ---
!
          PLWB(M,NS) = PLWB(2,NS)
!
!---  Simplified well model  ---
!
          PLW = PLWB(2,NS)
!
!---  Maximum pumping rate  ---
!
          IF( M.EQ.2 ) THEN
            QLW(3,NS) = QWSL
            QNW(3,NS) = QWSN
            QTW(3,NS) = (QWSL + QWSN)
          ENDIF
          DO 516 K = K1X,K2X
            N = ND(I,J,K)
            IZN = IZ(N)
            PWX = MAX( PLW+PATM,PGW(2,NS)+PATM )
            CALL WATLQD( T(2,N),PWX,RHOWL )
            CALL VOCLQD( T(2,N),PWX,PVO(2,N),RHOWN )
            RHOWT = (RMBX*RHOWL + RHOWN)/(1.D+0 + RMBX)
            BUOY = RHOWT*GRAVZ*DZGF(N)
            PLWP = MAX( PLW-5.D-1*BUOY,PGW(2,NS) )
!
!--- For flow into the formation use computed aqueous concentrations,
!    no NAPL injection is considered
!
            IF( HCFL(K,NS)/EPSL.GT.EPSL ) THEN
              PWX = PLWP+PATM
              CALL WATLQD( T(2,N),PWX,RHOLX )
              CALL VOCSP( T(2,N),PSOX )
              XMLOX = SRX(6)*PSOX/HCOW
              XMLWX = MAX( 1.D+0-XMLOX,0.D+0 )
              XLOX = XMLOX*WTMO/(XMLOX*WTMO+XMLWX*WTMW)
              XLWX = MAX( 1.D+0-XLOX,0.D+0 )
              XLOX = XLOX*RHOLX
              XLWX = XLWX*RHOLX
!
!--- For flow from the formation use local aqueous concentrations
!
            ELSE
              XLOX = XLO(M,N)*RHOL(M,N)
              XLWX = XLW(M,N)*RHOL(M,N)
              XNOX = RHON(M,N)
            ENDIF
            SRCW(M,N) = SRCW(M,N) + HCFL(K,NS)*XLWX
            SRCO(M,N) = SRCO(M,N) + (HCFL(K,NS)*XLOX + HCFN(K,NS)*XNOX)
            PWX = PLWP+PATM
            CALL WATLQD( T(2,N),PWX,RHOWL )
            CALL VOCLQD( T(2,N),PWX,PVO(2,N),RHOWN )
            RHOWT = (RMBX*RHOWL + RHOWN)/(1.D+0 + RMBX)
            BUOY = RHOWT*GRAVZ*DZGF(N)
            PLW = PLW-BUOY
  516     CONTINUE
  580     CONTINUE
  590     CONTINUE
        ENDIF
#endif
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORC4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
#ifdef awtap
      SUBROUTINE SORIT4( HCFL,HCFN,NSL )
#else
      SUBROUTINE SORIT4( NSL )
#endif
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Water-Oil Mode
!
!     Compute solute transport source integrals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, June 1994.
!     Last Modified by MD White, Battelle, June 9, 1994.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 SRX(8)
#ifdef awtap
      REAL*8 HCFL(LFZ,LSR),HCFN(LFZ,LSR)
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORIT4'
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
#ifdef future_coding
        IF( ISRT(NS).LT.20 .OR. ISRT(NS).GT.29 ) THEN
#endif
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
#ifdef future_coding
!
!---  Well source  ---
!
        ELSE
          I = ISRDM(1,NS)
          J = ISRDM(3,NS)
          K1X = ISRDM(5,NS)
          K2X = ISRDM(6,NS)
!          RDW = SRX(3)
!          N = ND(I,J,K1X)
!          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. (I.EQ.1) ) THEN
!            RDE = XGR(I+1,J,K1X)
!          ELSE
!            RDE = SQRT( AFZ(NSZ(N))/GPI )
!          ENDIF
!          PLW = PLWB(2,NS)
!          DO 520 K = K1X,K2X
!            N = ND(I,J,K)
!            IZN = IZ(N)
!            PX = MAX( PLW+PATM,PGW(2,NS)+PATM )
!            CALL WATLQD( T(2,N),PX,RHOWL )
!            BUOY = RHOWL*GRAVZ*DZGF(N)
!            AF = MIN( (PL(2,N)-PGW(2,NS))/BUOY,1.D+0 )
!            PLWP = MAX( PLW-5.D-1*BUOY,PGW(2,NS) )
!            DPLWPL = PLWP-PL(2,N)
!            IF( AF/EPSL.GT.EPSL ) THEN
!              PERMX = AF*SQRT( PERM(1,IZN)*PERM(2,IZN) )
!              DRD2 = (RDE**2-RDW**2)
!              HCFX = 2.D+0*GPI*PERMX*DRD2*DZGF(N)/
!     &          (VISL(2,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
!              QWS = QWS + HCFX*DPLWPL
!              PX = PLWP+PATM
!              CALL WATLQD( T(2,N),PX,RHOWL )
!              BUOY = RHOWL*GRAVZ*DZGF(N)
!              PLW = PLW-BUOY
!            ELSE
!              HCFX = 0.D+0
!            ENDIF
!            HCI = HCI + HCFX
!            HCF(K) = HCFX
!  520     CONTINUE
          PLW = PLWB(2,NS)
          DO 550 K = K1X,K2X
            N = ND(I,J,K)
            MP = IXP(N)
            IF( ILES.EQ.1 ) THEN
              MCOL = MP
              MROW = MDT
            ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
              MA = 1
              MCOL = KLUC(MP,MA)
              MA = MA + 1
            ENDIF
            IZN = IZ(N)
            PX = MAX( PLW+PATM,PGW(2,NS)+PATM )
            CALL WATLQD( T(2,N),PX,RHOWL )
            BUOY = RHOWL*GRAVZ*DZGF(N)
            PLWP = MAX( PLW-5.D-1*BUOY,PGW(2,NS) )
            IF( HCF(K,NS).LE.ZERO ) THEN
              SRCIC(N,NSL) = SRCIC(N,NSL) + C(N,NSL)*YL(N,NSL)*DT
     &          *HCF(K,NS)/(SL(2,N)*PORD(2,N))
            ELSE
              SRCIC(N,NSL) = SRCIC(N,NSL) + SRX(8+NS)*DT*HCF(K,NS)
            ENDIF
            PX = PLWP+PATM
            CALL WATLQD( T(2,N),PX,RHOWL )
            BUOY = RHOWL*GRAVZ*DZGF(N)
            PLW = PLW-BUOY
  550     CONTINUE
        ENDIF
#endif
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORIT4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
#ifdef awtap
      SUBROUTINE SORT4( HCFL,HCFN,NSL )
#else
      SUBROUTINE SORT4( NSL )
#endif
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Water-Oil Mode
!
!     Compute solute transport source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!
!----------------------Type Declarations-------------------------------!
!
#ifdef petsc
      INTEGER :: IROW
      REAL*8 :: BUFFER
      PetscErrorCode :: IERR
#endif
      REAL*8 SRX(8)
#ifdef awtap
      REAL*8 HCFL(LFZ,LSR),HCFN(LFZ,LSR)
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORT4'
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
#ifdef petsc
!
!---    Flush MatSetValues to switch between ADD_VALUES and 
!       INSERT_VALUES  ---
!
        CALL MatAssemblyBegin(T_MAT,MAT_FLUSH_ASSEMBLY,IERR)
        CALL MatAssemblyEnd(T_MAT,MAT_FLUSH_ASSEMBLY,IERR)
#endif
#ifdef future_coding
        IF( ISRT(NS).LT.20 .OR. ISRT(NS).GT.29 ) THEN
#endif
!
!---  Loop over source domain  ---
!
        DO 500 I = ISRDM(1,NS),ISRDM(2,NS)
        DO 500 J = ISRDM(3,NS),ISRDM(4,NS)
        DO 500 K = ISRDM(5,NS),ISRDM(6,NS)
          N = ND(I,J,K)
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
!---  Aqueous Volumetric Sink  ---
!
          IF( ISRT(NS).EQ.3 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YL(N,NSL)/(PORD(2,N)*SL(2,N))
!
!---  NAPL Volumetric Sink  ---
!
          ELSEIF( ISRT(NS).EQ.6 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YN(N,NSL)/(PORD(2,N)*SN(2,N))
!
!---  Aqueous Mass Sink  ---
!
          ELSEIF( ISRT(NS).EQ.7 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YL(N,NSL)/(SL(2,N)*PORD(2,N)*RHOL(2,N))
!
!---  NAPL Mass Sink  ---
!
          ELSEIF( ISRT(NS).EQ.10 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YN(N,NSL)/(RHON(2,N)*PORD(2,N)*SN(2,N))
!
!---  Solute source  ---
!
          ELSEIF( ISRT(NS).EQ.-NSL ) THEN
            BLU(MP) = BLU(MP) + SRX(4)
!
!---  Solute density source  ---
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
     &                      ADD_VALUES,IERR)
#endif
          ENDIF
  500   CONTINUE
#ifdef future_coding
!
!---  Well source  ---
!
        ELSE
          I = ISRDM(1,NS)
          J = ISRDM(3,NS)
          K1X = ISRDM(5,NS)
          K2X = ISRDM(6,NS)
!          RDW = SRX(3)
!          N = ND(I,J,K1X)
!          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. (I.EQ.1) ) THEN
!            RDE = XGR(I+1,J,K1X)
!          ELSE
!            RDE = SQRT( AFZ(NSZ(N))/GPI )
!          ENDIF
!          PLW = PLWB(2,NS)
!          DO 520 K = K1X,K2X
!            N = ND(I,J,K)
!            IZN = IZ(N)
!            PX = MAX( PLW+PATM,PGW(2,NS)+PATM )
!            CALL WATLQD( T(2,N),PX,RHOWL )
!            BUOY = RHOWL*GRAVZ*DZGF(N)
!            AF = MIN( (PL(2,N)-PGW(2,NS))/BUOY,1.D+0 )
!            PLWP = MAX( PLW-5.D-1*BUOY,PGW(2,NS) )
!            DPLWPL = PLWP-PL(2,N)
!            IF( AF/EPSL.GT.EPSL ) THEN
!              PERMX = AF*SQRT( PERM(1,IZN)*PERM(2,IZN) )
!              DRD2 = (RDE**2-RDW**2)
!              HCFX = 2.D+0*GPI*PERMX*DRD2*DZGF(N)/
!     &          (VISL(2,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
!              QWS = QWS + HCFX*DPLWPL
!              PX = PLWP+PATM
!              CALL WATLQD( T(2,N),PX,RHOWL )
!              BUOY = RHOWL*GRAVZ*DZGF(N)
!              PLW = PLW-BUOY
!            ELSE
!              HCFX = 0.D+0
!            ENDIF
!            HCI = HCI + HCFX
!            HCF(K) = HCFX
!  520     CONTINUE
          PLW = PLWB(2,NS)
          DO 550 K = K1X,K2X
            N = ND(I,J,K)
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
            IZN = IZ(N)
            PX = MAX( PLW+PATM,PGW(2,NS)+PATM )
            CALL WATLQD( T(2,N),PX,RHOWL )
            BUOY = RHOWL*GRAVZ*DZGF(N)
            PLWP = MAX( PLW-5.D-1*BUOY,PGW(2,NS) )
            SORTX = 0.D+0
            IF( HCF(K,NS).LE.ZERO ) THEN
              SORTX = -YL(N,NSL)*HCF(K,NS)/(SL(2,N)*PORD(2,N))
            ELSE
              BLU(MP) = BLU(MP) + SRX(8+NS)*HCF(K,NS)
            ENDIF
            PX = PLWP+PATM
            CALL WATLQD( T(2,N),PX,RHOWL )
            BUOY = RHOWL*GRAVZ*DZGF(N)
            PLW = PLW-BUOY
!
!---        Load Jacobian  ---
!
            IF( ILES.EQ.1 ) THEN
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + SORTX
            ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
              DLU(MCOL) = DLU(MCOL) + SORTX
#ifdef petsc
            ELSEIF( ILES.EQ.5 ) THEN
              CALL MatSetValues(T_MAT,1,IROW,1,IROW,SORTX,
     &                      ADD_VALUES,IERR)
#endif
            ENDIF
  550     CONTINUE
        ENDIF
#endif
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORT4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPRP4( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Calculates the aqueous-phase and NAPL solute
!     mole fractions from user-specified partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January 1995.
!     Last Modified by MD White, Battelle, PNL, January 19, 1995.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/SPRP4'
!
!---  Loop over all nodes  ---
!
      DO 900 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 900
        IZN = IZ(N)
        SVLX = PORD(2,N)*SL(2,N)+SMALL
        SVNX = PORD(2,N)*SN(2,N)+SMALL
!
!-------  Linear Kd sorption  ---
!
        IF( IPCSL(IZN,NSL).EQ.1 ) THEN
          PCSLX = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
!
!-------  Linear, Freundlich, or Langmuir sorption  ---
!
        ELSE
          PCSLX = PCSL(1,IZN,NSL)
        ENDIF
!
!---------  Aqueous-NAPL adsorption linear isotherm and
!           Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
        CTX = MAX( C(N,NSL),EPSL )
        IF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
          CALL LN_LN( CTX,CNX,CLX,SVLX,SVNX,N,IZN,NSL )
!
!---------  Aqueous-NAPL adsorption Freundlich isotherm and
!           Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
        ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
          CALL FR_LN( CTX,CNX,CLX,SVLX,SVNX,N,IZN,NSL )
!
!---------  Aqueous-NAPL adsorption Langmuir isotherm and
!           Solid-Aqueous adsorption linear or linear Kd isotherm  ---
!
        ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).LE.1 ) THEN
          CALL LG_LN( CTX,CNX,CLX,SVLX,SVNX,N,IZN,NSL )
!
!---------  Aqueous-NAPL adsorption linear isotherm and
!           Solid-Aqueous adsorption Freundlich isotherm  ---
!
        ELSEIF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
          CALL LN_FR( CTX,CNX,CLX,SVLX,SVNX,N,IZN,NSL )
!
!---------  Aqueous-NAPL adsorption Freundlich isotherm and
!           Solid-Aqueous adsorption Freundlich isotherm  ---
!
        ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
          CALL FR_FR( CTX,CNX,CLX,SVLX,SVNX,N,IZN,NSL )
!
!---------  Aqueous-NAPL adsorption Langmuir isotherm and
!           Solid-Aqueous adsorption Fruendlich isotherm  ---
!
        ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).EQ.2 ) THEN
          CALL LG_FR( CTX,CNX,CLX,SVLX,SVNX,N,IZN,NSL )
!
!---------  Aqueous-NAPL adsorption linear isotherm and
!           Solid-Aqueous adsorption Langmuir isotherm  ---
!
        ELSEIF( IPCLN(NSL).EQ.0 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
          CALL LN_LG( CTX,CNX,CLX,SVLX,SVNX,N,IZN,NSL )
!
!---------  Aqueous-NAPL adsorption Freundlich isotherm and
!           Solid-Aqueous adsorption Langmuir isotherm  ---
!
        ELSEIF( IPCLN(NSL).EQ.2 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
          CALL FR_LG( CTX,CNX,CLX,SVLX,SVNX,N,IZN,NSL )
!
!---------  Aqueous-NAPL adsorption Langmuir isotherm and
!           Solid-Aqueous adsorption Langmuir isotherm  ---
!
        ELSEIF( IPCLN(NSL).EQ.3 .AND. IPCSL(IZN,NSL).EQ.3 ) THEN
          CALL LG_LG( CTX,CNX,CLX,SVLX,SVNX,N,IZN,NSL )
        ENDIF
        YL(N,NSL) = SVLX*CLX/CTX
        YN(N,NSL) = SVNX*CNX/CTX
        YG(N,NSL) = MAX( 1.D+0-YL(N,NSL)-YN(N,NSL),0.D+0 )
  900 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SPRP4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TMPR4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
!     Compute thermo-physical properties for aqueous, napl, and gas
!     phases.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     Last Modified by MD White, PNNL, November 3, 1999.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/TMPR4'
!
!---  Loop over nodes and secondary variable indices  ---
!
      DO 300 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 300
        IZN = IZ(N)
        DO 200 M = 2,ISVC+2
          PLX = PL(M,N) + PATM
          PNX = PN(M,N) + PATM
          PGX = PG(2,N) + PATM
          TX = T(2,N)
!
!---      Aqueous phase w/ dissolved oil
!         Unknowns: aqueous pressure, dissolved-oil mole fraction  ---
!
          IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.2 ) THEN
            PCAP = MAX( PGX-PLX,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),TX,RHORL,RCW,PVW(M,N) )
            ENDIF
            PVO(M,N) = XMLO(M,N)*HCOW
            PX = MAX( PLX,PGX,PNX )
!
!---      Aqueous phase w/entrapped NAPL
!         Unknowns: aqueous pressure, trapped NAPL saturation  ---
!
          ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.4 ) THEN
            PCAP = MAX( PGX-PLX,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),TX,RHORL,RCW,PVW(M,N) )
            ENDIF
            PVO(M,N) = PSO(2,N)
            XMLO(M,N) = PVO(M,N)/HCOW
!            IF( XMLO(M,N).LT.EPSL ) XMLO(M,N) = 0.D+0
!            IF( SN(M,N).LT.EPSL ) XMLO(M,N) = 0.D+0
            PX = MAX( PLX,PGX,PNX )
!
!---      Aqueous and NAPL phases
!         Unknowns: aqueous pressure, NAPL pressure  ---
!
          ELSEIF( NPHAZ(2,N).EQ.5 .OR. NPHAZ(2,N).EQ.6 ) THEN
            PCAP = MAX( PNX-PLX,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),TX,RHORL,RCW,PVW(M,N) )
            ENDIF
            PVO(M,N) = PSO(2,N)
            XMLO(M,N) = PVO(M,N)/HCOW
!            IF( XMLO(M,N).LT.EPSL ) XMLO(M,N) = 0.D+0
!            IF( SN(M,N).LT.EPSL ) XMLO(M,N) = 0.D+0
            PX = MAX( PLX,PGX,PNX )
          ENDIF
!
!---      Compute the aqueous-phase mass fractions  ---
!
          XMLW(M,N) = MAX( 1.D+0-XMLO(M,N),ZERO )
          IF( XMLW(M,N).LT.EPSL ) XMLW(M,N) = 0.D+0
          XLO(M,N) = XMLO(M,N)*WTMO/(XMLO(M,N)*WTMO + XMLW(M,N)*WTMW)
          IF( XLO(M,N).LT.EPSL ) XLO(M,N) = 0.D+0
          XLW(M,N) = MAX( 1.D+0-XLO(M,N),ZERO )
          IF( XLW(M,N).LT.EPSL ) XLW(M,N) = 0.D+0
!
!---      Compute the liquid water density and viscosity  ---
!
          CALL WATLQD( TX,PX,RHOL(M,N) )
          RHOML(M,N) = RHOL(M,N)/(XMLW(M,N)*WTMW + XMLO(M,N)*WTMO)
          CALL WATLQV( TX,PX,PSW(2,N),VISL(M,N) )
!
!---      Compute soil mass fraction of sorbed oil  ---
!
          CLOX = XLO(M,N)*RHOL(M,N)
          RHOBX = RHOS(IZN)*(1.D+0-PORT(M,N))
          CALL SORB( CLOX,RHOBX,XSO(M,N),IZN )
!
!---      Compute the NAPL density and viscosity  ---
!
          CALL VOCLQD( TX,PX,PVO(M,N),RHON(M,N) )
          CALL VOCLQV( TX,VISN(M,N),VISCO,IVISO )
!
!---      Compute liquid diffusion coefficients  ---
!
          CALL VOCDFL( TX,VISL(M,N),DFLO(M,N) )
  200   CONTINUE
  300 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TMPR4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRPNL4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil Mode
!
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
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      IF( (IEO.EQ.2) .AND. (NSTEP-NRST.EQ.0) ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TRPNL4'
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
        DFMLX = SQRT( ABS(UL(1,NPX))*ABS(UL(1,NQX)) +
     &    ABS(VL(1,NPY))*ABS(VL(1,NQY)) +
     &    ABS(WL(1,NPZ))*ABS(WL(1,NQZ)))
        ULX = 5.D-1*(UL(1,NPX)+UL(1,NQX))
        VLX = 5.D-1*(VL(1,NPY)+VL(1,NQY))
        WLX = 5.D-1*(WL(1,NPZ)+WL(1,NQZ))
        ULGX = 5.D-1*(UL(1,NPX)*GRVX(NPX)+UL(1,NQX)*GRVX(NQX))
        VLGX = 5.D-1*(VL(1,NPY)*GRVY(NPY)+VL(1,NQY)*GRVY(NQY))
        WLGX = 5.D-1*(WL(1,NPZ)*GRVZ(NPZ)+WL(1,NQZ)*GRVZ(NQZ))
        DFALX = ((ULGX + VLGX + WLGX)/GRAV)/
     &    ( SQRT( ULX**2 + VLX**2 + WLX**2 ) + SMALL )
        SKL = SQRT((PERM(1,IZN)*(5.D-1*(GRVX(NPX)+GRVX(NQX))))**2 + 
     &    (PERM(2,IZN)*5.D-1*(GRVY(NPY)+GRVY(NQY)))**2 +
     &    (PERM(3,IZN)*5.D-1*(GRVZ(NPZ)+GRVZ(NQZ)))**2)/GRAV
        BNL = SCHR(8,IZN)/(SCHR(6,IZN)*CA_NL)
        RKLMX = (RKL(1,2,N)*RKL(2,2,N)*RKL(3,2,N))**(1./3.)
        BNDX = (RHOL(2,N)-RHON(2,N))*GRAV*SKL*RKLMX*BNL/SCHR(8,IZN)
        CAPX = DFMLX*VISL(2,N)*BNL/SCHR(8,IZN)
        TRPNL(2,N) = SQRT(CAPX**2 + 2.D+0*CAPX*BNDX*DFALX + BNDX**2)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TRPNL4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT4
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Water-Oil Mode
!
!     Update the primary variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_wo.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
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
      SUB_LOG(ISUB_LOG) = '/UPDT4'
      IF( ICNV.EQ.1 ) GOTO 300
      IERR = 0
!
!---  Update primary variables
!
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 200
        NMD = IXP(N)
        MPL = IM(IEQW,NMD)
        MPN = IM(IEQO,NMD)
        IZN = IZ(N)
!
!---    Aqueous w/ dissolved oil system prior to iteration  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          IF( XMLO(2,N)+BLU(MPN).LT.ZERO ) THEN
            XMLO(2,N) = MAX( XMLO(2,N)+6.D-1*BLU(MPN),ZERO )
          ELSE
            XMLO(2,N) = MAX( XMLO(2,N)+BLU(MPN),ZERO )
          ENDIF
          IF( XMLO(2,N).LT.EPSL ) XMLO(2,N) = 0.D+0
!
!---    Aqueous-gas w/ dissolved oil system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          IF( XMLO(2,N)+BLU(MPN).LT.ZERO ) THEN
            XMLO(2,N) = MAX( XMLO(2,N)+6.D-1*BLU(MPN),ZERO )
          ELSE
            XMLO(2,N) = MAX( XMLO(2,N)+BLU(MPN),ZERO )
          ENDIF
          IF( XMLO(2,N).LT.EPSL ) XMLO(2,N) = 0.D+0
        ENDIF
!
!---    Fluid-pair interfacial tension and scaling factor   ---
!
        XMLOX = MIN( XMLO(2,N),PSO(2,N)/HCOW )
        XMLWX = MAX( 1.D+0-XMLOX,0.D+0 )
        XLOX = XMLOX*WTMO/(XMLOX*WTMO + XMLWX*WTMW)
        CALL SFT4( XLOX,IZN )
        ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---    Aqueous w/ dissolved oil system prior to iteration  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          IF( ABS(BLU(MPL)).GT.1.D+6 ) IERR = 1
          IF( BGL*(PG(2,N)-PL(2,N)-BLU(MPL)).GT.ENPR ) THEN
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MPL)
          ELSE
            PL(2,N) = PL(2,N) + BLU(MPL)
          ENDIF
!
!---    Aqueous-gas w/ dissolved oil system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          IF( ABS(BLU(MPL)).GT.1.D+6 ) IERR = 1
          IF( BGL*(PG(2,N)-PL(2,N)-BLU(MPL)).LT.ENPR ) THEN
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MPL)
          ELSE
            PL(2,N) = PL(2,N) + BLU(MPL)
          ENDIF
!
!---    Aqueous w/ trapped NAPL system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
          IF( ABS(BLU(MPL)).GT.1.D+6 ) IERR = 1
          IF( BGL*(PG(2,N)-PL(2,N)-BLU(MPL)).GT.ENPR ) THEN
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MPL)
          ELSE
            PL(2,N) = PL(2,N) + BLU(MPL)
          ENDIF
          IF( SN(2,N)+BLU(MPN).LT.ZERO ) THEN
            SN(2,N) = MAX( SN(2,N)+6.D-1*BLU(MPN),ZERO )
          ELSEIF( SN(2,N)+BLU(MPN).GT.(1.D+0-SCHR(4,IZN)) ) THEN
            SN(2,N) = MIN( SN(2,N)+6.D-1*BLU(MPN),(1.D+0-SCHR(4,IZN)) )
          ELSE
            SN(2,N) = MAX( SN(2,N)+BLU(MPN),ZERO )
            SN(2,N) = MIN( SN(2,N),(1.D+0-SCHR(4,IZN)) )
          ENDIF
          IF( SN(2,N).LT.EPSL ) SN(2,N) = 0.D+0
!
!---    Aqueous-gas w/ trapped NAPL system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
          IF( ABS(BLU(MPL)).GT.1.D+6 ) IERR = 1
          IF( BGL*(PG(2,N)-PL(2,N)-BLU(MPL)).LT.ENPR ) THEN
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MPL)
          ELSE
            PL(2,N) = PL(2,N) + BLU(MPL)
          ENDIF
          IF( SN(2,N)+BLU(MPN).LT.ZERO ) THEN
            SN(2,N) = MAX( SN(2,N)+6.D-1*BLU(MPN),ZERO )
          ELSEIF( SN(2,N)+BLU(MPN).GT.(1.D+0-SCHR(4,IZN)) ) THEN
            SN(2,N) = MIN( SN(2,N)+6.D-1*BLU(MPN),(1.D+0-SCHR(4,IZN)) )
          ELSE
            SN(2,N) = MAX( SN(2,N)+BLU(MPN),ZERO )
            SN(2,N) = MIN( SN(2,N),(1.D+0-SCHR(4,IZN)) )
          ENDIF
          IF( SN(2,N).LT.EPSL ) SN(2,N) = 0.D+0
!
!---    Aqueous-NAPL system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
          IF( ABS(BLU(MPL)).GT.1.D+6 ) BLU(MPL) = 6.D-1*BLU(MPL)
          IF( ABS(BLU(MPN)).GT.1.D+6 ) BLU(MPN) = 6.D-1*BLU(MPN)
          IF( ABS(BLU(MPL)).GT.1.D+6 ) IERR = 1
          IF( ABS(BLU(MPN)).GT.1.D+6 ) IERR = 1
          PLX = PL(2,N) + BLU(MPL)
          PNX = PN(2,N) + BLU(MPN)
          PNC = MAX( (BNL*PLX+BGN*PG(2,N))/(BNL+BGN),
     &      (PLX*BNL+ENPR)/BNL )
          IF( PNX.LT.PNC ) THEN
            PN(2,N) = PN(2,N) + 6.D-1*BLU(MPN)
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MPL)
          ELSE
            PN(2,N) = PNX
            PL(2,N) = PLX
          ENDIF
          PL(2,N) = MAX( PL(2,N),PN(2,N)-1.D+5 )
!
!---    Aqueous-NAPL-Gas system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
          IF( ABS(BLU(MPL)).GT.1.D+6 ) BLU(MPL) = 6.D-1*BLU(MPL)
          IF( ABS(BLU(MPN)).GT.1.D+6 ) BLU(MPN) = 6.D-1*BLU(MPN)
          IF( ABS(BLU(MPL)).GT.1.D+6 ) IERR = 1
          IF( ABS(BLU(MPN)).GT.1.D+6 ) IERR = 1
          PLX = PL(2,N) + BLU(MPL)
          PNX = PN(2,N) + BLU(MPN)
          PNC = MAX( (BNL*PLX+BGN*PG(2,N))/(BNL+BGN),
     &      (PLX*BNL+ENPR)/BNL )
          IF( PNX.LT.PNC ) THEN
            PN(2,N) = PN(2,N) + 6.D-1*BLU(MPN)
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MPL)
          ELSE
            PN(2,N) = PNX
            PL(2,N) = PLX
          ENDIF
!          PL(2,N) = MAX( PL(2,N),PN(2,N)-1.D+5 )
        ENDIF
!
!---  Reduce time step for excessive changes in primary variables   ---
!
        IF( IERR.EQ.1 ) THEN
          ICNV = 1
          WRITE(ISC,'(A)') 'Execution Note: Excessive change'
          WRITE(IWR,'(A)') 'Execution Note: Excessive change'
          WRITE(ISC,'(A,I6,$)') '  Node: ',N
          WRITE(IWR,'(A,I6,$)') '  Node: ',N
          WRITE(ISC,'(A,I1)') '  Phase: ',NPHAZ(2,N)
          WRITE(IWR,'(A,I1)') '  Phase: ',NPHAZ(2,N)
          WRITE(ISC,'(A,1PE12.5,$)') '  PL: ',PL(2,N)+PATM
          WRITE(IWR,'(A,1PE12.5,$)') '  PL: ',PL(2,N)+PATM
          WRITE(ISC,'(A,1PE12.5,$)') '  DPL: ',BLU(MPL)
          WRITE(IWR,'(A,1PE12.5,$)') '  DPL: ',BLU(MPL)
          WRITE(ISC,'(A,1PE12.5,$)') '  PN: ',PN(2,N)+PATM
          WRITE(IWR,'(A,1PE12.5,$)') '  PN: ',PN(2,N)+PATM
          WRITE(ISC,'(A,1PE12.5)') '  DPN: ',BLU(MPN)
          WRITE(IWR,'(A,1PE12.5)') '  DPN: ',BLU(MPN)
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
       WRITE(ISC,'(A,1PE11.4,1X,$)') '  Time step reduced from ',VARX
       WRITE(ISC,'(2A,1PE11.4,1X,A)') UNTM(1:NCH),' to ',VAR,UNTM(1:NCH)
       WRITE(IWR,'(A,1PE11.4,1X,$)') '  Time step reduced from ',VARX
       WRITE(IWR,'(2A,1PE11.4,1X,A)') UNTM(1:NCH),' to ',VAR,UNTM(1:NCH)
          DO 400 N = 1,NFLD
            PL(2,N) = PL(1,N)
            PN(2,N) = PN(1,N)
            XMLO(2,N) = XMLO(1,N)
            SN(2,N) = SN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  400     CONTINUE
#ifdef awtap
          DO 402 N = 1,NSR
            PLWB(2,N) = PLWB(1,N)
            PGW(2,N) = PGW(1,N)
  402     CONTINUE
#endif
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
          DO 410 N = 1,NFLD
            PL(2,N) = PL(1,N)
            PN(2,N) = PN(1,N)
            XMLO(2,N) = XMLO(1,N)
            SN(2,N) = SN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  410     CONTINUE
#ifdef awtap
          DO 412 N = 1,NSR
            PLWB(2,N) = PLWB(1,N)
            PGW(2,N) = PGW(1,N)
  412     CONTINUE
#endif
          NSTEP = NSTEP - 1
          TM = TM - DT
          WRITE(ISC,'(A)') 'Execution Note:  Reduction Limit Exceeded'
          WRITE(IWR,'(A)') 'Execution Note:  Reduction Limit Exceeded'
          ICNV = 4
        ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT4 group
!
      RETURN
      END

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
!---------------------------------------------------------!
!
!     STOMP: Subsurface Transport Over Multiple Phases
!
!     Water-NComponent-NaCl-Energy Mode (STOMP-WNSE)
!
!     This engineering program numerically simulates the transport
!     of H2O, NaCl and CO2 through multifluid subsurface environments
!     under nonisothermal conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(1) = 'STOMP_COMP'
      ICODE = 40
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
      CALL RDINPT_COMP
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
      CALL CHK_COMP
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
      CALL INCRM_COMP
!
!---  Initial hydrologic and thermodynamic properties on interior nodes
!     and boundary surfaces  ---
!
       CALL PROP_COMP
       CALL BCP_COMP
!
!---  Compute initial solute concentrations  ---
!
      CALL CISC_COMP
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
!---  Compute initial fluxes on interior and boundary surfaces  ---
!
      ISVF = 1
      CALL DRCVL
      CALL DRCVG
      CALL DFFGW_COMP
      DO IGC = 1,NGC 
       CALL DFFLC(IGC)
       CALL DFFGC(IGC)
      ENDDO
!
!---  Isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS_COMP
!
!---  Nonisothermal simulations  ---
!
      IF( ISLC(30).EQ.0 ) THEN
        CALL THD_COMP
        CALL THAL_COMP
        CALL THAG_COMP
        CALL THDG_COMP
      ENDIF
      CALL BCF_COMP
      ISVF = 2*ISVC+1
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
      CALL LDO_COMP
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
        CHMSG = 'Simulation Stopped:  Simulation Time Limit'
        CALL WRMSGS( INDX )
        GOTO 900
      ENDIF
!
!---  Stop simulation file "stop_stomp" exists  ---
!
      INQUIRE( FILE="stop_stomp", EXIST=HALT )
      IF( HALT ) THEN
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
!---  Generate plot if "plot_stomp" exists  ---
!
      INQUIRE( FILE="plot_stomp", EXIST=PLOT )
      IF( PLOT ) THEN
        OPEN( UNIT=19, FILE="plot_stomp" )
        CLOSE( UNIT=19, STATUS='DELETE' )
        CALL WRPLOT
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
        CALL BCP_COMP
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
      CALL BCP_COMP
!
!---  Compute source contributions  ---
!
      CALL SORC_COMP
!
!---  Compute aqueous-phase volumetric flux (interior surfaces)  ---
!
      CALL DRCVL
!
!---  Compute gas-phase volumetric flux (interior surfaces)  ---
!
      CALL DRCVG
!
!---  Compute water vapor diffusion flux through the gas phase
!     (interior surfaces)  ---
!
      CALL DFFGW_COMP
!
!---  Compute dissolved-CO2 diffusion flux through the aqueous phase
!     (interior surfaces)  ---
!
      DO IGC = 1, NGC
       CALL DFFLC(IGC)
       CALL DFFGC(IGC)
      ENDDO
!
!---  Compute dissolved-NaCl diffusion flux through the aqueous phase
!     (interior surfaces), isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS_COMP
!
!---  Nonisothermal simulations  ---
!
      IF( ISLC(30).EQ.0 ) THEN
!
!---    Compute diffusive thermal flux (interior surfaces)  ---
!
        CALL THD_COMP
!      STOP
!
!---    Compute aqueous-phase advective thermal flux
!       (interior surfaces)  ---
!
        CALL THAL_COMP
!
!---    Compute gas-phase advective thermal flux
!       (interior surfaces)  ---
!
        CALL THAG_COMP
!
!---    Compute vapor diffusive thermal flux
!       (interior surfaces)  ---
!
        CALL THDG_COMP
      ENDIF
!
!---  Compute aqueous-phase volumetric flux, gas-phase volumetric flux,
!     water vapor mass flux, diffusive thermal flux, aqueous-phase
!     advective thermal flux, gas-phase advective thermal flux,
!     vapor diffusive thermal flux, aqueous-phase salt flux
!     (boundary surfaces)  ---
!
      CALL BCF_COMP
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
      IF( ISLC(30).EQ.0 ) CALL JCBT_COMP
!
!---  Load Jacobian matrix for the water equation
!     (zero flux boundary)  ---
!
      CALL JCBW_COMP
!
!---  Load Jacobian matrix for the salt equation
!     (zero flux boundary), isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL JCBS_COMP
!
!---  Load Jacobian matrix for the CO2 equation
!     (zero flux boundary)  ---
!
      DO IGC = 1, NGC
       CALL JCBC_COMP(IGC)
      ENDDO
!
!---  Modify the Jacobian matrix for boundary conditions  ---
!
      CALL BCJ_COMP
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
      CALL UPDT_COMP
!
!---  Compute convergence from maximum relative residuals  ---
!
      CALL RSDL_COMP
!
!---  Compute primary variable increments, saturation,
!     relative permeability, porosity, tortuosity,
!     thermodynamic properties for interior nodes,
!     except immediately after a new time step  ---
!
      CALL INCRM_COMP
      CALL PROP_COMP
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
!---  Integrate coupled-equation sources ---
!
      CALL SORIC_COMP
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
!---  Compute water vapor diffusion flux through the gas phase
!     (interior surfaces)  ---
!
      CALL DFFGW_COMP
!
!---  Compute dissolved CO2 diffusion flux through the aqueous phase
!     (interior surfaces)  ---
!
      DO IGC = 1, NGC
        CALL DFFLC(IGC)
        CALL DFFGC(IGC)
      ENDDO
!
!---  Compute dissolved-NaCl diffusion flux through the aqueous phase
!     (interior surfaces), isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS_COMP
!
!---  Nonisothermal simulations  ---
!
      IF( ISLC(30).EQ.0 ) THEN
!
!---    Compute diffusive thermal flux (interior surfaces)  ---
!
        CALL THD_COMP
!
!---    Compute aqueous-phase advective thermal flux
!       (interior surfaces)  ---
!
        CALL THAL_COMP
!
!---    Compute gas-phase advective thermal flux
!       (interior surfaces)  ---
!
        CALL THAG_COMP
!
!---    Compute vapor diffusive thermal flux
!       (interior surfaces)  ---
!
        CALL THDG_COMP
      ENDIF
!
!---  Compute aqueous-phase volumetric flux, gas-phase volumetric flux,
!     water vapor mass flux, diffusive thermal flux, aqueous-phase
!     advective thermal flux, gas-phase advective thermal flux,
!     vapor diffusive thermal flux (boundary surfaces)  ---
!
      CALL BCF_COMP
!
!---  Compute Local Courant Numbers  ---
!
      IF( ICRNT.EQ.1 ) CALL CRNTNB
      ISVF = 2*ISVC+1
!
!---  Beginning of transport equation solution  ---
!
      IF( IEQC.EQ.0 .AND. ISLC(40).EQ.0 ) GOTO 800
!
!---  Loop over number of solutes  ---
!
      DO 700 NSL = 1,NSOLU
!
!---  Courant number limiting  ---
!
        N_CRN(NSL) = 1
        IF( ISLC(17).NE.0 ) CALL CRN_LIM( NSL )
!
!---    Sub-time step loop  ---
!
        DO 690 NC = 1,N_CRN(NSL)
          IF( ISLC(17).NE.0 ) TM = MIN( TM+DT,TM_CRN )
!
!---      Compute solute mole fractions ---
!
          CALL SPRP_COMP( NSL )
!
!---      Solute transport ---
!
          CALL TPORT_COMP( NSL )
!
!---      Load old sub-time-step concentrations  ---
!
          IF( ISLC(17).NE.0 ) CALL UPDTCO( NSL)
!
!---    Bottom of sub-time step loop  ---
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
        N_CRN(NSOLU+1) = 1
        IF( ISLC(17).NE.0 ) CALL CRN_LIM( NSOLU+1 )
!
!---    Courant-limiting sub-time step loop  ---
!
        DO 792 NCR = 1,N_CRN(NSOLU+1)
          IF( ISLC(17).NE.0 ) TM = MIN( TM+DT,TM_CRN )
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
!---        Skip transport for linked aqueous component   ---
!
            DO IGC = 1, NGC
            IF( ISPLK(14+NSPLK+IGC).EQ.NSL ) GOTO 720
            END DO
!
!---        Mobile conservation component fractions   ---
!
            CALL MOBCF( NEQ )
!
!---        Solute transport ---
!
            CALL TPORT_COMP( NSL )
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
!---        Skip transport for linked aqueous CO2   ---
!
            DO IGC = 1, NGC
            IF( ISPLK(14+NSPLK+IGC).EQ.NSL ) GOTO 740
            END DO
! 
!---        Mobile kinetic component fractions   ---
!
            CALL MOBKF( NEQ )
!
!---        Solute transport ---
!
            CALL TPORT_COMP( NSL )
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
      SUBROUTINE BCF_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Compute boundary surface fluxes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXS
      USE FLUXP
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/BCF_COMP'
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
            WDGW(M,NPZ) = 0.D+0
            WQ(M,NPZ) = 0.D+0
            WS(M,NPZ) = 0.D+0
            WDS(M,NPZ) = 0.D+0
!
           DO IGC = 1, NGC
            WGC(IGC,M,NPZ) = 0.D+0
            WLC(IGC,M,NPZ) = 0.D+0
            WDGC(IGC,M,NPZ) = 0.D+0
            WDLC(IGC,M,NPZ) = 0.D+0
           ENDDO
   10     CONTINUE
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          DO 20 M = 1,ISVF
            VL(M,NPY) = 0.D+0
            VG(M,NPY) = 0.D+0
            VDGW(M,NPY) = 0.D+0
            VQ(M,NPY) = 0.D+0
            VS(M,NPY) = 0.D+0
            VDS(M,NPY) = 0.D+0
!
           DO IGC = 1, NGC
            VGC(IGC,M,NPY) = 0.D+0
            VLC(IGC,M,NPY) = 0.D+0
            VDGC(IGC,M,NPY) = 0.D+0
            VDLC(IGC,M,NPY) = 0.D+0
           ENDDO
   20     CONTINUE
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          DO 30 M = 1,ISVF
            UL(M,NPX) = 0.D+0
            UG(M,NPX) = 0.D+0
            UDGW(M,NPX) = 0.D+0
            UQ(M,NPX) = 0.D+0
            US(M,NPX) = 0.D+0
            UDS(M,NPX) = 0.D+0
!
           DO IGC = 1, NGC
            UGC(IGC,M,NPX) = 0.D+0
            ULC(IGC,M,NPX) = 0.D+0
            UDGC(IGC,M,NPX) = 0.D+0
            UDLC(IGC,M,NPX) = 0.D+0
           ENDDO
   30     CONTINUE
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          DO 40 M = 1,ISVF
            UL(M,NQX) = 0.D+0
            UG(M,NQX) = 0.D+0
            UDGW(M,NQX) = 0.D+0
            UQ(M,NQX) = 0.D+0
            US(M,NQX) = 0.D+0
            UDS(M,NQX) = 0.D+0
!
           DO IGC = 1, NGC
            UGC(IGC,M,NQX) = 0.D+0
            ULC(IGC,M,NQX) = 0.D+0
            UDGC(IGC,M,NQX) = 0.D+0
            UDLC(IGC,M,NQX) = 0.D+0
           ENDDO
   40     CONTINUE
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          DO 50 M = 1,ISVF
            VL(M,NQY) = 0.D+0
            VG(M,NQY) = 0.D+0
            VDGW(M,NQY) = 0.D+0
            VQ(M,NQY) = 0.D+0
            VS(M,NQY) = 0.D+0
!
           DO IGC = 1, NGC
            VGC(IGC,M,NQY) = 0.D+0
            VLC(IGC,M,NQY) = 0.D+0
            VDGC(IGC,M,NQY) = 0.D+0
            VDLC(IGC,M,NQY) = 0.D+0
            VDS(M,NQY) = 0.D+0
           ENDDO
   50     CONTINUE
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          DO 60 M = 1,ISVF
            WL(M,NQZ) = 0.D+0
            WG(M,NQZ) = 0.D+0
            WDGW(M,NQZ) = 0.D+0
            WQ(M,NQZ) = 0.D+0
            WS(M,NQZ) = 0.D+0
            WDS(M,NQZ) = 0.D+0
!
           DO IGC = 1, NGC
            WGC(IGC,M,NQZ) = 0.D+0
            WLC(IGC,M,NQZ) = 0.D+0
            WDGC(IGC,M,NQZ) = 0.D+0
            WDLC(IGC,M,NQZ) = 0.D+0
           ENDDO
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
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 110 M = 1,ISVF
              WL(M,NSZ(N)) = BCX(3)
  110       CONTINUE
            DO IGC = 1, NGC
              CALL DFFLCB( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLB( N,NB )
            DO IGC = 1, NGC
              CALL DFFLCB( N,NB,IGC )
            ENDDO
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(4,NB).EQ.2 ) THEN
            DO 111 M = 1,ISVF
              WG(M,NSZ(N)) = BCX(4)
  111       CONTINUE
            CALL DFFGWB_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCB( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(4,NB).NE.3 ) THEN
            CALL DRCVGB( N,NB )
            CALL DFFGWB_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCB( N,NB,IGC )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(IEQT,NB).EQ.2 ) THEN
              DO 112 M = 1,ISVF
                WQ(M,NSZ(N)) = BCX(2)
  112         CONTINUE
!
!---        Energy outflow  ---
!
            ELSEIF( IBCT(IEQT,NB).EQ.7 ) THEN
              CALL THALB_COMP( N,NB )
              CALL THAGB_COMP( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL THDB_COMP( N,NB )
              CALL THALB_COMP( N,NB )
              CALL THAGB_COMP( N,NB )
              CALL THDGB_COMP( N,NB )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSB_COMP( N,NB )
            ENDIF
          ENDIF
!
!---  South boundary
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 120 M = 1,ISVF
              VL(M,NSY(N)) = BCX(3)
  120       CONTINUE
            DO IGC = 1, NGC
              CALL DFFLCS( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLS( N,NB )
            DO IGC = 1, NGC
              CALL DFFLCS( N,NB,IGC )
            ENDDO
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(4,NB).EQ.2 ) THEN
            DO 121 M = 1,ISVF
              VG(M,NSY(N)) = BCX(4)
  121       CONTINUE
            CALL DFFGWS_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCS( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(4,NB).NE.3 ) THEN
            CALL DRCVGS( N,NB )
            CALL DFFGWS_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCS( N,NB,IGC )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(IEQT,NB).EQ.2 ) THEN
              DO 122 M = 1,ISVF
                VQ(M,NSY(N)) = BCX(2)
  122         CONTINUE
!
!---        Energy outflow  ---
!
            ELSEIF( IBCT(IEQT,NB).EQ.7 ) THEN
              CALL THALS_COMP( N,NB )
              CALL THAGS_COMP( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL THDS_COMP( N,NB )
              CALL THALS_COMP( N,NB )
              CALL THAGS_COMP( N,NB )
              CALL THDGS_COMP( N,NB )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSS_COMP( N,NB )
            ENDIF
          ENDIF
!
!---    West boundary
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 130 M = 1,ISVF
              UL(M,NSX(N)) = BCX(3)
  130       CONTINUE
            DO IGC = 1, NGC
             CALL DFFLCW( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLW( N,NB )
            DO IGC = 1, NGC
             CALL DFFLCW( N,NB,IGC )
            ENDDO
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(4,NB).EQ.2 ) THEN
            DO 131 M = 1,ISVF
              UG(M,NSX(N)) = BCX(4)
  131       CONTINUE
            CALL DFFGWW_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCW( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(4,NB).NE.3 ) THEN
            CALL DRCVGW( N,NB )
            CALL DFFGWW_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCW( N,NB,IGC )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(IEQT,NB).EQ.2 ) THEN
              DO 132 M = 1,ISVF
                UQ(M,NSX(N)) = BCX(2)
  132         CONTINUE
!
!---        Energy outflow  ---
!
            ELSEIF( IBCT(IEQT,NB).EQ.7 ) THEN
              CALL THALW_COMP( N,NB )
              CALL THAGW_COMP( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL THDW_COMP( N,NB )
              CALL THALW_COMP( N,NB )
              CALL THAGW_COMP( N,NB )
              CALL THDGW_COMP( N,NB )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSW_COMP( N,NB )
            ENDIF
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---     Aqueous Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 140 M = 1,ISVF
              UL(M,NSX(N)+1) = BCX(3)
  140       CONTINUE
            DO IGC = 1, NGC
              CALL DFFLCE( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLE( N,NB )
            DO IGC = 1, NGC
              CALL DFFLCE( N,NB,IGC )
            ENDDO
          ENDIF
!
!---     Gas Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(4,NB).EQ.2 ) THEN
            DO 141 M = 1,ISVF
              UG(M,NSX(N)+1) = BCX(4)
  141       CONTINUE
            CALL DFFGWE_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCE( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(4,NB).NE.3 ) THEN
            CALL DRCVGE( N,NB )
            CALL DFFGWE_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCE( N,NB,IGC )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---       Energy Neumann  ---
!
            IF( IBCT(IEQT,NB).EQ.2 ) THEN
              DO 142 M = 1,ISVF
                UQ(M,NSX(N)+1) = BCX(2)
  142         CONTINUE
!
!---       Energy outflow  ---
!
            ELSEIF( IBCT(IEQT,NB).EQ.7 ) THEN
              CALL THALE_COMP( N,NB )
              CALL THAGE_COMP( N,NB )
!
!---       Energy not zero flux  ---
!
            ELSEIF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL THDE_COMP( N,NB )
              CALL THALE_COMP( N,NB )
              CALL THAGE_COMP( N,NB )
              CALL THDGE_COMP( N,NB )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSE_COMP( N,NB )
            ENDIF
          ENDIF
!
!---    North boundary
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 150 M = 1,ISVF
              VL(M,NSY(N)+IFLD) = BCX(3)
  150       CONTINUE
            DO IGC = 1, NGC
              CALL DFFLCN( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLN( N,NB )
            DO IGC = 1, NGC
              CALL DFFLCN( N,NB,IGC )
            ENDDO
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(4,NB).EQ.2 ) THEN
            DO 151 M = 1,ISVF
              VG(M,NSY(N)+IFLD) = BCX(4)
  151       CONTINUE
            CALL DFFGWN_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCN( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(4,NB).NE.3 ) THEN
            CALL DRCVGN( N,NB )
            CALL DFFGWN_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCN( N,NB,IGC )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(IEQT,NB).EQ.2 ) THEN
              DO 152 M = 1,ISVF
                VQ(M,NSY(N)+IFLD) = BCX(2)
  152         CONTINUE
!
!---        Energy outflow  ---
!
            ELSEIF( IBCT(IEQT,NB).EQ.7 ) THEN
              CALL THALN_COMP( N,NB )
              CALL THAGN_COMP( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL THDN_COMP( N,NB )
              CALL THALN_COMP( N,NB )
              CALL THAGN_COMP( N,NB )
              CALL THDGN_COMP( N,NB )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSN_COMP( N,NB )
            ENDIF
          ENDIF
!
!---    Top boundary
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 160 M = 1,ISVF
              WL(M,NSZ(N)+IJFLD) = BCX(3)
  160       CONTINUE
            DO IGC = 1, NGC
              CALL DFFLCT( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLT( N,NB )
            DO IGC = 1, NGC
              CALL DFFLCT( N,NB,IGC )
            ENDDO
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(4,NB).EQ.2 ) THEN
            DO 161 M = 1,ISVF
              WG(M,NSZ(N)+IJFLD) = BCX(4)
  161       CONTINUE
            CALL DFFGWT_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCT( N,NB,IGC )
            ENDDO
          ELSEIF( IBCT(4,NB).NE.3 ) THEN
            CALL DRCVGT( N,NB )
            CALL DFFGWT_COMP( N,NB )
            DO IGC = 1,NGC
              CALL DFFGCT( N,NB,IGC )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(IEQT,NB).EQ.2 ) THEN
              DO 162 M = 1,ISVF
                WQ(M,NSZ(N)+IJFLD) = BCX(2)
  162         CONTINUE
!
!---        Energy outflow  ---
!
            ELSEIF( IBCT(IEQT,NB).EQ.7 ) THEN
              CALL THALT_COMP( N,NB )
              CALL THAGT_COMP( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL THDT_COMP( N,NB )
              CALL THALT_COMP( N,NB )
              CALL THAGT_COMP( N,NB )
              CALL THDGT_COMP( N,NB )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLST_COMP( N,NB )
            ENDIF
          ENDIF
        ENDIF
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCF_COMP group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCJ_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Modify the Jacobian matrix for boundary conditions
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/BCJ_COMP'
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
!---     Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWB_COMP( N,NB,NPZ )
            DO IGC = 1, NGC
              CALL JCBLCB_COMP( IGC,N,NB,NPZ )
            ENDDO
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBGWB_COMP( N,NB,NPZ )
            DO IGC = 1, NGC
              CALL JCBGCB_COMP( IGC,N,NB,NPZ )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL JCBTB_COMP( N,NB,NPZ )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSB_COMP( N,NB,NPZ )
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWS_COMP( N,NB,NPY )
            DO IGC = 1, NGC
              CALL JCBLCS_COMP( IGC,N,NB,NPY )
            ENDDO
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBGWS_COMP( N,NB,NPY )
            DO IGC = 1, NGC
              CALL JCBGCS_COMP( IGC,N,NB,NPY )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL JCBTS_COMP( N,NB,NPY )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSS_COMP( N,NB,NPY )
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWW_COMP( N,NB,NPX )
            DO IGC = 1, NGC
              CALL JCBLCW_COMP( IGC,N,NB,NPX )
            ENDDO
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBGWW_COMP( N,NB,NPX )
            DO IGC = 1, NGC
              CALL JCBGCW_COMP( IGC,N,NB,NPX )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---       Energy  ---
!
            IF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL JCBTW_COMP( N,NB,NPX )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSW_COMP( N,NB,NPX )
          ENDIF
!
!---    East boundary
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWE_COMP( N,NB,NQX )
            DO IGC = 1, NGC
              CALL JCBLCE_COMP( IGC,N,NB,NQX )
            ENDDO
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBGWE_COMP( N,NB,NQX )
            DO IGC = 1, NGC
              CALL JCBGCE_COMP( IGC,N,NB,NQX )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL JCBTE_COMP( N,NB,NQX )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSE_COMP( N,NB,NQX )
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWN_COMP( N,NB,NQY )
            DO IGC = 1, NGC
              CALL JCBLCN_COMP( IGC,N,NB,NQY )
            ENDDO
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBGWN_COMP( N,NB,NQY )
            DO IGC = 1, NGC
              CALL JCBGCN_COMP( IGC,N,NB,NQY )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL JCBTN_COMP( N,NB,NQY )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSN_COMP( N,NB,NQY )
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWT_COMP( N,NB,NQZ )
            DO IGC = 1, NGC
              CALL JCBLCT_COMP( IGC,N,NB,NQZ )
            ENDDO
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(4,NB).NE.3 ) THEN
            CALL JCBGWT_COMP( N,NB,NQZ )
            DO IGC = 1, NGC
              CALL JCBGCT_COMP( IGC,N,NB,NQZ )
            ENDDO
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(IEQT,NB).NE.3 ) THEN
              CALL JCBTT_COMP( N,NB,NQZ )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBST_COMP( N,NB,NQZ )
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCJ_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCP_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Compute saturation, relative permeability and thermodynamic
!     properties for boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE JACOB
      USE GRID
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVT
      USE BCVS
      USE BCVP
      USE BCVGC
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
      REAL*8 BCX(LBCV),YLSX(2)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      REAL*8 RHOGC(LNGC),VISGC(LNGC),XMGCX(LNGC),XLCX(LNGC),XGCX(LNGC)
      REAL*8 VISCGX(LNGC),PGCX(LNGC),XMGX(LNGC),THKGCX(LNGC),XMLCX(LNGC)
      REAL*8 XMLX(LNGC),UEGCX(LNGC),HGCBX(LNGC)
      REAL*8 XMDGCBX(LNGC),XMDLCBX(LNGC)
      REAL*8 DFGCWX(LNGC),DFGC2X(LNGC,LNGC)
      REAL*8 DHCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCP_COMP'
      
      DO IGC = 1,NGC
        RHOGC(IGC) = 0.D+0
        VISGC(IGC) = 0.D+0
        XMGCX(IGC) = 0.D+0
        XLCX(IGC) = 0.D+0
        XGCX(IGC) = 0.D+0
        VISCGX(IGC) = 0.D+0
        PGCX(IGC) = 0.D+0
        XMGX(IGC) = 0.D+0
        THKGCX(IGC) = 0.D+0
        XMLCX(IGC) = 0.D+0
        XMLX(IGC) = 0.D+0
        UEGCX(IGC) = 0.D+0
        HGCBX(IGC) = 0.D+0
        XMDGCBX(IGC) = 0.D+0
        XMDLCBX(IGC) = 0.D+0
        DFGCWX(IGC) = 0.D+0
        DHCX(IGC) = 0.D+0
      DO JGC = 1,NGC
        DFGC2X(IGC,JGC) = 0.D+0
      END DO
      END DO
!
!---  Assign values for initial condition type boundary conditions  ---
!
      IF( NSTEP-NRST.LE.1 .AND. NITER.LE.1 ) THEN
       DO 50 NB = 1,NBC
        N = IBCN(NB)
        IF( IBCT(IEQW,NB).EQ.12 .OR. IBCT(4,NB).EQ.12 ) THEN
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
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 )  THEN
          IF( IBCT(IEQT,NB).EQ.12 ) TB(1,NB) = T(2,N)
        ENDIF
        IF( IBCT(IEQW,NB).EQ.12 ) PLB(1,NB) = PL(2,N) + RHOL(2,N)*GB
        IF( IBCT(4,NB).EQ.12 ) PGB(1,NB) = PG(2,N) + RHOG(2,N)*GB
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 )  THEN
          IF( IBCT(IEQS,NB).EQ.12 ) YLSB(1,NB) = YLS(2,N)
        ENDIF
!
!---    Aqueous component intial conditions  ---
!
        IF( IBCT(5,NB).EQ.12 ) THEN
          DO 10 IGC = 1,NGC
            XLCB(IGC,1,NB) = XLC(IGC,2,N)
   10     CONTINUE
        ENDIF
!
!---    Gas component intial conditions  ---
!
        IF( IBCT(6,NB).EQ.12 ) THEN
          DO 20 IGC = 1,NGC
            XGCB(IGC,1,NB) = XGC(IGC,2,N)
   20     CONTINUE
        ENDIF
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
!---         Nonisothermal simulations  ---
!
             IF( ISLC(30).EQ.0 )  THEN
               IF( IBCT(IEQT,NB).EQ.2 ) THEN
                 BCX(2) = BCX(2)-5.D-1*DTBC*
     &             (BC(2,M,MB)-BC(2,M-1,MB))/TDBC
               ENDIF
             ELSE
               BCX(2) = T(2,N)
             ENDIF
             IF( IBCT(IEQW,NB).EQ.2 ) THEN
               BCX(3) = BCX(3)-5.D-1*DTBC*(BC(3,M,MB)-BC(3,M-1,MB))/TDBC
             ENDIF
             IF( IBCT(4,NB).EQ.2 ) THEN
               BCX(4) = BCX(4)-5.D-1*DTBC*(BC(4,M,MB)-BC(4,M-1,MB))/TDBC
             ENDIF
!
!---         Isobrine option  ---
!
             IF( ISLC(32).EQ.0 )  THEN
               IF( IBCT(IEQS,NB).EQ.2 ) THEN
                 BCX(6) = BCX(6)-5.D-1*DTBC*
     &             (BC(6,M,MB)-BC(6,M-1,MB))/TDBC
               ENDIF
             ELSE
               BCX(6) = YLS(2,N)
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
        IF( IBCT(IEQW,NB).EQ.12 ) BCX(3) = PLB(1,NB)
        IF( IBCT(4,NB).EQ.12 ) BCX(4) = PGB(1,NB)
!
!---    Aqueous component intial conditions  ---
!
        IF(IBCT(5,NB).EQ.12) THEN
         DO IGC = 1, NGC
           BCX(6+IGC) = XLCB(IGC,1,NB)
         ENDDO
        ENDIF 
!
!---    Gas component intial conditions  ---
!
        IF(IBCT(6,NB).EQ.12) THEN
         DO IGC = 1, NGC
           BCX(6+LNGC+IGC) = XGCB(IGC,1,NB)
         ENDDO
        ENDIF 
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 )  THEN
          IF( IBCT(IEQT,NB).EQ.12 ) BCX(2) = TB(1,NB)
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 )  THEN
          IF( IBCT(IEQS,NB).EQ.12 ) BCX(6) = YLSB(1,NB)
        ENDIF
        N = IBCN(NB)
        IBD = ABS(IBCD(NB))
        IZN = IZ(N)
#ifdef ecke
        POR0(1,N) = POR0(1,N)
        POR0(2,N) = POR0(2,N)
#endif
!
!---    Assign gas-entry pressure for non Brooks-Corey;
!       Brooks-Corey; Brooks-Corey, Dual Porosity; and
!       Brooks-Corey, Entrapment  ---
!
        IF( ISCHR(IZN).EQ.2 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.102 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
        ELSE
          ENPR = 0.D+0
        ENDIF
!
!---    Initial trapped gas saturation for the van Genuchten or
!       Brooks/Corey entrapment model  ---
!
        IF( ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 ) THEN
          SGRMX = SCHR(15,IZN)
        ELSE
          SGRMX = 0.D+0
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
!---    Loop over secondary variable indices  ---
!
        DO 300 M = 2,ISVC+2
          TX = T(M,N)
          PLX = PL(M,N)
          PGX = PG(M,N)
!
!---      Gas Dirichlet  ---
!
          IF( IBCT(4,NB).EQ.1 .OR. IBCT(4,NB).EQ.12 ) THEN
            PGX = BCX(4)
!
!---      Gas Neumann  ---
!
          ELSEIF( IBCT(4,NB).EQ.2 ) THEN
            PGX = PGX + BCX(4)*DB*VISG(M,N)/PERM(IBD,IZN)
     &        + RHOG(M,N)*GB
!
!---      Gas Zero Flux  ---
!
          ELSEIF( IBCT(4,NB).EQ.3 ) THEN
            IF( ABS(BCX(4)+PATM).GT.EPSL ) THEN
              PGX = BCX(4)
            ELSE
              PGX = PGX + RHOG(M,N)*GB
            ENDIF
!
!---      Gas Hydraulic Gradient Base Node  ---
!
          ELSEIF( IBCT(4,NB).EQ.-11 ) THEN
            PGX = BCX(4)
!
!---      Gas Hydraulic Gradient  ---
!
          ELSEIF( IBCT(4,NB).EQ.11 ) THEN
            IF( M.EQ.2 ) THEN
              PGX = HGBG_COMP( NBHG(4,NB),NB,M )
            ELSE
              PGX = PGB(2,NB)
            ENDIF
          ENDIF
!
!---      Aqueous Dirichlet  ---
!
          IF( IBCT(IEQW,NB).EQ.1 .OR. IBCT(IEQW,NB).EQ.12 ) THEN
            PLX = BCX(3)
!
!---      Aqueous Neumann  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.2 ) THEN
            PLX = PLX + BCX(3)*DB*VISL(M,N)/PERM(IBD,IZN)
     &        + RHOL(M,N)*GB
!
!---      Aqueous Zero Flux  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.3 ) THEN
            IF( ABS(BCX(3)+PATM).GT.EPSL ) THEN
              PLX = BCX(3)
            ELSE
              PLX = PLX + RHOL(M,N)*GB
            ENDIF
!
!---      Aqueous Saturated  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.4 ) THEN
            PLX = PGX
!
!---      Aqueous Hydraulic Gradient Base Node  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.-11 ) THEN
            PLX = BCX(3)
!
!---      Aqueous Hydraulic Gradient  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.11 ) THEN
            IF( M.EQ.2 ) THEN
              PLX = HGBL_COMP( NBHG(1,NB),NB,M )
            ELSE
              PLX = PLB(2,NB)
            ENDIF
!
!---      Aqueous Seepage Face Base Node  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.-11 ) THEN
            PLX = MAX( PGX,BCX(3) )
!
!---      Aqueous Seepage Face  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.11 ) THEN
            IF( M.EQ.2 ) THEN
              PLX = MAX( PGX,HGBL_COMP( NBHG(1,NB),NB,M ) )
            ELSE
              PLX = PLB(2,NB)
            ENDIF
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Dirichlet  ---
!
            IF( IBCT(IEQT,NB).EQ.1 .OR. IBCT(IEQT,NB).EQ.12 ) THEN
              TX = BCX(2)
!
!---        Energy Neumann  ---
!
            ELSEIF( IBCT(IEQT,NB).EQ.2 ) THEN
              INDX = ABS( IBCD(NB) )
              TKP = MAX(1.D+0-PORD(M,N),0.D+0)*THKS(INDX,IZ(N)) + 
     &          PORD(M,N)*(THKL(M,N)*SL(M,N) + 
     &          THKG(M,N)*SG(M,N))
              TX = TX + BCX(2)*DB/TKP
!
!---        Energy zero flux or outflow  ---
!
            ELSEIF( IBCT(IEQT,NB).EQ.3 .OR. IBCT(IEQT,NB).EQ.7 ) THEN
              TX = TX
            ENDIF
          ENDIF
!
!---      Convert pressure to absolute prior to computing physical
!         properties  ---
!
          PLX = PLX + PATM
          PGX = PGX + PATM
          PX = MAX( PLX,PGX )
          PCX = MAX( PGX-PLX,0.D+0 )
!
!---      Aqueous component properties defined as mass fractions,
!         mole fractions, or initial conditions  ---
!
          IF( IBCT(5,NB).EQ.36 .OR. IBCT(5,NB).EQ.37 .OR.
     &      IBCT(5,NB).EQ.43 .OR. IBCT(5,NB).EQ.44 .OR.
     &      IBCT(5,NB).EQ.12 ) THEN
!
!---        Gas component properties defined as zero flux or outflow,
!           use aqueous boundary values to calculate aqueous and
!           gas properties, assuming equilibrium conditions  ---
!
            IF( IBCT(6,NB).EQ.3 .OR. IBCT(6,NB).EQ.7 ) THEN
!
!---          Aqueous component compositions defined as aqueous dry-gas
!             component mass fractions  ---
!
              IF( IBCT(5,NB).EQ.36 .OR. IBCT(5,NB).EQ.37 ) THEN
!
! --- Calculate dry gas mole fraction
!
              XMVA = 0.
              DO IGC = 1, NGC
                XLCB(IGC,M,NB) = BCX(6+IGC)
                XMVA = XMVA + XLCB(IGC,M,NB)/GCPP(1,IGC)
              ENDDO
              IF ( XMVA.EQ.0 ) XMVA = 1.D+0
              DO IGC = 1,NGC
                XMDGCBX(IGC) = XLCB(IGC,M,NB)/GCPP(1,IGC)/XMVA
                XMDLCBX(IGC) = 0.0
              ENDDO
              XMDGWBX = 0.0
              XMDLWBX = 0.0
!               CALL FLH_1( BCX(6+1),
!    &            PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!    &            XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!    &            XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!    &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!    &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB))
!
!---          Aqueous component compositions defined as initial
!             conditions (i.e., aqueous component mass fractions)  ---
!
              ELSEIF( IBCT(5,NB).EQ.12 ) THEN
!               CALL FLH_2( PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!    &            XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!    &            XLAB(M,NB),BCX(6+1),XLWB(M,NB),
!    &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!    &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
!
!---          Aqueous component compositions defined as aqueous dry-gas
!             component mole frations  ---
!
              ELSEIF( IBCT(5,NB).EQ.43 .OR. IBCT(5,NB).EQ.44 ) THEN
!               CALL FLH_3( BCX(6+1),PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!    &            XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!    &            XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!    &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!    &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
              ENDIF
!
!---        Gas component compositions defined as gas dry-gas
!           component mass fractions, use aqueous and gas boundary
!           values independently to calculate aqueous and
!           gas properties, assuming nonequilibrium conditions  ---
!
            ELSEIF( IBCT(6,NB).EQ.36 .OR. IBCT(6,NB).EQ.37 ) THEN
!
!---          Aqueous component compositions defined as aqueous dry-gas
!             component mass fractions  ---
!
              IF( IBCT(5,NB).EQ.36 .OR. IBCT(5,NB).EQ.37 ) THEN
!                CALL FLH_4( BCX(6+LNGC+1),BCX(6+1),
!     &            PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &            XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!     &            XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!     &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
!
!---          Aqueous component compositions defined as initial
!             conditions (i.e., aqueous component mass fractions)  ---
!
              ELSEIF( IBCT(5,NB).EQ.12 ) THEN
!                CALL FLH_5( BCX(6+LNGC+1),
!     &            PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &            XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!     &            XLAB(M,NB),BCX(6+1),XLWB(M,NB),
!     &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
!
!---          Aqueous component compositions defined as aqueous dry-gas
!             component mole frations  ---
!
              ELSEIF( IBCT(5,NB).EQ.43 .OR. IBCT(5,NB).EQ.44 ) THEN
!                CALL FLH_6( BCX(6+LNGC+1),BCX(6+1),
!     &            PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &            XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!     &            XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!     &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
              ENDIF
!
!---        Gas component properties defined as initial conditions,
!           use aqueous and gas boundary values independently to
!           calculate aqueous and gas properties,
!           assuming nonequilibrium conditions  ---
!
            ELSEIF( IBCT(6,NB).EQ.12 ) THEN
!
!---          Aqueous component compositions defined as aqueous dry-gas
!             component mass fractions  ---
!
              IF( IBCT(5,NB).EQ.36 .OR. IBCT(5,NB).EQ.37 ) THEN
!                CALL FLH_7( BCX(6+1),
!     &            PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &            XGAB(M,NB),BCX(6+LNGC+1),XGWB(M,NB),
!     &            XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!     &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
!
!---          Aqueous component compositions defined as initial
!             conditions (i.e., aqueous component mass fractions)  ---
!
              ELSEIF( IBCT(5,NB).EQ.12 ) THEN
!                CALL FLH_8( PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &            XGAB(M,NB),BCX(6+LNGC+1),XGWB(M,NB),
!     &            XLAB(M,NB),BCX(6+1),XLWB(M,NB),
!     &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
!
!---          Aqueous component compositions defined as aqueous dry-gas
!             component mole frations  ---
!
              ELSEIF( IBCT(5,NB).EQ.43 .OR. IBCT(5,NB).EQ.44 ) THEN
!                CALL FLH_9( BCX(6+1),
!     &            PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &            XGAB(M,NB),BCX(6+LNGC+1),XGWB(M,NB),
!     &            XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!     &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
              ENDIF
!
!---        Gas component compositions defined as gas dry-gas
!           component mole fractions, use aqueous and gas boundary
!           values independently to calculate aqueous and
!           gas properties, assuming nonequilibrium conditions  ---
!
            ELSEIF( IBCT(6,NB).EQ.43 .OR. IBCT(6,NB).EQ.44 ) THEN
!
!---          Aqueous component compositions defined as aqueous dry-gas
!             component mass fractions  ---
!
              IF( IBCT(5,NB).EQ.36 .OR. IBCT(5,NB).EQ.37 ) THEN
!                CALL FLH_10( BCX(6+LNGC+1),BCX(6+1),
!     &            PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &            XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!     &            XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!     &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
!
!---          Aqueous component compositions defined as initial
!             conditions (i.e., aqueous component mass fractions)  ---
!
              ELSEIF( IBCT(5,NB).EQ.12 ) THEN
!                CALL FLH_11( BCX(6+LNGC+1),
!     &            PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &            XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!     &            XLAB(M,NB),BCX(6+1),XLWB(M,NB),
!     &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
!
!---          Aqueous component compositions defined as aqueous dry-gas
!             component mole frations  ---
!
              ELSEIF( IBCT(5,NB).EQ.43 .OR. IBCT(5,NB).EQ.44 ) THEN
!                CALL FLH_12( BCX(6+LNGC+1),BCX(6+1),
!     &            PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &            XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!     &            XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!     &            XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &            XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
              ENDIF
            ENDIF
!
!---      Aqueous component properties defined as zero flux,
!         outflow, or relative saturation  ---
!
          ELSEIF( IBCT(5,NB).EQ.3 .OR. IBCT(5,NB).EQ.7 .OR.
     &      IBCT(5,NB).EQ.34 .OR. IBCT(5,NB).EQ.35 ) THEN
!
!---        Gas component compositions defined as gas dry-gas
!           component mass fractions, use gas boundary values
!           to calculate aqueous and gas properties, assuming
!           equilibrium conditions  ---
!
            IF( IBCT(6,NB).EQ.36 .OR. IBCT(6,NB).EQ.37 ) THEN
!              CALL FLH_13( BCX(6+LNGC+1),BCX(6+1),
!     &          PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &          XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!     &          XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!     &          XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &          XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
!
!---        Gas component properties defined as initial conditions,
!           use gas boundary values to calculate aqueous and gas
!           properties, assuming equilibrium conditions  ---
!
            ELSEIF( IBCT(6,NB).EQ.12 ) THEN
!              CALL FLH_14( BCX(6+1),
!     &          PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &          XGAB(M,NB),BCX(6+LNGC+1),XGWB(M,NB),
!     &          XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!     &          XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &          XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
!
!---        Gas component compositions defined as gas dry-gas
!           component mole fractions, use gas boundary values
!           to calculate aqueous and gas properties, assuming
!           equilibrium conditions  ---
!
            ELSEIF( IBCT(6,NB).EQ.43 .OR. IBCT(6,NB).EQ.44 ) THEN
              DO IGC = 1,NGC
                XMDGCBX(IGC)=BCX(6+LNGC+IGC)
                XMDLCBX(IGC) = 0.0
              ENDDO
              XMDGWBX = 0.0
              XMDLWBX = 0.0
!              CALL FLH_15( BCX(6+LNGC+1),BCX(6+1),
!     &          PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &          XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!     &          XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!     &          XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &          XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
!
!---        Gas component compositions defined as zero flux
!           or outflow, assign aqueous properties assuming
!           zero component concentrations and gas properties
!           assuming uniform component concentrations  ---
!
            ELSEIF( IBCT(6,NB).EQ.3 .OR. IBCT(6,NB).EQ.7 ) THEN
              DO IGC = 1,NGC
                XMDGCBX(IGC) = XMVGC(IGC,2,N)
                XMDLCBX(IGC) = 0.0
              ENDDO
              XMDGWBX = 0.0
              XMDLWBX = 0.0
!              CALL FLH_16( PVAB(M,NB),PVCB(1,M,NB),PVWB(M,NB),
!     &          XGAB(M,NB),XGCB(1,M,NB),XGWB(M,NB),
!     &          XLAB(M,NB),XLCB(1,M,NB),XLWB(M,NB),
!     &          XMGAB(M,NB),XMGCB(1,M,NB),XMGWB(M,NB),
!     &          XMLAB(M,NB),XMLCB(1,M,NB),XMLWB(M,NB) )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Dissolved-salt aqueous concentration boundary condition  ---
!
            IF( IBCT(IEQS,NB).EQ.40 .OR. IBCT(IEQS,NB).EQ.41 ) THEN
              RHOLSX = BCX(6)
              ISRX = 1
              CALL DENS_W( TX,PX,RHOLWX,RHOGWX,ISRX )
!
!---          Guess dissolved-salt aqueous mass fraction  ---
!
              YLSX(2) = RHOLSX/(RHOLWX+RHOLSX)
              CALL SOL_LS( TX,XLSMX )
              DYLSX = 1.D-4*XLSMX
              NC = 0
!
!---          Single-variable Newton-Raphson loop (YLS)  ---
!
  200         CONTINUE
              NC = NC + 1
              IF( NC.GT.32 ) THEN
                INDX = 17
                N_DB = N
                CHMSG = 'Unconverged Boundary Condition: Dissolved-Salt'
     &            // ' Aqu. Conc.  @Node: '
                CALL WRMSGS( INDX )
              ENDIF
              DO 210 L = 2,3
                YLSX(L) = YLSX(2)
                IF( L.EQ.3 ) YLSX(L) = YLSX(2) + DYLSX
                XLSX = MIN( YLSX(L),XLSMX )
                CALL SP_B( TX,XLSX,PSBX )
                PX = MAX( PGX,PSBX )
                CALL DENS_B( TX,PX,XLSX,RHOBX )
                IF( ISLC(44).EQ.1 ) THEN
                  PVBX = PSBX
                ELSE
                  CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,XLSX )
                ENDIF
                PVBX = BCX(5)*PVBX
                PVAX = MAX( PGX-PVBX,0.D+0 )
                PGAX = 0.D+0
                DO IGC = 1,NGC
                  PGAX = PGAX + PGCX(IGC)
                ENDDO
                CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &            XGCX,XGWX,XLCX,XLSX,XLWX,XMGCX,XMGWX,
     &            XMLCX,XMLSX,XMLWX,ZGX,WTMGCX,XMDGCBX,XMDGWBX,XMDLCBX,
     &            XMDLWBX )
                DO IGC = 1, NGC
                 XLCX(IGC) = BCX(6+IGC)*XLCX(IGC)
                ENDDO
                CALL DENS_L( TX,RHOBX,XLCX,RHOLX )
                GX(L,1) = XLSX - RHOLSX/RHOLX
  210         CONTINUE
              RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
              RPX(1) = -GX(2,1)
              CYLSX = RPX(1)/RX(1,1)
              YLSX(2) = YLSX(2) + CYLSX
              IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 200
              YLSB(M,NB) = YLSX(2)
              XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
!
!---        Dissolved-salt aqueous relative saturation
!           boundary condition  ---
!
            ELSEIF( IBCT(IEQS,NB).EQ.34 .OR. IBCT(IEQS,NB).EQ.35 ) THEN
              PHILSX = BCX(6)
              CALL SOL_LS( TX,XLSMX )
              YLSB(M,NB) = PHILSX*XLSMX
              XLSB(M,NB) = MIN(YLSB(M,NB),XLSMX)
!
!---        Dissolved-salt aqueous mass fraction boundary condition  ---
!
            ELSEIF( IBCT(IEQS,NB).EQ.36 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
              YLSB(M,NB) = BCX(6)
              CALL SOL_LS( TX,XLSMX )
              XLSB(M,NB) = MIN(YLSB(M,NB),XLSMX)
!
!---        Dissolved-salt aqueous concentration boundary condition  ---
!
            ELSEIF( IBCT(IEQS,NB).EQ.38 .OR. IBCT(IEQS,NB).EQ.39 ) THEN
              SBX = BCX(6)
              CALL PORSTY( N,PGX,PCMP(N),PORDX,PORTX )
              PORDX = MAX( PORDX,EPSL )
              PORTX = MAX( PORTX,PORDX )
              RHOLSX = SBX/PORDX
              ISRX = 1
              CALL DENS_W( TX,PGX,RHOLWX,RHOGWX,ISRX )
!
!---          Guess dissolved-salt aqueous mass fractions  ---
!
              YLSX(2) = RHOLSX/(RHOLWX+RHOLSX)
              CALL SOL_LS( TX,XLSMX )
              DYLSX = 1.D-4*XLSMX
              NC = 0
!
!---          Single-variable Newton-Raphson loop (YLS)  ---
!
  220         CONTINUE
              NC = NC + 1
              IF( NC.GT.32 ) THEN
                INDX = 17
                N_DB = N
                CHMSG = 'Unconverged Initial Conditions: Salt' //
     &            ' Vol. Conc., Gas Pressure, Aqueous Pressure  @Node: '
                CALL WRMSGS( INDX )
              ENDIF
              DO 230 L = 2,3
                YLSX(L) = YLSX(2)
                IF( L.EQ.3 ) YLSX(L) = YLSX(2) + DYLSX
                XLSX = MIN( YLSX(L),XLSMX )
                CALL SP_B( TX,XLSX,PSBX )
                PX = MAX( PGX,PSBX )
                CALL DENS_B( TX,PX,XLSX,RHOBX )
                IF( ISLC(44).EQ.1 ) THEN
                  PVBX = PSBX
                ELSE
                  CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,XLSX )
                ENDIF
                PVBX = BCX(5)*PVBX
                PVAX = MAX( PGX-PVBX,0.D+0 )
                PGAX = 0.D+0
                DO IGC = 1,NGC
                  PGAX = PGAX + PGCX(IGC)
                ENDDO
                CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &            XGCX,XGWX,XLCX,XLSX,XLWX,XMGCX,XMGWX,
     &            XMLCX,XMLSX,XMLWX,ZGX,WTMGCX,XMDGCBX,XMDGWBX,
     &            XMDLCBX,XMDLWBX)
                DO IGC = 1,NGC
                  XLCX(IGC) = BCX(6+IGC)*XLCX(IGC)
                ENDDO
                 CALL DENS_L( TX,RHOBX,XLCX,RHOLX )
                GX(L,1) = XLSX - RHOLSX/RHOLX
  230         CONTINUE
              RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
              RPX(1) = -GX(2,1)
              CYLSX = RPX(1)/RX(1,1)
              YLSX(2) = YLSX(2) + CYLSX
              IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 220
              YLSB(M,NB) = YLSX(2)
              XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
            ENDIF
          ELSE
            YLSB(M,NB) = 0.D+0
            CALL SOL_LS( TX,XLSMX )
            XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
          ENDIF
!
!---      Dissolved-CO2 mass fraction, water-vapor partial pressure,
!         CO2 partial pressure  ---
!
          CALL SP_B( TX,XLSB(M,NB),PSBX )
          PX = MAX( PGX,PSBX )
          CALL DENS_B( TX,PX,XLSB(M,NB),RHOBX )
          IF( ISLC(44).EQ.1 ) THEN
            PVBX = PSBX
          ELSE
            CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,XLSB(M,NB) )
          ENDIF
          PVBX = BCX(5)*PVBX
          PVWB(M,NB) = PVBX
          PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
          IF( PVAB(M,NB).GT.EPSL ) THEN
            CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &        XGCX,XGWX,XLCX,XLSB(M,NB),XLWX,XMGCX,XMGWX,XMLCX,
     &        XMLSX,XMLWX,ZGX,WTMGCX,XMDGCBX,XMDGWBX,
     &        XMDLCBX,XMDLWBX)
            XLAB(M,NB) = 0.D+0
            DO IGC = 1, NGC
              IF( IBCT(5,NB).EQ.34 )
     &        XLCB(IGC,M,NB) = BCX(6+IGC)*XLCX(IGC)
              XLAB(M,NB) = XLAB(M,NB) + XLCB(IGC,M,NB)
            ENDDO
            XLWB(M,NB) = 1.D+0-XLAB(M,NB)-XLSB(M,NB)
            XGAB(M,NB) = 0.D+0
            XMGAB(M,NB) = 0.D+0
            DO IGC = 1, NGC
              XGCB(IGC,M,NB) = XGCX(IGC)
              XGAB(M,NB) = XGAB(M,NB) + XGCB(IGC,M,NB)
              XMGCB(IGC,M,NB) = XMGCX(IGC)
              XMGAB(M,NB) = XMGAB(M,NB) + XMGCB(IGC,M,NB)
            ENDDO
            XGWB(M,NB) = XGWX
            XMGWB(M,NB) = XMGWX
          ELSE
            DO IGC = 1, NGC
             XLCB(IGC,M,NB) = 0.D+0
             XLAB(M,NB) = 0.0
            ENDDO
            XLWB(M,NB) = 1.D+0-XLAB(M,NB)-XLSB(M,NB)
            DO IGC = 1, NGC
             XGAB(M,N) = 0.D+0
            ENDDO
            XGWB(M,N) = 1.D+0
            XMGWB(M,N) = 1.D+0
          ENDIF
!
!---      Porous-media porosity  ---
!
          CALL PORSTY( N,PX,PCMP(N),PORDB(M,NB),PORTB(M,NB) )
          PORDB(M,NB) = MAX( PORDB(M,NB),EPSL )
          PORTB(M,NB) = MAX( PORTB(M,NB),PORDB(M,NB) )
!
!---      Surface tension and saturation  ---
!
          CALL SFT_L( TX,XLSB(M,NB),SFTLX )
          BTGLB(M,NB) = 1.D+0
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &      BTGLB(M,NB) = SCHR(16,IZN)/SFTLX
          INDX = 0
!
!---      No trapped gas on the boundary surface  ---
!
          ASLMINX = 1.D+0
          CALL KSP_COMP( N,IZN,M,PGX,PLX,SGB(M,NB),SLB(M,NB),
     &      RKLB(1,M,NB),RKGB(M,NB),ASLX,ASLMINX,ASGTX,SGRMX,
     &      SLRX,BTGLB(M,NB),INDX )
!
!---      Aqueous and gas tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SLB(M,NB),SGB(M,NB),ZERO,
     &      PORDB(M,NB),TORLB(M,NB),TORGB(M,NB),TORNX )

!
!---      Gas density and component fractions  ---
!
          ISRX = 2
          CALL DENS_W( TX,PVBX,RHOLWX,RHOGWX,ISRX )
          IF( PVAB(M,NB).GT.EPSL ) THEN
            CALL DENS_GC( TX,PVAB(M,NB),XMGCX,RHOGAX,ISRX,ZGX )
                Z1X = ZGX
          ELSE
            RHOGAX = 0.D+0
          ENDIF
          RHOGB(M,NB) = RHOGWX+RHOGAX
          WTMGX = 1.D+0/(XGAB(M,NB)/WTMGCX + XGWB(M,NB)/WTMW)
          RHOMGB(M,NB) = RHOGB(M,NB)/WTMGX
!
!---      Gas viscosity  ---
!
          DO IGC = 1, NGC
           XMGX(IGC) = XMGCB(IGC,M,NB)
          ENDDO
          CALL VISC_GC( TX,PX,ZGX,XMGX,RHOGAX,VISCGX,VISGCX )
          CALL VISC_W( TX,PVBX,RHOGWX,VISGWX )
          DO IGC = 1, NGC
           XMGX(IGC) = XMGCB(IGC,M,NB)
          ENDDO
          CALL VISC_G( VISCGX,VISGWX,XMGX,XMGWB(M,NB),
     &      VISGB(M,NB) )
          VISGB(M,NB)=VISGB(M,NB)+VISGCX
!
!---      Water-vapor diffusion coefficient  ---
!
          CALL DIFC_GW( TX,PX,DFGWB(M,NB) )
!
!---      Aqueous component fractions and density  ---
!
          XLWB(M,NB) = MAX( 1.D+0-XLSB(M,NB)-XLAB(M,NB),0.D+0 )
          XLAB(M,NB) = 0.D+0
          DO IGC = 1, NGC
           XLAB(M,NB) = XLAB(M,NB) + XLCB(IGC,M,NB)
          ENDDO
          WTMLX = 1.D+0/(XLAB(M,NB)/WTMGX + XLWB(M,NB)/WTMW +
     &      XLSB(M,NB)/WTMS)
          DO IGC = 1, NGC
            XMLCB(IGC,M,NB) = XLCB(IGC,M,NB)*WTMLX/GCPP(1,IGC)
            XLCX(IGC) = XLCB(IGC,M,NB)
          ENDDO
          XMLWB(M,NB) = XLWB(M,NB)*WTMLX/WTMW
          XMLSB(M,NB) = XLSB(M,NB)*WTMLX/WTMS
           CALL DENS_L( TX,RHOBX,XLCX,RHOLB(M,NB) )
          RHOMLB(M,NB) = RHOLB(M,NB)/WTMLX
!
!---      Aqueous viscosity  ---
!
          ISRX = 1
          CALL DENS_W( TX,PX,RHOLWX,RHOX,ISRX )
          CALL VISC_W( TX,PX,RHOLWX,VISLWX )
          CALL VISC_B( TX,XLSB(M,NB),VISLWX,VISBX )
          XMLAB(M,NB) = 0.D+0
          DO IGC = 1, NGC
            XMLAB(M,NB) = XMLAB(M,NB) + XMLCB(IGC,M,NB)
            XMLX(IGC) = XMLCB(IGC,M,NB)
          ENDDO
          CALL VISC_L( XMLX,VISBX,VISCGX,VISLB(M,NB) )
!
!---      Gas diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            IF( NSTEP-NRST.EQ.0 ) THEN
              DO IGC = 1,NGC
                DFGCB(IGC,M,NB) = DFGAC
              END DO
            ENDIF
          ELSEIF( ISLC(2).EQ.2 ) THEN
            DO IGC = 1,NGC
              CALL BNDFGCW( TX,PGX,DFGCWX(IGC),IGC )
            END DO
            DO IGC = 1,NGC-1
              DO JGC = IGC+1,NGC
                CALL BNDFGC2( TX,PGX,DFGC2X(IGC,JGC),IGC,JGC )
              END DO
            END DO
            DO IGC = 1,NGC
              IF( (1.D+0-XMGC(IGC,M,N)).LT.EPSL ) THEN
                DFGCNX = 2.D+0*DFGCWX(IGC)
                DFGCDX = DFGCWX(IGC)
                DO JGC = 1,NGC
                  IF( JGC.NE.IGC ) THEN
                    MNGC = MIN(IGC,JGC)
                    MXGC = MAX(IGC,JGC)
                    DFGCNX = DFGCNX*DFGC2X(MNGC,MXGC)
                    DFGCDX = DFGCDX + DFGC2X(MNGC,MXGC)
                  ENDIF
                ENDDO
                DFGCB(IGC,M,NB) = DFGCNX/DFGCDX
              ELSE
                DFGCNX = 1.D+0-XMGCB(IGC,M,NB)
                DFGCDX = XMGWB(M,NB)/DFGCWX(IGC)
                DO JGC = 1,NGC
                  IF( JGC.NE.IGC ) THEN
                    MNGC = MIN(IGC,JGC)
                    MXGC = MAX(IGC,JGC)
                    DFGCDX = DFGCDX + XMGCB(JGC,M,NB)/DFGC2X(MNGC,MXGC)
                  ENDIF
                ENDDO
                DFGCB(IGC,M,NB) = DFGCNX/DFGCDX
              ENDIF
            ENDDO
          ENDIF
!
!---      Dissolved-component and -NaCl diffusion coefficient  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DO IGC = 1, NGC
              DFLCB(IGC,M,NB) = DFLAC
            ENDDO
            DFLSB(M,NB) = DFLSC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            DO IGC = 1, NGC
              CALL BNDFLC( TX,VISLB(M,NB),VISGC(IGC),
     &          DFLCB(IGC,M,NB),IGC )
            ENDDO
            CALL DIFC_LS( TX,XLSB(M,NB),VISLB(M,NB),DFLSB(M,NB) )
          ENDIF
!
!---      Nonisothermal simulation  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Gas enthalpy and internal energy  ---
!
            ISRX = 2
            WTMX = 0.D+0
            DO IGC = 1, NGC
             WTMX = WTMX + XMGCB(IGC,M,NB)*GCPP(1,IGC)
            ENDDO
            PVX = PVAB(M,NB)
            DO IGC = 1, NGC
              XMGX(IGC) = XMGCB(IGC,M,NB)
            ENDDO
            PVX=PX
            CALL ENTH_GC( TX,PVX,XMGX,XMGWB(M,NB),
     &        WTMX,HGAB(M,NB),UEGAX,ISRX )
            CALL ENTH_W( TX,PVWB(M,NB),HLWX,HGWB(M,NB),ISRX )
            UGWX = HGWB(M,NB) - PVWB(M,NB)/RHOGWX
            UEGB(M,NB) = XGWB(M,NB)*UGWX+XGAB(M,NB)*UEGAX
            HGB(M,NB) = XGWB(M,NB)*HGWB(M,NB)+XGAB(M,NB)*HGAB(M,NB)
!
!---        Gas thermal conductivity  ---
!
            DO IGC = 1, NGC
              XMGX(IGC) = XMGCB(IGC,M,NB)
            ENDDO
            CALL THK_GC( PGX,RHOGAX,THKGAX,TX,XMGX,Z1X )
            CALL THK_W( TX,PGX,RHOGWX,THKGWX )
            DO IGC = 1, NGC
              XMGX(IGC) = XMGCB(IGC,M,NB)
            ENDDO
            CALL THK_G( TX,THKGAX,THKGWX,
     &        XMGX,XMGWB(M,NB),THKGB(M,NB) )
!
!---        Aqueous enthalpy and internal energy  ---
!
            ISRX = 1
            CALL ENTH_W( TX,PX,HLWX,HGWX,ISRX )
            CALL ENTH_B( TX,XLSB(M,NB),HLWX,HBX )
!
!---        Change in solubility with temperature ---
!
! values for source or node?
!            DO IGC = 1, NGC
!              DHCX(IGC) = XMLCX(IGC)
!            END DO
            DTX = 1.D-6
            CALL EQUIL( TX+DTX,PX,PGAX,PGWX,PSBX,PVBX,
     &        XGCX,XGWX,XLCX,XLSB(M,NB),XLWX,XMGCX,XMGWX,XMLCX,
     &        XMLSX,XMLWX,ZGX,WTMGCX,XMDGCBX,XMDGWBX,
     &        XMDLCBX,XMDLWBX)
            DO IGC = 1, NGC
!              IF( XMLCX(IGC).GT.EPSL ) THEN
!                DHCX(IGC) = LOG(DHCX(IGC)/XMLCX(IGC))/DTX
!              ELSE
                DHCX(IGC) = 0.D+0
!              ENDIF
            END DO
            CALL ENTH_L( TX,XLSB(M,NB),XLCX,DHCX,
     &        HBX,HGAB(M,NB),HLB(M,NB) )
!  
!---        Aqueous thermal conductivity  ---
!
            CALL THK_W( TX,PX,RHOLWX,THKLWX )
            CALL THK_B( TX,XLSB(M,NB),THKLWX,THKLB(M,NB) )
          ENDIF
!
!---      Assign boundary primary variables  ---
!
          TB(M,NB) = TX
          PLB(M,NB) = PLX - PATM
          PGB(M,NB) = PGX - PATM
          PVWB(M,NB) = PVBX
  300   CONTINUE
  400 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCP_COMP group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CAP_COMP( IZN,SLX,SGTX,BTGLX,CPGL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2 Mode
!
!     Compute the gas/aqueous capillary pressure from the aqueous
!     saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/CAP_COMP'
!
!---  van Genuchten saturation function with the Webb extension  ---
!
      IF( ISCHR(IZN).EQ.1 .AND. ISM(IZN).EQ.2 ) THEN
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
!---    Aqueous saturation below the matching point,
!       use Webb extension  ---
!
        IF( SLX.LT.SMPX ) THEN
          HMPX = SCHR(9,IZN)
          DMPX = -(LOG10(HDOD)-LOG10(HMPX))/SMPX
          HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---    Aqueous saturation at or above the matching point,
!       use van Genuchten function
!
        ELSE
          ASLX = (SLX-SCHR(4,IZN))/(1.-SCHR(4,IZN))
          HDGL = (((1./ASLX)**(1./CM)-1.)**(1./CN))/SCHR(1,IZN)
        ENDIF
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  van Genuchten saturation function w/ or w/o
!     the Fayer-Simmons extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.1 ) THEN
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
        IF( SLX.GT.SCHR(4,IZN ) ) THEN
          ASLX = (SLX-SCHR(4,IZN))/(1.-SCHR(4,IZN))
          HDGL = (((1./ASLX)**(1./CM)-1.)**(1./CN))/SCHR(1,IZN)
        ELSEIF( ISM(IZN).EQ.1 ) THEN
          HDGL = EXP((1.D+0-SLX/SCHR(4,IZN))*LOG(HDOD))
        ELSE
          HDGL = HDOD
          SLX = SCHR(4,IZN) + 1.D-6
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  100   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        DSLRX = -SCHR(4,IZN)/(HDGL*LOG(HDOD))*REALX
        ASLX = 1.D+0/((1.D+0 + (SCHR(1,IZN)*HDGL)**CN)**CM)
        DASLX = -CM*SCHR(1,IZN)*CN*((SCHR(1,IZN)*HDGL)**(CN-1.D+0))
     &  /((1.D+0 + (SCHR(1,IZN)*HDGL)**CN)**(CM+1.D+0))
        SLZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
        F = SLX - SLZ
        DF = -DSLZ
        DH = -F/(DF+SMALL)
        HDGL = HDGL + DH
!
!---    No convergence on van Genuchten 
!       capillary pressure  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 14
          RLMSG = SLX
          CHMSG = 'No Convergence on van Genuchten '
     &      // 'Capillary Pressure: Aqueous Saturation = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 100
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Brooks and Corey saturation function with the Webb extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 .AND. ISM(IZN).EQ.2 ) THEN
        CL = MAX( SCHR(3,IZN),SMALL )
!
!---    Aqueous saturation below the matching point,
!       use Webb extension  ---
!
        IF( SLX.LT.SMPX ) THEN
          HMPX = SCHR(9,IZN)
          DMPX = -(LOG10(HDOD)-LOG10(HMPX))/SMPX
          HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---    Aqueous saturation at or above the matching point,
!       use Brooks and Corey function
!
        ELSE
          ASLX = (SLX-SCHR(4,IZN))/(1.-SCHR(4,IZN))
          HDGL = SCHR(1,IZN)*(1./ASLX)**(1./CL)
        ENDIF
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Brooks and Corey saturation function w/ or w/o
!     the Fayer-Simmons extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 ) THEN
        IF( (1.D+0-SLX)/EPSL.LT.EPSL ) THEN
          HDGL = SCHR(1,IZN)
          GOTO 210
        ENDIF
        CL = MAX( SCHR(3,IZN),SMALL )
        IF( SLX.GT.SCHR(4,IZN ) ) THEN
          ASLX = (SLX-SCHR(4,IZN))/(1.-SCHR(4,IZN))
          HDGL = SCHR(1,IZN)*(1./ASLX)**(1./CL)
        ELSEIF( ISM(IZN).EQ.1 ) THEN
          HDGL = EXP((1.D+0-SLX/SCHR(4,IZN))*LOG(HDOD))
        ELSE
          HDGL = HDOD
          SLX = SCHR(4,IZN) + 1.D-9
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  200   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        DSLRX = -SCHR(4,IZN)*REALX/(HDGL*LOG(HDOD))
        HDGLX = MAX( SCHR(1,IZN),HDGL )
        ASLX = (SCHR(1,IZN)/HDGLX)**CL
        DASLX = -CL*(SCHR(1,IZN)/(HDGLX**2))
     &    *(SCHR(1,IZN)/HDGLX)**(CL-1.D+0)
        SLZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
        F = SLX - SLZ
        DF = -DSLZ
        DH = -F/(DF+SMALL)
        HDGL = HDGL + DH
!
!---    No convergence on Brooks and Corey 
!       capillary pressure  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 14
          RLMSG = SLX
          CHMSG = 'No Convergence on Brooks & Corey '
     &      // 'Capillary Pressure: Aqueous Saturation = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 200
  210   CONTINUE
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Dual porosity van Genuchten saturation function
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
        IF( SLX.GT.SCHR(4,IZN ) ) THEN
          ASLX = (SLX-SCHR(4,IZN))/(1.-SCHR(4,IZN))
          HDGL = (((1./ASLX)**(1./CMM)-1.)**(1./CNM))/SCHR(1,IZN)
        ELSEIF( ISM(IZN).EQ.1 ) THEN
          HDGL = EXP((1.D+0-SLX/SCHR(4,IZN))*LOG(HDOD))
        ELSE
          HDGL = HDOD
          SLX = SCHR(4,IZN) + 1.D-6
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  300   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
!
!---    Matrix saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        DSLRX = -SCHR(4,IZN)/(HDGL*LOG(HDOD))*REALX
        ASLX = 1.D+0/((1.D+0 + (SCHR(1,IZN)*HDGL)**CNM)**CMM)
        DASLX = -CMM*SCHR(1,IZN)*CNM*((SCHR(1,IZN)*HDGL)**(CNM-1.D+0))
     &  /((1.D+0 + (SCHR(1,IZN)*HDGL)**CNM)**(CMM+1.D+0))
        SLMZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLMZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
!
!---    Fracture saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(7,IZN),ZERO )
        DSLRX = -SCHR(7,IZN)/(HDGL*LOG(HDOD))*REALX
        ASLX = 1.D+0/((1.D+0 + (SCHR(5,IZN)*HDGL)**CNF)**CMF)
        DASLX = -CMF*SCHR(5,IZN)*CNF*((SCHR(5,IZN)*HDGL)**(CNF-1.D+0))
     &  /((1.D+0 + (SCHR(5,IZN)*HDGL)**CNF)**(CMF+1.D+0))
        SLFZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLFZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
        F = SLX - SLMZ*PORD_MX - SLFZ*PORD_FX
        DF = -DSLMZ*PORD_MX -DSLFZ*PORD_FX
        DH = -F/(DF+SMALL)
        HDGL = HDGL + DH
!
!---    No convergence on dual porosity van Genuchten 
!       capillary pressure  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 14
          RLMSG = SLX
          CHMSG = 'No Convergence on Dual Porosity van Genuchten '
     &      // 'Capillary Pressure: Aqueous Saturation = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 300
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Dual Porosity Brooks and Corey saturation function  ---
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
        IF( SLX.GT.SCHR(4,IZN ) ) THEN
          ASLX = (SLX-SCHR(4,IZN))/(1.-SCHR(4,IZN))
          HDGL = SCHR(1,IZN)*(1./ASLX)**(1./CLM)
        ELSEIF( ISM(IZN).EQ.1 ) THEN
          HDGL = EXP((1.D+0-SLX/SCHR(4,IZN))*LOG(HDOD))
        ELSE
          HDGL = HDOD
          SLX = SCHR(4,IZN) + 1.D-9
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  400   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
!
!---    Matrix saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        DSLRX = -SCHR(4,IZN)*REALX/(HDGL*LOG(HDOD))
        HDGLX = MAX( SCHR(1,IZN),HDGL )
        ASLX = (SCHR(1,IZN)/HDGLX)**CLM
        DASLX = -CLM*(SCHR(1,IZN)/(HDGLX**2))
     &    *(SCHR(1,IZN)/HDGLX)**(CLM-1.D+0)
        SLMZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLMZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
!
!---    Fracture saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(7,IZN),ZERO )
        DSLRX = -SCHR(7,IZN)*REALX/(HDGL*LOG(HDOD))
        HDGLX = MAX( SCHR(5,IZN),HDGL )
        ASLX = (SCHR(5,IZN)/HDGLX)**CLF
        DASLX = -CLF*(SCHR(5,IZN)/(HDGLX**2))
     &    *(SCHR(5,IZN)/HDGLX)**(CLF-1.D+0)
        SLFZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLFZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
        F = SLX - SLMZ*PORD_MX - SLFZ*PORD_FX
        DF = -DSLMZ*PORD_MX -DSLFZ*PORD_FX
        DH = -F/(DF+SMALL)
        HDGL = HDGL + DH
!
!---    No convergence Brooks-Corey capillary pressure  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 14
          RLMSG = SLX
          CHMSG = 'No Convergence on Brooks & Corey '
     &      // 'Capillary Pressure: Aqueous Saturation = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 400
  410   CONTINUE
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Haverkamp saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.5 ) THEN
        IF( SLX.GT.SCHR(4,IZN ) ) THEN
          ASLX = (SLX-SCHR(4,IZN))/(1.-SCHR(4,IZN))
          HDGL = SCHR(1,IZN) + SCHR(5,IZN)*
     &      ((SCHR(2,IZN)/ASLX)-SCHR(2,IZN))**(1.D+0/SCHR(3,IZN))
        ELSEIF( ISM(IZN).EQ.1 ) THEN
          HDGL = EXP((1.D+0-SLX/SCHR(4,IZN))*LOG(HDOD))
        ELSE
          HDGL = HDOD
          SLX = SCHR(4,IZN) + 1.D-9
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  500   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        DSLRX = -SCHR(4,IZN)*REALX/(HDGL*LOG(HDOD))
        ASLX = SCHR(2,IZN)/(SCHR(2,IZN)+((HDGL-SCHR(1,IZN))/
     &    SCHR(5,IZN))**SCHR(3,IZN))
        DASLX = -(SCHR(2,IZN)*SCHR(3,IZN)*
     &    (((HDGL-SCHR(1,IZN))/SCHR(5,IZN))**(SCHR(3,IZN)-1.D+0))
     &    /SCHR(5,IZN))/((SCHR(2,IZN)+((HDGL-SCHR(1,IZN))/SCHR(5,IZN))
     &    **SCHR(3,IZN))**2)
        SLZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
        F = SLX - SLZ
        DF = -DSLZ
        DH = -F/(DF+SMALL)
        HDGL = HDGL + DH
!
!---    No convergence Haverkamp capillary pressure  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 14
          RLMSG = SLX
          CHMSG = 'No Convergence on Haverkamp '
     &      // 'Capillary Pressure: Aqueous Saturation = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 500
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Linear interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.10 ) THEN
        ITBX = 0
        HDGL = FNTBLX( SLX,ISLTBL(1,IZN),ISLTBL(2,IZN),ITBX )
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Cubic spline interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.11 ) THEN
        HDGL = FSPLNX( SLX,ISLTBL(1,IZN),ISLTBL(2,IZN) )
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  van Genuchten saturation function w/ gas entrapment  ---
!
      ELSEIF( ISCHR(IZN).EQ.101 ) THEN
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
        ASLX = SLX + SGTX
        IF( ASLX.GT.SCHR(4,IZN ) ) THEN
          ASLX = (SLX-SCHR(4,IZN))/(1.-SCHR(4,IZN))
          ASGTX = SGTX/(1.D+0-SCHR(4,IZN))
          ASLX = ASLX + ASGTX
          HDGL = (((1./ASLX)**(1./CM)-1.)**(1./CN))/SCHR(1,IZN)
        ELSEIF( ISM(IZN).EQ.1 ) THEN
          HDGL = EXP((1.D+0-ASLX/SCHR(4,IZN))*LOG(HDOD))
        ELSE
          HDGL = HDOD
          SLX = SCHR(4,IZN) + SGTX + 1.D-6
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  600   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        DSLRX = -SCHR(4,IZN)/(HDGL*LOG(HDOD))*REALX
        ASLX = 1.D+0/((1.D+0 + (SCHR(1,IZN)*HDGL)**CN)**CM)
        DASLX = -CM*SCHR(1,IZN)*CN*((SCHR(1,IZN)*HDGL)**(CN-1.D+0))
     &  /((1.D+0 + (SCHR(1,IZN)*HDGL)**CN)**(CM+1.D+0))
        SLZ = ASLX*(1.D+0-SLRX) - SGTX + SLRX
        DSLZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
        F = SLX - SLZ
        DF = -DSLZ
        DH = -F/(DF+SMALL)
        HDGL = HDGL + DH
!
!---    No convergence van Genuchten w/ gas entrapment 
!       capillary pressure  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 14
          RLMSG = SLX
          CHMSG = 'No Convergence on van Genuchten w/ Gas Entrapment '
     &      // 'Capillary Pressure: Aqueous Saturation = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 600
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Brooks and Corey saturation function w/ gas entrapment  ---
!
      ELSEIF( ISCHR(IZN).EQ.102 ) THEN
        ASLX = SLX + SGTX        
        IF( (1.D+0-ASLX)/EPSL.LT.EPSL ) THEN
          HDGL = SCHR(1,IZN)
          GOTO 710
        ENDIF
        CL = MAX( SCHR(3,IZN),SMALL )
        IF( ASLX.GT.SCHR(4,IZN ) ) THEN
          ASLX = (SLX-SCHR(4,IZN))/(1.-SCHR(4,IZN))
          ASGTX = SGTX/(1.D+0-SCHR(4,IZN))
          ASLX = ASLX + ASGTX
          HDGL = SCHR(1,IZN)*(1./ASLX)**(1./CL)
        ELSEIF( ISM(IZN).EQ.1 ) THEN
          HDGL = EXP((1.D+0-ASLX/SCHR(4,IZN))*LOG(HDOD))
        ELSE
          HDGL = HDOD
          SLX = SCHR(4,IZN) + SGTX + 1.D-9
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  700   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        DSLRX = -SCHR(4,IZN)*REALX/(HDGL*LOG(HDOD))
        ASLX = (SCHR(1,IZN)/HDGL)**CL
        DASLX = -CL*(SCHR(1,IZN)/(HDGL**2))
     &    *(SCHR(1,IZN)/HDGL)**(CL-1.D+0)
        SLZ = ASLX*(1.D+0-SLRX) - SGTX + SLRX
        DSLZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
        F = SLX - SLZ
        DF = -DSLZ
        DH = -F/(DF+SMALL)
        HDGL = HDGL + DH
!
!---    No convergence Brooks and Corey w/ gas entrapment 
!       capillary pressure  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 14
          RLMSG = SLX
          CHMSG = 'No Convergence on Brooks & Corey w/ Gas Entrapment '
     &      // 'Capillary Pressure: Aqueous Saturation = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-8 ) GOTO 700
  710   CONTINUE
        CPGL = HDGL*RHORL*GRAV/BTGLX
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Check the thermodynamic and hydrologic states declared through
!     user inputs.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE JACOB
      USE HYST
      USE GRID
      USE FDVS
      USE FDVP
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
      REAL*8 GX(4,2),RX(2,2),RPX(2),GXX(LNGC,4,2)
      REAL*8 XGCX(LNGC),XLCX(LNGC), XLCMX(LNGC),DXLCX(LNGC),XLCZ(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),XLAX(LSV,LFD),HGCX(LNGC)
      REAL*8 UEGCX(LNGC)
      REAL*8 XMDGCX(LNGC),RHOLCX(LNGC)
      REAL*8 XMDLCX(LNGC)
      CHARACTER(LEN=12) :: STRING
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHK_COMP'
      EPSLX = 1.D+0-(EPSL*1.D+2)
      XMDGWX = 0.D+0
!
!---  Read Component property file  ---
!
      CALL RDPF_GC
!
!---  Set oven-dried head  ---
!
      HDOD = 1.D+8/RHORL/GRAV
!
!---  Check initial temperature, aqueous pressure, gas pressure,
!     and aqueous saturation  ---
!
      INDX = 0
      DO 100 N = 1,NFLD
        IZN = IZ(N)
        N_DB = N
        IF( IXP(N).EQ.0 ) GOTO 100
        IF( T(2,N).GT.374.14D+0 .OR. T(2,N).LT.0.01D+0 ) THEN
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
        IF( PG(2,N).GT.8.D+8-PATM .OR. PG(2,N).LT.6.1125D+2-PATM ) THEN
          INDX = 16
          IMSG = N
          CHMSG = 'Out of Range Initial Gas Pressure(Pa) @ Node'
          RLMSG = PG(2,N)+PATM
          CALL WRMSGS( INDX )
        ENDIF
        IF( SL(2,N).GT.1.D+0 .OR. SL(2,N).LT.0.D+0 ) THEN
          INDX = 16
          IMSG = N
          CHMSG = 'Out of Range Initial Aqueous Saturation @ Node'
          RLMSG = SL(2,N)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Webb saturation and capillary pressure matching points  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          IF( ISCHR(IZN).EQ.3 ) THEN
            IF( SCHR(14,IZN).LE.0.D+0 ) THEN
              CNX = MAX( SCHR(3,IZN),SMALL )
              IF( IRPN(IZN).EQ.2 ) THEN
                SCHR(14,IZN) = 1.D+0 - 2.D+0/CNX
              ELSE
                SCHR(14,IZN) = 1.D+0 - 1.D+0/CNX
              ENDIF
            ENDIF
            IF( SCHR(15,IZN).LE.0.D+0 ) THEN
              CNX = MAX( SCHR(6,IZN),SMALL )
              IF( IRPN(IZN).EQ.2 ) THEN
                SCHR(15,IZN) = 1.D+0 - 2.D+0/CNX
              ELSE
                SCHR(15,IZN) = 1.D+0 - 1.D+0/CNX
              ENDIF
            ENDIF
            CALL WEBB_VG( IZN )
          ELSEIF( ISCHR(IZN).EQ.4 ) THEN
            CALL WEBB_BC( IZN )
          ENDIF
        ENDIF
#ifdef ecke
!
!---    Load reactive transport total and diffusive porosity  ---
!
        POR0(1,N) = POR(1,IZ(N))
        POR0(2,N) = POR(2,IZ(N))
#endif
  100 CONTINUE
      IF( INDX.GT.0 ) STOP
!
!---  Restart simulations  ---
!
      IF( (IEO.EQ.2)
     &  .AND. (ISLC(21).EQ.0 .OR. ISLC(21).EQ.IOM) ) GOTO 3200
!
!---    Establish reference pressure for soil compressibility  ---
!
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 200
        TCMP(N) = T(2,N)
        IZN = IZ(N)
        IF( CMP(3,IZN).GT.PATM ) THEN
          PCMP(N) = CMP(3,IZN)
        ELSE
          PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
        ENDIF
  200 CONTINUE
!
!---  Convert initial dissolved-salt and dissolved-component into
!     aqueous mass fractions, considering initial saturation
!     and surface tension affects  ---
!
      DO 3000 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 3000
        IZN = IZ(N)
        N_DB = N
!
!---    Assign gas-entry pressure for non Brooks-Corey;
!       Brooks-Corey; Brooks-Corey, Dual Porosity; and
!       Brooks-Corey, Entrapment  ---
!
        IF( ISCHR(IZN).EQ.2 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.102 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
        ELSE
          ENPR = 0.D+0
        ENDIF
!
!---    Initial trapped gas saturation for the van Genuchten or
!       Brooks/Corey entrapment model  ---
!
        IF( ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 ) THEN
          SGRMX = SCHR(15,IZN)
        ELSE
          SGRMX = 0.D+0
        ENDIF
!
!---    Component initialization check  ---
!
        XMGCTX = 0.D+0
        XGCTX = 0.D+0
        XMLCTX = 0.D+0
        XLCTX = 0.D+0
        IMGCX = 0
        IGCX = 0
        IMLCX = 0
        ILCX = 0
        DO IGC = 1,NGC
          XMGCTX = XMGCTX + XMGC(IGC,2,N)
          XGCTX = XGCTX + XGC(IGC,2,N)
          XMLCTX = XMLCTX + XMLC(IGC,2,N)
          XLCTX = XLCTX + XLC(IGC,2,N)
          IF( ABS(XMGC(IGC,2,N))/EPSL.GT.EPSL ) IMGCX = IMGCX+1
          IF( ABS(XGC(IGC,2,N))/EPSL.GT.EPSL ) IGCX = IGCX+1
          IF( ABS(XMLC(IGC,2,N))/EPSL.GT.EPSL ) IMLCX = IMLCX+1
          IF( ABS(XLC(IGC,2,N))/EPSL.GT.EPSL ) ILCX = ILCX+1
        END DO
!
!---    Check for multiple specifications of initial component
!       concentrations  ---
!
        IF( IMGCX.GE.1 .AND.
     &    (IGCX.GE.1 .OR. IMLCX.GE.1 .OR. ILCX.GE.1) ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Multiple Specifications of Component ' //
     &      'Concentration at Node: '
          CALL WRMSGS( INDX )
        ELSEIF( IGCX.GE.1 .AND.
     &    (IMGCX.GE.1 .OR. IMLCX.GE.1 .OR. ILCX.GE.1) ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Multiple Specifications of Component ' //
     &      'Concentration at Node: '
          CALL WRMSGS( INDX )
        ELSEIF( IMLCX.GE.1 .AND.
     &    (IGCX.GE.1 .OR. IMGCX.GE.1 .OR. ILCX.GE.1) ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Multiple Specifications of Component ' //
     &      'Concentration at Node: '
          CALL WRMSGS( INDX )
        ELSEIF( ILCX.GE.1 .AND.
     &    (IGCX.GE.1 .OR. IMLCX.GE.1 .OR. IMGCX.GE.1) ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Multiple Specifications of Component ' //
     &      'Concentration at Node: '
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for under specifications of initial component
!       concentrations  ---
!
        IF( IMGCX.GT.0 .AND. IMGCX.LT.(NGC-1) ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Under Specification of Gas Component ' //
     &      'Mole Fractions at Node: '
          CALL WRMSGS( INDX )
        ELSEIF( IGCX.GT.0 .AND. IGCX.LT.(NGC-1) ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Under Specification of Gas Component ' //
     &      'Mass Fractions at Node: '
          CALL WRMSGS( INDX )
        ELSEIF( IMLCX.GT.0 .AND. IMLCX.LT.(NGC-1) ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Under Specification of Aqueous Component ' //
     &      'Mole Fractions at Node: '
          CALL WRMSGS( INDX )
        ELSEIF( ILCX.GT.0 .AND. ILCX.LT.(NGC-1) ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Under Specification of Aqueous Component ' //
     &      'Mass Fractions at Node: '
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for out-of-range initial component
!       concentrations  ---
!
        IF( XMGCTX.GT.1.D+0 .OR. XMGCTX.LT.0.D+0 ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Initial Gas Component Mole Fractions ' //
     &      'at Node: '
          CALL WRMSGS( INDX )
        ELSEIF( XGCTX.GT.1.D+0 .OR. XGCTX.LT.0.D+0 ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Initial Gas Component Mass Fractions ' //
     &      'at Node: '
          CALL WRMSGS( INDX )
        ELSEIF( XMLCTX.GT.1.D+0 .OR. XMLCTX.LT.0.D+0 ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Initial Aqueous Component Mole Fractions ' //
     &      'at Node: '
          CALL WRMSGS( INDX )
        ELSEIF( XLCTX.GT.1.D+0 .OR. XLCTX.LT.0.D+0 ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Initial Aqueous Component Mass Fractions ' //
     &      'at Node: '
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for out-of-range initial component
!       concentrations  ---
!
        IF( IMGCX.EQ.NGC .AND. ABS(1.D+0-XMGCTX).GT.EPSL ) THEN
          INDX = 7
          IMSG = N
          WRITE(STRING,'(1PE12.5)') XMGCTX
          CHMSG = 'Summation of Initial Gas Component ' //
     &      'Mole Fractions,' // STRING //' at Node'
          CALL WRMSGS( INDX )
        ELSEIF( IGCX.EQ.NGC .AND. ABS(1.D+0-XGCTX).GT.EPSL ) THEN
          INDX = 7
          IMSG = N
          WRITE(STRING,'(1PE12.5)') XGCTX
          CHMSG = 'Summation of Initial Gas Component ' //
     &      'Mass Fractions,' // STRING //' at Node'
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Unsaturated conditions  ---
!
        IF( (ISIC.EQ.1 .AND. (1.D+0-SL(2,N)).GT.1.D-6) .OR.
     &      (ISIC.EQ.2 .AND. (1.D+0-SL(2,N)).GT.1.D-6) .OR.
     &      (ISIC.EQ.3 .AND. (PG(2,N)-PL(2,N)-ENPR).GT.1.D-6) ) THEN
!
!---    Assign initial phase condition  ---
!
        NPHAZ(2,N) = 2
!
!---    Gas pressure and aqueous saturation  ---
!
        IF( ISIC.EQ.1 ) THEN
!
!---    Dissolved-salt aqueous concentration  ---
!
          IF( ICBRN(N).EQ.1 ) THEN
           RHOLSX = TMS(2,N)
           PGX = PG(2,N)+PATM
           ISRX = 1
           CALL DENS_W( T(2,N),PGX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           NC = 0
!
!---       Single-variable Newton-Raphson loop (YLS)  ---
!
  410      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Dissolved-Salt' //
     &         ' Aqu. Conc., Gas Pressure, Aqueous Saturation  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 420 M = 2,3
             YLS(M,N) = YLS(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             CALL SFT_L( T(2,N),XLS(M,N),SFTLX )
             BTGL(M,N) = 1.D+0
             IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &         BTGL(M,N) = SCHR(16,IZN)/SFTLX
             CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(M,N),PCX )
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             PX = MAX( PGX,PSBX )
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             IF( ISLC(44).EQ.1 ) THEN
               PVBX = PSBX
             ELSE
               CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(M,N) )
             ENDIF
             IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX
               CHMSG = 'Initial Condition Transition: Vapor ' //
     &           'Pressure > Gas Pressure @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
               PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             XLA(M,N) = 0.D+0
             DO IGC = 1, NGC
               XLC(IGC,M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
               XLA(M,N) = XLA(M,N) + XLC(IGC,M,N)
               XLCX(IGC) = XLC(IGC,M,N)
             ENDDO
              CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
  420      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RPX(1) = -GX(2,1)
           CYLSX = RPX(1)/RX(1,1)
           YLS(2,N) = YLS(2,N) + CYLSX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 410
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(2,N),PCX )
           PL(2,N) = PG(2,N)-PCX
!
!---    Dissolved-salt relative concentration  ---
!
          ELSEIF( ICBRN(N).EQ.2 ) THEN
           PHILSX = TMS(2,N)
           PGX = PG(2,N)+PATM
           CALL SOL_LS( T(2,N),XLSMX )
           YLS(2,N) = PHILSX*XLSMX
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           PX = MAX( PGX,PSBX )
           CALL DENS_B( T(2,N),PX,XLS(2,N),RHOBX )
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(2,N),PCX )
           PL(2,N) = PG(2,N)-PCX
           IF( ISLC(44).EQ.1 ) THEN
             PVBX = PSBX
           ELSE
             CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(2,N) )
           ENDIF
           IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX
             CHMSG = 'Initial Condition Transition: Vapor ' //
     &         'Pressure > Gas Pressure @ Node'
             CALL WRMSGS( INDX )
           ENDIF
           PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
             PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLA(2,N) = 0.D+0
           DO IGC = 1, NGC
            XLC(IGC,2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
            XLA(2,N) = XLA(2,N) + XLC(IGC,2,N)
           ENDDO
!
!---    Dissolved-salt aqueous mass fraction  ---
!
          ELSEIF( ICBRN(N).EQ.3 .OR. ICBRN(N).EQ.0 ) THEN
           YLS(2,N) = TMS(2,N)
           PGX = PG(2,N)+PATM
           PCX = PG(2,N)-PL(2,N)
           CALL SOL_LS( T(2,N),XLSMX )
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           PX = MAX( PGX,PSBX )
           CALL DENS_B( T(2,N),PX,XLS(2,N),RHOBX )
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(2,N),PCX )
           PL(2,N) = PG(2,N)-PCX
           IF( ISLC(44).EQ.1 ) THEN
             PVBX = PSBX
           ELSE
             CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(2,N) )
           ENDIF
           IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX
             CHMSG = 'Initial Condition Transition: Vapor ' //
     &         'Pressure > Gas Pressure @ Node'
             CALL WRMSGS( INDX )
           ENDIF
           PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
             PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLA(2,N) = 0.D+0
           DO IGC = 1, NGC
            XLC(IGC,2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
            XLA(2,N) = XLA(2,N) + XLC(IGC,2,N)
           ENDDO
!
!---    Dissolved-salt volumetric concentration  ---
!
          ELSEIF( ICBRN(N).EQ.4 ) THEN
           PGX = PG(2,N)+PATM
           IF( CMP(3,IZN).GT.PATM ) THEN
             PCMP(N) = CMP(3,IZN)
           ELSE
             PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
           ENDIF
#ifdef ecke
           POR0(1,N) = POR0(1,N)
           POR0(2,N) = POR0(2,N)
#endif
           CALL PORSTY( N,PGX,PCMP(N),PORD(2,N),PORT(2,N) )
           PORD(2,N) = MAX( PORD(2,N),EPSL )
           PORT(2,N) = MAX( PORT(2,N),PORD(2,N) )
           RHOLSX = TMS(2,N)/PORD(2,N)
           ISRX = 1
           CALL DENS_W( T(2,N),PGX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           NC = 0
!
!---       Single-variable Newton-Raphson loop (YLS)  ---
!
  610      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Salt' //
     &         ' Vol. Conc., Gas Pressure, Aqueous Pressure  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 620 M = 2,3
             YLS(M,N) = YLS(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             PX = MAX( PGX,PSBX )
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             CALL SFT_L( T(2,N),XLS(M,N),SFTLX )
             BTGL(M,N) = 1.D+0
             IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &         BTGL(M,N) = SCHR(16,IZN)/SFTLX
             CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(M,N),PCX )
             PL(M,N) = PG(2,N)-PCX
             IF( ISLC(44).EQ.1 ) THEN
               PVBX = PSBX
             ELSE
               CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(M,N) )
             ENDIF
             IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX
               CHMSG = 'Initial Condition Transition: Vapor ' //
     &           'Pressure > Gas Pressure @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
               PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
           XLA(M,N) = 0.D+0
           DO IGC = 1, NGC
             XLC(IGC,M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
             XLA(M,N) = XLA(M,N) + XLC(IGC,M,N)
             XLCX(IGC) = XLC(IGC,M,N)
           ENDDO
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
  620      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RPX(1) = -GX(2,1)
           CYLSX = RPX(1)/RX(1,1)
           YLS(2,N) = YLS(2,N) + CYLSX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 610
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(2,N),PCX )
           PL(2,N) = PG(2,N)-PCX
          ENDIF
!
!---    Aqueous pressure and aqueous saturation  ---
!
        ELSEIF( ISIC.EQ.2 ) THEN

!---    Dissolved-salt aqueous concentration  ---
!
          IF( ICBRN(N).EQ.1 ) THEN
           RHOLSX = TMS(2,N)
           CALL SP_W( T(2,N),PSWX )
           ISRX = 1
           CALL DENS_W( T(2,N),PSWX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           NC = 0
!
!---       Single-variable Newton-Raphson loop (YLS)  ---
!
  710      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Dissolved-Salt' //
     &         ' Aqu. Conc., Aqueous Pressure, Aqu. Saturation  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 720 M = 2,3
             YLS(M,N) = YLS(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             CALL SFT_L( T(2,N),XLS(M,N),SFTLX )
             BTGL(M,N) = 1.D+0
             IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &         BTGL(M,N) = SCHR(16,IZN)/SFTLX
             CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(M,N),PCX )
             PGX = PL(2,N)+PCX+PATM
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             PX = MAX( PGX,PSBX )
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             IF( ISLC(44).EQ.1 ) THEN
               PVBX = PSBX
             ELSE
               CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(M,N) )
             ENDIF
             IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX
               CHMSG = 'Initial Condition Transition: Vapor ' //
     &           'Pressure > Gas Pressure @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
               PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
            XLA(M,N) = 0.D+0
            DO IGC = 1, NGC
             XLC(IGC,M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
             XLA(M,N) = XLA(M,N) + XLC(IGC,M,N)
             XLCX(IGC) = XLC(IGC,M,N)
            ENDDO
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
  720      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RPX(1) = -GX(2,1)
           CYLSX = RPX(1)/RX(1,1)
           YLS(2,N) = YLS(2,N) + CYLSX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 710
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(2,N),PCX )
           PG(2,N) = PL(2,N)+PCX
!
!---    Dissolved-salt relative concentration  ---
!
          ELSEIF( ICBRN(N).EQ.2 ) THEN
           PHILSX = TMS(2,N)
           CALL SOL_LS( T(2,N),XLSMX )
           YLS(2,N) = PHILSX*XLSMX
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(2,N),PCX )
           PG(2,N) = PL(2,N)+PCX
           PGX = PG(2,N)+PATM
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           PX = MAX( PGX,PSBX )
           CALL DENS_B( T(2,N),PX,XLS(2,N),RHOBX )
           IF( ISLC(44).EQ.1 ) THEN
             PVBX = PSBX
           ELSE
             CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(2,N) )
           ENDIF
           IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX
             CHMSG = 'Initial Condition Transition: Vapor ' //
     &         'Pressure > Gas Pressure @ Node'
             CALL WRMSGS( INDX )
           ENDIF
           PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
             PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLA(2,N) = 0.D+0
           DO IGC = 1, NGC
             XLC(IGC,2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
             XLA(2,N) = XLA(2,N) + XLC(IGC,2,N)
           ENDDO
!
!---    Dissolved-salt aqueous mass fraction  ---
!
          ELSEIF( ICBRN(N).EQ.3 .OR. ICBRN(N).EQ.0 ) THEN
           YLS(2,N) = TMS(2,N)
           CALL SOL_LS( T(2,N),XLSMX )
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(2,N),PCX )
           PG(2,N) = PL(2,N)+PCX
           PGX = PG(2,N)+PATM
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           PX = MAX( PGX,PSBX )
           CALL DENS_B( T(2,N),PX,XLS(2,N),RHOBX )
           IF( ISLC(44).EQ.1 ) THEN
             PVBX = PSBX
           ELSE
             CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(2,N) )
           ENDIF
           IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX
             CHMSG = 'Initial Condition Transition: Vapor ' //
     &         'Pressure > Gas Pressure @ Node'
             CALL WRMSGS( INDX )
           ENDIF
           PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
             PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
           XLA(2,N) = 0.D+0
           DO IGC = 1, NGC
             XLC(IGC,2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
             XLA(2,N) = XLA(2,N) + XLC(IGC,2,N)
           ENDDO
        
!
!---    Dissolved-salt volumetric concentration  ---
!
          ELSEIF( ICBRN(N).EQ.4 ) THEN
           BTGL(2,N) = 1.D+0
           CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(2,N),PCX )
           PG(2,N) = PL(2,N)+PCX
           PGX = PG(2,N)+PATM
           IF( CMP(3,IZN).GT.PATM ) THEN
             PCMP(N) = CMP(3,IZN)
           ELSE
             PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
           ENDIF
#ifdef ecke
           POR0(1,N) = POR0(1,N)
           POR0(2,N) = POR0(2,N)
#endif
           CALL PORSTY( N,PGX,PCMP(N),PORD(2,N),PORT(2,N) )
           PORD(2,N) = MAX( PORD(2,N),EPSL )
           PORT(2,N) = MAX( PORT(2,N),PORD(2,N) )
           RHOLSX = TMS(2,N)/PORD(2,N)
           ISRX = 1
           CALL DENS_W( T(2,N),PGX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           NC = 0
!
!---       Single-variable Newton-Raphson loop (YLS)  ---
!
 1010      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Salt' //
     &         ' Vol. Conc., Gas Pressure, Aqueous Pressure  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 1020 M = 2,3
             YLS(M,N) = YLS(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             PX = MAX( PGX,PSBX )
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             CALL SFT_L( T(2,N),XLS(M,N),SFTLX )
             BTGL(M,N) = 1.D+0
             IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &         BTGL(M,N) = SCHR(16,IZN)/SFTLX
             CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(M,N),PCX )
             PG(M,N) = PL(2,N)+PCX
             PGX = PG(M,N)+PATM
             IF( ISLC(44).EQ.1 ) THEN
               PVBX = PSBX
             ELSE
               CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(M,N) )
             ENDIF
             IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX
               CHMSG = 'Initial Condition Transition: Vapor ' //
     &           'Pressure > Gas Pressure @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
               PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             XLA(M,N) = 0.D+0
             DO IGC = 1, NGC
              XLC(IGC,M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
              XLA(M,N) = XLA(M,N) + XLC(IGC,M,N)
              XLCX(IGC) = XLC(IGC,M,N)
             ENDDO
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
 1020      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RPX(1) = -GX(2,1)
           CYLSX = RPX(1)/RX(1,1)
           YLS(2,N) = YLS(2,N) + CYLSX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 1010
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           CALL CAP_COMP( IZN,SL(2,N),SGT(2,N),BTGL(2,N),PCX )
           PL(2,N) = PG(2,N)-PCX
          ENDIF
!
!---    Gas pressure and aqueous pressure  ---
!
        ELSEIF( ISIC.EQ.3 ) THEN
!
!---    Dissolved-salt aqueous concentration  ---
!
          IF( ICBRN(N).EQ.1 ) THEN
           RHOLSX = TMS(2,N)
           PGX = PG(2,N)+PATM
           PCX = PG(2,N)-PL(2,N)
           ISRX = 1
           CALL DENS_W( T(2,N),PGX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           NC = 0
!
!---       Single-variable Newton-Raphson loop (YLS)  ---
!
 1110      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Dissolved-Salt' //
     &         ' Aqu. Conc., Gas Pressure, Aqueous Pressure  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 1120 M = 2,3
             YLS(M,N) = YLS(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             PX = MAX( PGX,PSBX )
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             IF( ISLC(44).EQ.1 ) THEN
               PVBX = PSBX
             ELSE
               CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(M,N) )
             ENDIF
             IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX
               CHMSG = 'Initial Condition Transition: Vapor ' //
     &           'Pressure > Gas Pressure @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
               PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             XLA(M,N) = 0.D+0
             DO IGC = 1, NGC
              XLC(IGC,M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
              XLA(M,N) = XLA(M,N) + XLC(IGC,M,N)
              XLCX(IGC) = XLC(IGC,M,N)
             ENDDO
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
 1120      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RPX(1) = -GX(2,1)
           CYLSX = RPX(1)/RX(1,1)
           YLS(2,N) = YLS(2,N) + CYLSX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 1110
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           M = 2
           INDX = 0
           IF( SGT(2,N).GT.EPSL ) THEN
             SG(2,N) = SGT(2,N)
             INDX = 1
           ENDIF
           CALL KSP_COMP( N,IZN,M,PG(2,N),PL(2,N),SG(2,N),SL(2,N),
     &       RKL(1,2,N),RKG(2,N),ASL(N),ASLMIN(2,N),ASGT(N),SGRMX,
     &       SLRX,BTGL(2,N),INDX )
!
!---    Dissolved-salt relative concentration  ---
!
          ELSEIF( ICBRN(N).EQ.2 ) THEN
           PHILSX = TMS(2,N)
           PGX = PG(2,N)+PATM
           PCX = PG(2,N)-PL(2,N)
           CALL SOL_LS( T(2,N),XLSMX )
           YLS(2,N) = PHILSX*XLSMX
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           PX = MAX( PGX,PSBX )
           CALL DENS_B( T(2,N),PX,XLS(2,N),RHOBX )
           IF( ISLC(44).EQ.1 ) THEN
             PVBX = PSBX
           ELSE
             CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(2,N) )
           ENDIF
           IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX
             CHMSG = 'Initial Condition Transition: Vapor ' //
     &         'Pressure > Gas Pressure @ Node'
             CALL WRMSGS( INDX )
           ENDIF
           PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
             PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLA(2,N) = 0.D+0
           DO IGC = 1, NGC
            XLC(IGC,2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
            XLA(M,N) = XLA(M,N) + XLC(IGC,2,N)
           ENDDO
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           M = 2
           INDX = 0
           IF( SGT(2,N).GT.EPSL ) THEN
             SG(2,N) = SGT(2,N)
             INDX = 1
           ENDIF
           CALL KSP_COMP( N,IZN,M,PG(2,N),PL(2,N),SG(2,N),SL(2,N),
     &       RKL(1,2,N),RKG(2,N),ASL(N),ASLMIN(2,N),ASGT(N),SGRMX,
     &       SLRX,BTGL(2,N),INDX )
!
!---    Dissolved-salt aqueous mass fraction  ---
!
          ELSEIF( ICBRN(N).EQ.3 .OR. ICBRN(N).EQ.0 ) THEN
           YLS(2,N) = TMS(2,N)
           PGX = PG(2,N)+PATM
           PCX = PG(2,N)-PL(2,N)
           CALL SOL_LS( T(2,N),XLSMX )
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           PX = MAX( PGX,PSBX )
           CALL DENS_B( T(2,N),PX,XLS(2,N),RHOBX )
           IF( ISLC(44).EQ.1 ) THEN
             PVBX = PSBX
           ELSE
             CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(2,N) )
           ENDIF
           IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX
             CHMSG = 'Initial Condition Transition: Vapor ' //
     &         'Pressure > Gas Pressure @ Node'
             CALL WRMSGS( INDX )
           ENDIF
           PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine here?
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
             PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLA(2,N) = 0.D+0
           DO IGC = 1, NGC
             XLC(IGC,2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC) 
             XLA(2,N) =  XLA(2,N) + XLC(IGC,2,N)
           ENDDO
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           M = 2
           INDX = 0
           IF( SGT(2,N).GT.EPSL ) THEN
             SG(2,N) = SGT(2,N)
             INDX = 1
           ENDIF
           CALL KSP_COMP( N,IZN,M,PG(2,N),PL(2,N),SG(2,N),SL(2,N),
     &       RKL(1,2,N),RKG(2,N),ASL(N),ASLMIN(2,N),ASGT(N),SGRMX,
     &       SLRX,BTGL(2,N),INDX )
!
!---    Salt volumetric concentration  ---
!
          ELSEIF( ICBRN(N).EQ.4 ) THEN
           PGX = PG(2,N)+PATM
           PCX = PG(2,N)-PL(2,N)
           IF( CMP(3,IZN).GT.PATM ) THEN
             PCMP(N) = CMP(3,IZN)
           ELSE
             PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
           ENDIF
#ifdef ecke
           POR0(1,N) = POR0(1,N)
           POR0(2,N) = POR0(2,N)
#endif
           CALL PORSTY( N,PGX,PCMP(N),PORD(2,N),PORT(2,N) )
           PORD(2,N) = MAX( PORD(2,N),EPSL )
           PORT(2,N) = MAX( PORT(2,N),PORD(2,N) )
           RHOLSX = TMS(2,N)/PORD(2,N)
           ISRX = 1
           CALL DENS_W( T(2,N),PGX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           NC = 0
!
!---       Single-variable Newton-Raphson loop (YLS)  ---
!
 1410      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Salt' //
     &         ' Vol. Conc., Gas Pressure, Aqueous Pressure  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 1420 M = 2,3
             YLS(M,N) = YLS(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             PX = MAX( PGX,PSBX )
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             IF( ISLC(44).EQ.1 ) THEN
               PVBX = PSBX
             ELSE
               CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(M,N) )
             ENDIF
             IF( PVBX .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX
               CHMSG = 'Initial Condition Transition: Vapor ' //
     &           'Pressure > Gas Pressure @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
               PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N) 
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             XLA(M,N) = 0.D+0  
             DO IGC = 1, NGC
               XLC(IGC,M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
               XLA(M,N) = XLA(M,N) + XLC(IGC,M,N)  
               XLCX(IGC) = XLC(IGC,M,N)
             ENDDO
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
 1420      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RPX(1) = -GX(2,1)
           CYLSX = RPX(1)/RX(1,1)
           YLS(2,N) = YLS(2,N) + CYLSX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 1410
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
           BTGL(2,N) = 1.D+0
           IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &       BTGL(2,N) = SCHR(16,IZN)/SFTLX
           M = 2
           INDX = 0
           IF( SGT(2,N).GT.EPSL ) THEN
             SG(2,N) = SGT(2,N)
             INDX = 1
           ENDIF
           CALL KSP_COMP( N,IZN,M,PG(2,N),PL(2,N),SG(2,N),SL(2,N),
     &       RKL(1,2,N),RKG(2,N),ASL(N),ASLMIN(2,N),ASGT(N),SGRMX,
     &       SLRX,BTGL(2,N),INDX )
          ENDIF
        ENDIF
!
!---  Saturated conditions  ---
!
        ELSE
!
!---    Assign initial phase condition  ---
!
        NPHAZ(2,N) = 1
!
!---    Dissolved-salt aqueous concentration  ---
!
        IF( ICBRN(N).EQ.1 ) THEN
!
!---      Dissolved-co2 aqueous concentration  ---
!
          IF( ICCOMP(1,N).EQ.1 ) THEN
           RHOLSX = TMS(2,N)
           RHOLAX = PVA(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PLX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           PX = MAX( PGX,PLX )
           ISRX = 1
           CALL DENS_W( T(2,N),PX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt and -co2 aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX+RHOLAX)
           XLA(2,N) = RHOLAX/(RHOLWX+RHOLSX+RHOLAX)
           XLSX = 0.D+0
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           XLAMX = 0.D+0
           DO IGC = 1, NGC
              XLAMX = XLAMX + XLCMX(IGC)
           ENDDO
           DXLAX = 1.D-4*XLAMX
           NC = 0
!
!---       Two-variable Newton-Raphson loop (YLS, XLA)  ---
!
1510      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Dissolved-Salt' //
     &         ' Aqu. Conc. and Dissolved-CO2 Aqu. Conc.  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 1520 M = 2,4
             YLS(M,N) = YLS(2,N)
             XLA(M,N) = XLA(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             IF( M.EQ.4 ) XLA(M,N) = XLA(2,N) + DXLAX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             XLAX(M,N) = 0.D+0
             DO IGC = 1, NGC
              XLAX(M,N) = XLAX(M,N) + XLCX(IGC)
              XLC(IGC,M,N) = XLCX(IGC)
             ENDDO
             PVA(M,N) = 0.D+0
             DO IGC = 1, NGC
               PVC(IGC,M,N) = XLC(IGC,M,N)*PGAX/XLCX(IGC)
               PVA(M,N) = PVA(M,N) + PVC(IGC,M,N)
             ENDDO
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             IF( PSBX+PVA(2,N) .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX+PVA(2,N)
               CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &           'Conditions @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
             GX(M,2) = XLA(M,N) - RHOLAX/RHOLX
 1520      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RX(1,2) = (GX(4,1)-GX(2,1))/DXLAX
           RX(2,1) = (GX(3,2)-GX(2,2))/DYLSX
           RX(2,2) = (GX(4,2)-GX(2,2))/DXLAX
           RPX(1) = -GX(2,1)
           RPX(2) = -GX(2,2)
           CYLSX = (RPX(2)-RPX(1)*RX(2,2)/(RX(1,2)+SMALL))/
     &       (RX(2,1)-RX(1,1)*RX(2,2)/(RX(1,2)+SMALL))
           CXLAX = (RPX(2)-RPX(1)*RX(2,1)/(RX(1,1)+SMALL))/
     &       (RX(2,2)-RX(1,2)*RX(2,1)/(RX(1,1)+SMALL))
           YLS(2,N) = YLS(2,N) + CYLSX
           XLA(2,N) = XLA(2,N) + CXLAX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) .OR.
     &       ABS(CXLAX).GT.(1.D-6*XLAMX) ) GOTO 1510
!
!---      Dissolved-co2 relative saturation  ---
!
          ELSEIF( ICCOMP(1,N).EQ.2 ) THEN
           RHOLSX = TMS(2,N)
           PHILAX = PVA(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           ISRX = 1
           CALL DENS_W( T(2,N),PX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved salt aqueous mass fraction  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX+RHOLAX)
           XLSX = 0.D+0
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           XLAMX = 0.D+0
           DO IGC = 1, NGC
              XLAMX = XLAMX + XLCMX(IGC)
           ENDDO
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           DXLAX = 1.D-4*XLAMX
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLA(2,N) = 0.D+0
           XLAX(2,N) = 0.D+0
           DO IGC = 1, NGC
             XLC(IGC,2,N) = PHILAX*XLCX(IGC)
             XLA(2,N) = XLA(2,N) + XLC(IGC,2,N)
             XLAX(2,N) = XLAX(2,N) + XLCX(IGC)
            ENDDO
           PVA(2,N) = 0.D+0
           DO IGC = 1, NGC
            PVC(IGC,2,N) = MIN( XLC(IGC,2,N)/XLCX(IGC),1.D+0 )*PGAX
            PVA(2,N) = PVA(2,N) + PVC(IGC,2,N)
           ENDDO
           IF( PSBX+PVA(2,N) .GT. PGX ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX+PVA(2,N)
             CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &         'Conditions @ Node'
             CALL WRMSGS( INDX )
           ENDIF
           NC = 0
!
!---       Two-variable Newton-Raphson loop (YLS, XLA)  ---
!
 1610      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Dissolved-Salt' //
     &         ' Aqu. Conc. and Dissolved-CO2 Rel. Sat.  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 1620 M = 2,4
             YLS(M,N) = YLS(2,N)
             XLA(M,N) = XLA(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             IF( M.EQ.4 ) XLA(M,N) = XLA(2,N) + DXLAX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             PVBX = PSBX
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCZ,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
                XLAZ = 0.D+0
                DO IGC = 1, NGC
                 XLAZ = XLAZ + XLCZ(IGC)
                ENDDO
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             PVA(M,N) = MIN( XLA(M,N)/XLAX(M,N),1.D+0 )*PGAX
             DO IGC = 1, NGC
               XLC(IGC,M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
             ENDDO
             IF( PSBX+PVA(2,N) .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX+PVA(2,N)
               CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &           'Conditions @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
             GX(M,2) = XLA(M,N) - PHILAX*XLAX(M,N)
 1620      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RX(1,2) = (GX(4,1)-GX(2,1))/DXLAX
           RX(2,1) = (GX(3,2)-GX(2,2))/DYLSX
           RX(2,2) = (GX(4,2)-GX(2,2))/DXLAX
           RPX(1) = -GX(2,1)
           RPX(2) = -GX(2,2)
           CYLSX = (RPX(2)-RPX(1)*RX(2,2)/(RX(1,2)+SMALL))/
     &       (RX(2,1)-RX(1,1)*RX(2,2)/(RX(1,2)+SMALL))
           CXLAX = (RPX(2)-RPX(1)*RX(2,1)/(RX(1,1)+SMALL))/
     &       (RX(2,2)-RX(1,2)*RX(2,1)/(RX(1,1)+SMALL))
           YLS(2,N) = YLS(2,N) + CYLSX
           XLA(2,N) = XLA(2,N) + CXLAX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) .OR.
     &       ABS(CXLAX).GT.(1.D-6*XLAMX) ) GOTO 1610
!
!---      Dissolved-co2 aqueous mass fraction  ---
!
          ELSEIF( ICCOMP(1,N).EQ.3 ) THEN
           RHOLSX = TMS(2,N)
           XLA(2,N) = PVA(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PLX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           PX = MAX( PGX,PLX )
           ISRX = 1
           CALL DENS_W( T(2,N),PX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           NC = 0
!
!---       Single-variable Newton-Raphson loop (YLS)  ---
!
 1710      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Dissolved-Salt' //
     &         ' Aqu. Conc. and Dissolved-CO2 Aqu. Mass Frac.  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 1720 M = 2,3
             YLS(M,N) = YLS(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             PVBX = PSBX
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             PVA(2,N) = MIN( XLA(2,N)/XLAX(2,N),1.D+0 )*PGAX
             DO IGC = 1, NGC
               XLC(IGC,M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
             ENDDO
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             IF( PSBX+PVA(2,N) .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX+PVA(2,N)
               CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &           'Conditions @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
 1720      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RPX(1) = -GX(2,1)
           CYLSX = RPX(1)/RX(1,1)
           YLS(2,N) = YLS(2,N) + CYLSX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 1710
!
!---      Dissolved-co2 equilibrium partial pressure  ---
!
          ELSEIF( ICCOMP(1,N).EQ.4 ) THEN
           RHOLSX = TMS(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           ISRX = 1
           CALL DENS_W( T(2,N),PX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
           CALL SOL_LS( T(2,N),XLSMX )
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           NC = 0
!
!---       Single-variable Newton-Raphson loop (YLS)  ---
!
 1760      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Dissolved-Salt' //
     &         ' Aqu. Conc. and Dissolved-CO2 Aqu. Mass Frac.  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 1770 M = 2,3
             YLS(M,N) = YLS(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             PVBX = PSBX
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             XLA(M,N) = 0.D+0
             DO IGC = 1, NGC
              XLC(IGC,M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLCX(IGC)
              XLA(M,N) = XLA(M,N) + XLC(IGC,M,N)
              XLCX(IGC) = XLC(IGC,M,N)
             ENDDO
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             IF( PSBX+PVA(2,N) .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX+PVA(2,N)
               CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &           'Conditions @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
 1770      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RPX(1) = -GX(2,1)
           CYLSX = RPX(1)/RX(1,1)
           YLS(2,N) = YLS(2,N) + CYLSX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 1760
          ENDIF
!
!---    Dissolved-salt relative saturation  ---
!
        ELSEIF( ICBRN(N).EQ.2 ) THEN
!
!---      Dissolved-co2 aqueous concentration  ---
!
          IF( ICCOMP(1,N).EQ.1 ) THEN
           PHILSX = TMS(2,N)
           RHOLAX = PVA(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PLX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           PX = MAX( PGX,PLX )
           CALL SOL_LS( T(2,N),XLSMX )
           YLS(2,N) = PHILSX*XLSMX
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL DENS_B( T(2,N),PX,XLS(2,N),RHOBX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
!
!---       Guess dissolved co2 aqueous mass fraction  ---
!
           XLA(2,N) = 0.D+0
           DO IGC = 1, NGC
            XLC(IGC,2,N) = PVC(IGC,2,N)/(RHOBX+PVA(2,N))
            XLA(2,N) = XLA(2,N) + XLC(IGC,2,N)
           ENDDO
           XLSX = 0.D+0
           CALL SP_B( T(2,N),XLSX,PSBX )
           PVBX = PSBX
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           DO IGC = 1, IGC
            DXLCX(IGC) = 1.D-4*XLCMX(IGC)
           ENDDO
           NC = 0
!
!---       Single-variable Newton-Raphson loop (XLA)  ---
!
 1810      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Dissolved-Salt' //
     &         ' Rel. Sat. and Dissolved-CO2 Aqu. Conc.  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 1820 M = 2,3
             DO IGC = 1, NGC
               XLC(IGC,M,N) = XLC(IGC,2,N)
             ENDDO
             IF( M.EQ.3 ) THEN
               XLA(M,N) = 0.D+0
               DO IGC = 1, NGC
                 XLC(IGC,M,N) = XLC(IGC,2,N) + DXLCX(IGC)
                 XLA(M,N) = XLA(M,N) + XLC(IGC,M,N)
               ENDDO
             ENDIF
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
                XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             PVA(M,N) = 0.D+0
             DO IGC = 1, NGC
              PVC(IGC,M,N) = MIN( XLC(IGC,M,N)/XLCX(IGC),
     &                1.D+0 )*PGAX
              PVA(M,N) = PVA(M,N) + PVC(IGC,M,N)
             ENDDO
             IF( PSBX+PVA(2,N) .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX+PVA(2,N)
               CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &           'Conditions @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLA(M,N) - RHOLAX/RHOLX
 1820      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DXLAX
           RPX(1) = -GX(2,1)
           CXLAX = RPX(1)/RX(1,1)
           XLA(2,N) = XLA(2,N) + CXLAX
           IF( ABS(CXLAX).GT.(1.D-6*XLAMX) ) GOTO 1810
!
!---      Dissolved-co2 relative saturation  ---
!
          ELSEIF( ICCOMP(1,N).EQ.2 ) THEN
           PHILSX = TMS(2,N)
           PHILAX = PVA(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           CALL SOL_LS( T(2,N),XLSMX )
           YLS(2,N) = PHILSX*XLSMX
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           PVBX = PSBX
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLA(2,N) = 0.D+0
           DO IGC = 1, NGC
             XLC(IGC,2,N) = PVC(IGC,2,N)*XLCX(IGC)
             XLA(2,N) = XLA(2,N) + XLC(IGC,2,N)
           ENDDO
           PVA(2,N) = 0.D+0
           DO IGC = 1, NGC
             PVC(IGC,2,N) = MIN( XLC(IGC,2,N)/XLCX(IGC),1.D+0 )*PGAX
             PVA(2,N) = PVA(2,N) + PVC(IGC,2,N)
           ENDDO
           IF( PSBX+PVA(2,N) .GT. PGX ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX+PVA(2,N)
             CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &         'Conditions @ Node'
             CALL WRMSGS( INDX )
           ENDIF
!
!---      Dissolved-co2 aqueous mass fraction  ---
!
          ELSEIF( ICCOMP(1,N).EQ.3 ) THEN
           PHILSX = TMS(2,N)
           XLA(2,N) = PVA(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           CALL SOL_LS( T(2,N),XLSMX )
           YLS(2,N) = PHILSX*XLSMX
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           PVBX = PSBX
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLAX(2,N) = 0.D+0
           DO IGC = 1, NGC
            XLAX(2,N) = XLAX(2,N) + XLCX(IGC)
           ENDDO
           PVA(2,N) = MIN( XLA(2,N)/XLAX(2,N),1.D+0 )*PGAX
           IF( PSBX+PVA(2,N) .GT. PGX ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX+PVA(2,N)
             CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &         'Conditions @ Node'
             CALL WRMSGS( INDX )
           ENDIF
!
!---      Dissolved-co2 partial pressure  ---
!
          ELSEIF( ICCOMP(1,N).EQ.4 ) THEN
           PHILSX = TMS(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           CALL SOL_LS( T(2,N),XLSMX )
           YLS(2,N) = PHILSX*XLSMX
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           PVBX = PSBX
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLAX(M,N) = 0.D+0
           DO IGC = 1, NGC
            XLAX(M,N) = XLAX(M,N) + XLCX(IGC)
           ENDDO
           XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX(M,N)
           IF( PSBX+PVA(2,N) .GT. PGX ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX+PVA(2,N)
             CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &         'Conditions @ Node'
             CALL WRMSGS( INDX )
           ENDIF
          ENDIF
!
!---    Dissolved-salt aqueous mass fraction  ---
!
        ELSEIF( ICBRN(N).EQ.3 .OR. ICBRN(N).EQ.0 ) THEN
!
!---      Dissolved-co2 aqueous concentration  ---
!
          IF( ICCOMP(1,N).EQ.1 ) THEN
           YLS(2,N) = TMS(2,N)
           RHOLAX = PVA(2,N)
           DO IGC = 1, NGC
             RHOLCX(IGC) = PVC(IGC,2,N)
           ENDDO
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PLX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           PX = MAX( PGX,PLX )
           CALL SOL_LS( T(2,N),XLSMX )
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL DENS_B( T(2,N),PX,XLS(2,N),RHOBX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
!
!---       Guess dissolved co2 aqueous mass fraction  ---
!
           XLA(2,N) = RHOLAX/(RHOBX+RHOLAX)
           DO IGC = 1, NGC
             XLC(IGC,2,N) = RHOLCX(IGC)/(RHOBX+RHOLAX)
           ENDDO
           XLSX = 0.D+0
           CALL SP_B( T(2,N),XLSX,PSBX )
           PVBX = PSBX
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           DXLAX = 0.D+0
           XLAMX = 0.D+0
           DO IGC = 1, NGC
            DXLCX(IGC) = 1.D-4*XLCMX(IGC)
            DXLAX = DXLAX + DXLCX(IGC)
            XLAMX = XLAMX + XLCMX(IGC)
           ENDDO
           NC = 0
!
!---       Single-variable Newton-Raphson loop (XLA)  ---
!
 2110      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Dissolved-Salt' //
     &         ' Aqu. Mass Frac. and Dissolved-CO2 Aqu. Conc.  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 2120 M = 2,3
             DO IGC = 1, NGC
               XLC(IGC,M,N) = XLC(IGC,2,N)
             IF( M.EQ.3 ) XLC(IGC,M,N) = XLC(IGC,2,N) + DXLCX(IGC)
             ENDDO
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
              XLAX(M,N) = 0.D+0
             DO IGC = 1, NGC
               XLAX(M,N) = XLAX(M,N) + XLCX(IGC)
             ENDDO
             PVA(M,N) = MIN( XLA(M,N)/XLAX(M,N),1.D+0 )*PGAX
             IF( PSBX+PVA(2,N) .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX+PVA(2,N)
               CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &           'Conditions @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             DO IGC = 1, NGC
              GXX(IGC,M,1) = XLC(IGC,M,N) - RHOLCX(IGC)/RHOLX
             ENDDO
 2120      CONTINUE
           DO IGC = 1, NGC
             RX(1,1) = (GXX(IGC,3,1)-GXX(IGC,2,1))/DXLCX(IGC)
             RPX(1) = -GXX(IGC,2,1)
             CXLAX = RPX(1)/RX(1,1)
             XLC(IGC,2,N) = XLC(IGC,2,N) + CXLAX
             IF( ABS(CXLAX).GT.(1.D-6*XLCMX(IGC)) ) GOTO 2110
            ENDDO
!
!---      Dissolved-co2 relative saturation  ---
!
          ELSEIF( ICCOMP(1,N).EQ.2 ) THEN
           YLS(2,N) = TMS(2,N)
           PHILAX = PVA(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           CALL SOL_LS( T(2,N),XLSMX )
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           PVBX = PSBX
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLAX(2,N) = 0.D+0
           DO IGC = 1, NGC
            XLAX(2,N) = XLAX(2,N) + XLCX(IGC)
           ENDDO
           XLA(2,N) = PHILAX*XLAX(2,N)
           PVA(2,N) = MIN( XLA(2,N)/XLAX(2,N),1.D+0 )*PGAX
           IF( PSBX+PVA(2,N) .GT. PGX ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX+PVA(2,N)
             CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &         'Conditions @ Node'
             CALL WRMSGS( INDX )
           ENDIF
!
!---      Dissolved-co2 aqueous mass fraction  ---
!
          ELSEIF( ICCOMP(1,N).EQ.3 ) THEN
           YLS(2,N) = TMS(2,N)
           XLA(2,N) = 0.D+0
           DO IGC = 1,NGC
             XLA(2,N) = XLA(2,N) + XLC(IGC,2,N)
           ENDDO
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PLX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           PX = MAX( PGX,PLX )
           CALL SOL_LS( T(2,N),XLSMX )
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           PVBX = PSBX
!
! --- Calculate dry gas mole fraction
!
           XMVA = 0.
           DO IGC = 1, NGC
             XMVA = XMVA + XLC(IGC,2,N)/GCPP(1,IGC)
           ENDDO
           IF ( XMVA.EQ.0 ) XMVA = 1.D+0
           DO IGC = 1, NGC
             XMVGC(IGC,2,N) = XLC(IGC,2,N)/GCPP(1,IGC)/XMVA
           END DO
!
! --- Call equilibrium routine here!
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
! --- Check
!
           XLAX(2,N) = 0.D+0
           DO IGC = 1, NGC
             XLAX(2,N) = XLAX(2,N) + XLCX(IGC)
           ENDDO
           PVA(2,N) = MIN( XLA(2,N)/XLAX(2,N),1.D+0 )*PGAX
           IF( PSBX+PVA(2,N) .GT. PGX ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX+PVA(2,N)
             CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &         'Conditions @ Node'
             CALL WRMSGS( INDX )
           ENDIF
!
!---      Dissolved-co2 partial pressure  ---
!
          ELSEIF( ICCOMP(1,N).EQ.4 ) THEN
           YLS(2,N) = TMS(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           CALL SOL_LS( T(2,N),XLSMX )
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
           IF( PSBX+PVA(2,N) .GT. PGX ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX+PVA(2,N)
             CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &         'Conditions @ Node'
             CALL WRMSGS( INDX )
           ENDIF
          ENDIF
!
!---    Dissolved-salt volumetric concentration  ---
!
        ELSEIF( ICBRN(N).EQ.4 ) THEN
!
!---      Dissolved-co2 aqueous concentration  ---
!
          IF( ICCOMP(1,N).EQ.1 ) THEN
           RHOLAX = PVA(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PLX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           PX = MAX( PGX,PLX )
           IF( CMP(3,IZN).GT.PATM ) THEN
             PCMP(N) = CMP(3,IZN)
           ELSE
             PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
           ENDIF
#ifdef ecke
           POR0(1,N) = POR0(1,N)
           POR0(2,N) = POR0(2,N)
#endif
           CALL PORSTY( N,PX,PCMP(N),PORD(2,N),PORT(2,N) )
           PORD(2,N) = MAX( PORD(2,N),EPSL )
           PORT(2,N) = MAX( PORT(2,N),PORD(2,N) )
           RHOLSX = TMS(2,N)/PORD(2,N)
           ISRX = 1
           CALL DENS_W( T(2,N),PX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt and -co2 aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX+RHOLAX)
           XLA(2,N) = RHOLAX/(RHOLWX+RHOLSX+RHOLAX)
           XLSX = 0.D+0
           CALL SP_B( T(2,N),XLSX,PSBX )
           PVBX = PSBX
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLAMX = 0.D+0
           DO IGC = 1, NGC
            XLAMX = XLAMX + XLCMX(IGC)
           ENDDO
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           DXLAX = 1.D-4*XLAMX
           NC = 0
!
!---       Two-variable Newton-Raphson loop (YLS, XLA)  ---
!
 2410      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Salt' //
     &         ' Vol. Conc. and Dissolved-CO2 Aqu. Conc.  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 2420 M = 2,4
             YLS(M,N) = YLS(2,N)
             XLA(M,N) = XLA(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             IF( M.EQ.4 ) XLA(M,N) = XLA(2,N) + DXLAX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             XLAX(M,N) = 0.D+0
             DO IGC = 1, NGC
               XLAX(M,N) = XLAX(M,N) + XLCX(IGC)
             ENDDO
             PVA(M,N) = MIN( XLA(M,N)/XLAX(M,N),1.D+0 )*PGAX
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             IF( PSBX+PVA(2,N) .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX+PVA(2,N)
               CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &           'Conditions @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
             GX(M,2) = XLA(M,N) - RHOLAX/RHOLX
 2420      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RX(1,2) = (GX(4,1)-GX(2,1))/DXLAX
           RX(2,1) = (GX(3,2)-GX(2,2))/DYLSX
           RX(2,2) = (GX(4,2)-GX(2,2))/DXLAX
           RPX(1) = -GX(2,1)
           RPX(2) = -GX(2,2)
           CYLSX = (RPX(2)-RPX(1)*RX(2,2)/(RX(1,2)+SMALL))/
     &       (RX(2,1)-RX(1,1)*RX(2,2)/(RX(1,2)+SMALL))
           CXLAX = (RPX(2)-RPX(1)*RX(2,1)/(RX(1,1)+SMALL))/
     &       (RX(2,2)-RX(1,2)*RX(2,1)/(RX(1,1)+SMALL))
           YLS(2,N) = YLS(2,N) + CYLSX
           XLA(2,N) = XLA(2,N) + CXLAX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) .OR.
     &       ABS(CXLAX).GT.(1.D-6*XLAMX) ) GOTO 2410
!
!---      Dissolved-co2 relative saturation  ---
!
          ELSEIF( ICCOMP(1,N).EQ.2 ) THEN
           PHILAX = PVA(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PLX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           PX = MAX( PGX,PLX )
           IF( CMP(3,IZN).GT.PATM ) THEN
             PCMP(N) = CMP(3,IZN)
           ELSE
             PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
           ENDIF
#ifdef ecke
           POR0(1,N) = POR0(1,N)
           POR0(2,N) = POR0(2,N)
#endif
           CALL PORSTY( N,PX,PCMP(N),PORD(2,N),PORT(2,N) )
           PORD(2,N) = MAX( PORD(2,N),EPSL )
           PORT(2,N) = MAX( PORT(2,N),PORD(2,N) )
           RHOLSX = TMS(2,N)/PORD(2,N)
           ISRX = 1
           CALL DENS_W( T(2,N),PX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved salt aqueous mass fraction  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX+RHOLAX)
           XLSX = 0.D+0
           CALL SP_B( T(2,N),XLSX,PSBX )
           PVBX = PSBX
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLAMX = 0.D+0
           DO IGC = 1, NGC
            XLAMX = XLAMX + XLCMX(IGC)
           ENDDO
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           DXLAX = 1.D-4*XLAMX
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           CALL SP_B( T(2,N),XLS(2,N),PSBX )
!
! --- Call equilibrium routine
!
           DO IGC = 1, NGC
             XMDGCX(IGC) = XMVGC(IGC,2,N)
             XMDLCX(IGC) = XMVLC(IGC,2,N)
           ENDDO
           XMDGWX = XMVGW(2,N)
           XMDLWX = XMVLW(2,N)
           CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &       XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &       XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
           DO IGC = 1, NGC
             XMVLC(IGC,2,N) = XMDLCX(IGC) 
           ENDDO
           XMVGW(2,N) = XMDGWX
           XMVLW(2,N) = XMDLWX
!
           XLAX(2,N) = 0.D+0
           DO IGC = 1, NGC
            XLAX(2,N) = XLAX(2,N) + XLCX(IGC)
           ENDDO
           XLA(2,N) = PHILAX*XLAX(2,N)
           PVA(2,N) = MIN( XLA(2,N)/XLAX(2,N),1.D+0 )*PGAX
           IF( PSBX+PVA(2,N) .GT. PGX ) THEN
             INDX = 16
             IMSG = N
             RLMSG = PSBX+PVA(2,N)
             CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &         'Conditions @ Node'
             CALL WRMSGS( INDX )
           ENDIF
           NC = 0
!
!---       Two-variable Newton-Raphson loop (YLS, XLA)  ---
!
 2510      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Salt' //
     &         ' Vol. Conc. and Dissolved-CO2 Rel. Sat.  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 2520 M = 2,4
             YLS(M,N) = YLS(2,N)
             XLA(M,N) = XLA(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             IF( M.EQ.4 ) XLA(M,N) = XLA(2,N) + DXLAX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             PVBX = PSBX
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             XLAX(M,N) = 0.D+0
             DO IGC = 1, NGC
               XLAX(M,N) = XLAX(M,N) + XLCX(IGC)
             ENDDO
             PVA(M,N) = MIN( XLA(M,N)/XLAX(M,N),1.D+0 )*PGAX
             IF( PSBX+PVA(2,N) .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX+PVA(2,N)
               CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &           'Conditions @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             CALL DENS_L( T(2,N),RHOBX,XLCX,RHOLX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
             GX(M,2) = XLA(M,N) - PHILAX*XLAX(M,N)
 2520      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RX(1,2) = (GX(4,1)-GX(2,1))/DXLAX
           RX(2,1) = (GX(3,2)-GX(2,2))/DYLSX
           RX(2,2) = (GX(4,2)-GX(2,2))/DXLAX
           RPX(1) = -GX(2,1)
           RPX(2) = -GX(2,2)
           CYLSX = (RPX(2)-RPX(1)*RX(2,2)/(RX(1,2)+SMALL))/
     &       (RX(2,1)-RX(1,1)*RX(2,2)/(RX(1,2)+SMALL))
           CXLAX = (RPX(2)-RPX(1)*RX(2,1)/(RX(1,1)+SMALL))/
     &       (RX(2,2)-RX(1,2)*RX(2,1)/(RX(1,1)+SMALL))
           YLS(2,N) = YLS(2,N) + CYLSX
           XLA(2,N) = XLA(2,N) + CXLAX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) .OR.
     &       ABS(CXLAX).GT.(1.D-6*XLAMX) ) GOTO 2510
!
!---      Dissolved-co2 aqueous mass fraction  ---
!
          ELSEIF( ICCOMP(1,N).EQ.3 ) THEN
           XLA(2,N) = PVA(2,N)
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           IF( CMP(3,IZN).GT.PATM ) THEN
             PCMP(N) = CMP(3,IZN)
           ELSE
             PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
           ENDIF
#ifdef ecke
           POR0(1,N) = POR0(1,N)
           POR0(2,N) = POR0(2,N)
#endif
           CALL PORSTY( N,PX,PCMP(N),PORD(2,N),PORT(2,N) )
           PORD(2,N) = MAX( PORD(2,N),EPSL )
           PORT(2,N) = MAX( PORT(2,N),PORD(2,N) )
           RHOLSX = TMS(2,N)/PORD(2,N)
           ISRX = 1
           CALL DENS_W( T(2,N),PX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
           CALL SOL_LS( T(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           NC = 0
!
!---       Single-variable Newton-Raphson loop (YLS)  ---
!
 2610      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Salt' //
     &         ' Vol. Conc. and Dissolved-CO2 Aqu. Mass Frac.  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 2620 M = 2,3
             YLS(M,N) = YLS(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             PVBX = PSBX
!
! --- Call equilibrium routine
!
             DO IGC = 1, NGC
               XMDGCX(IGC) = XMVGC(IGC,2,N)
               XMDLCX(IGC) = XMVLC(IGC,2,N)
             ENDDO
             XMDGWX = XMVGW(2,N)
             XMDLWX = XMVLW(2,N)
             CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &         XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,XMLCX,
     &         XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
             DO IGC = 1, NGC
               XMVLC(IGC,2,N) = XMDLCX(IGC) 
             ENDDO
             XMVGW(2,N) = XMDGWX
             XMVLW(2,N) = XMDLWX
!
             XLAX(M,N) = 0.D+0
             DO IGC = 1, NGC
               XLAX(M,N) = XLAX(M,N) + XLCX(IGC)
             ENDDO
             PVA(M,N) = MIN( XLA(M,N)/XLAX(M,N),1.D+0 )*PGAX
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             IF( PSBX+PVA(2,N) .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX+PVA(2,N)
               CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &           'Conditions @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
 2620      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RPX(1) = -GX(2,1)
           CYLSX = RPX(1)/RX(1,1)
           YLS(2,N) = YLS(2,N) + CYLSX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 2610
!
!---      Dissolved-co2 partial pressure  ---
!
          ELSEIF( ICCOMP(1,N).EQ.4 ) THEN
           PCX = 0.D+0
           IF( ISIC.EQ.1 ) THEN
             PL(2,N) = PG(2,N)-ENPR+EPSLX
           ELSEIF( ISIC.EQ.2 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ELSEIF( ISIC.EQ.3 ) THEN
             PG(2,N) = PL(2,N)+ENPR-EPSLX
           ENDIF
           PX = PL(2,N)+PATM
           PGX = PG(2,N)+PATM
           IF( CMP(3,IZN).GT.PATM ) THEN
             PCMP(N) = CMP(3,IZN)
           ELSE
             PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
           ENDIF
#ifdef ecke
           POR0(1,N) = POR0(1,N)
           POR0(2,N) = POR0(2,N)
#endif
           CALL PORSTY( N,PX,PCMP(N),PORD(2,N),PORT(2,N) )
           PORD(2,N) = MAX( PORD(2,N),EPSL )
           PORT(2,N) = MAX( PORT(2,N),PORD(2,N) )
           RHOLSX = TMS(2,N)/PORD(2,N)
           ISRX = 1
           CALL DENS_W( T(2,N),PX,RHOLWX,RHOGWX,ISRX )
!
!---       Guess dissolved-salt aqueous mass fractions  ---
!
           YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
           CALL SOL_LS( T(2,N),XLSMX )
           XLS(2,N) = MIN( YLS(2,N),XLSMX )
           DYLSX = 1.D-4*XLSMX
           NC = 0
!
!---       Single-variable Newton-Raphson loop (YLS)  ---
!
 2660      CONTINUE
           NC = NC + 1
           IF( NC.GT.32 ) THEN
             INDX = 17
             N_DB = N
             CHMSG = 'Unconverged Initial Conditions: Salt' //
     &         ' Vol. Conc. and Dissolved-CO2 Aqu. Mass Frac.  @Node: '
             CALL WRMSGS( INDX )
           ENDIF
           DO 2670 M = 2,3
             YLS(M,N) = YLS(2,N)
             IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
             XLS(M,N) = MIN( YLS(M,N),XLSMX )
             CALL SP_B( T(2,N),XLS(M,N),PSBX )
             IF( PSBX+PVA(2,N) .GT. PGX .AND. M.EQ.2 ) THEN
               INDX = 16
               IMSG = N
               RLMSG = PSBX+PVA(2,N)
               CHMSG = 'Initial Condition Transition to Unsaturated ' //
     &           'Conditions @ Node'
               CALL WRMSGS( INDX )
             ENDIF
             CALL DENS_B( T(2,N),PX,XLS(M,N),RHOBX )
             GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
 2670      CONTINUE
           RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
           RPX(1) = -GX(2,1)
           CYLSX = RPX(1)/RX(1,1)
           YLS(2,N) = YLS(2,N) + CYLSX
           IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 2660
          ENDIF
        ENDIF
        ENDIF
!
!---    Assign dissolved-salt mass fractions,
!       isobrine option  ---
!
        IF( ISLC(32).EQ.1 ) THEN
          DO 2900 M = 3,ISVC+2
            XLS(M,N) = XLS(2,N)
            YLS(M,N) = YLS(2,N)
 2900     CONTINUE
        ENDIF
!
!---    Assign temperature,
!       isothermal option  ---
!
        IF( ISLC(30).EQ.1 ) THEN
          DO 2910 M = 3,ISVC+2
            T(M,N) = T(2,N)
 2910     CONTINUE
        ENDIF
 3000 CONTINUE
 3200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_COMP group  ---
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CISC_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Compute initial solute concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/CISC_COMP'
      DO 140 NSL = 1,NSOLU
        DO 110 N = 1,NFLD
          N_DB = N
          IF( IXP(N).LE.0 ) GOTO 110
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
!---      Convert reactive species from kmol to mol  ---
!
          IF( ICT(N,NSL).LT.0 ) THEN
            ICT(N,NSL) = ABS(ICT(N,NSL))
            C(N,NSL) = C(N,NSL)*1.D+3
          ENDIF
!
!---      Convert aqueous volumetric concentration to 
!         node volumetric concentration  ---
!
          IF( ICT(N,NSL).EQ.2 ) THEN
            C(N,NSL) = C(N,NSL)*(XVS + XVL + XVG*PCGLX)
!
!---      Convert gas volumetric concentration to 
!         node volumetric concentration  ---
!
          ELSEIF( ICT(N,NSL).EQ.3 ) THEN
            C(N,NSL) = C(N,NSL)*((XVS + XVL)/PCGLX + XVG)
!
!---      Convert aqueous molarity to 
!         node volumetric concentration  ---
!
          ELSEIF( ICT(N,NSL).EQ.4 ) THEN
            C(N,NSL) = C(N,NSL)*RHOL(2,N)*(XVS + XVL + XVG*PCGLX)
!
!---      Convert aqueous molality to 
!         node volumetric concentration  ---
!
          ELSEIF( ICT(N,NSL).EQ.5 ) THEN
            C(N,NSL) = C(N,NSL)*XLW(2,N)*RHOL(2,N)*
     &        (XVS + XVL + XVG*PCGLX)
          ENDIF
!
!---      Load old-time-step concentrations  ---
!
          CO(N,NSL) = C(N,NSL)
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
!---  End of CISC_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Compute primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 XGCX(LNGC),XLCX(LNGC), XLCMX(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),XCL(LNGC)
      REAL*8 XMDGCX(LNGC)
      REAL*8 XMDLCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM_COMP'

      PGAX = 0.
      PETA = 1.D-1
      EPSLX = 1.D+0-(EPSL*1.D+2)
      PETA = 1.D-1
      EPSLX = 1.D-4
!
!---  Phase options, compute phase condition   ---
!
      DO 100 N = 1,NFLD
        IF( IXP(N).LE.0 ) GOTO 100
        IZN = IZ(N)
        N_DB = N
        N_DB = N
!
!---    Assign gas-entry pressure for non Brooks-Corey;
!       Brooks-Corey; Brooks-Corey, Dual Porosity; and
!       Brooks-Corey, Entrapment  ---
!
        IF( ISCHR(IZN).EQ.2 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.102 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
        ELSE
          ENPR = 0.D+0
        ENDIF
!
!---    Initial trapped gas saturation for the van Genuchten or
!       Brooks/Corey entrapment model  ---
!
        IF( ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 ) THEN
          SGRMX = SCHR(15,IZN)
        ELSE
          SGRMX = 0.D+0
        ENDIF
!
!---    Saturated system w/o entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - dissolved-CO2 aqueous mass fraction
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          PX = PL(2,N)+PATM
          CALL SOL_LS( T(2,N),XLSMX )
          XLSX = MIN( YLS(2,N),XLSMX )
          CALL SP_B( T(2,N),XLSX,PSBX )
          PVBX = PSBX
!
! --- Calculate dry gas mole fraction
!
           XMVA = 0.
           DO IGC = 1, NGC
             XMVA = XMVA + XLC(IGC,2,N)/GCPP(1,IGC)
           ENDDO
           IF ( XMVA.EQ.0 ) XMVA = 1.D+0
           DO IGC = 1, NGC
             XMVGC(IGC,2,N) = XLC(IGC,2,N)/GCPP(1,IGC)/XMVA
           END DO
!
!---      Call gas/aqueous equilibrium routine
!
          DO IGC = 1, NGC
            XMDGCX(IGC) = XMVGC(IGC,2,N)
            XMDLCX(IGC) = XMVLC(IGC,2,N)
          ENDDO
          XMDGWX = XMVGW(2,N)
          XMDLWX = XMVLW(2,N)
          CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &      XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &      XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &      XMDLCX,XMDLWX)
          DO IGC = 1, NGC
            XMVLC(IGC,2,N) = XMDLCX(IGC) 
          ENDDO
          XMVGW(2,N) = XMDGWX
          XMVLW(2,N) = XMDLWX
!
!---      Water-vapor pressure plus CO2 bubbling pressure exceeds
!         gas-entry pressure, transition to unsaturated condition  ---
!     
          XLAX = 0.D+0
          DO IGC = 1, NGC
           XLAX = XLAX + XLCX(IGC)
          ENDDO
          IF( XLA(2,N).GT.XLAX ) THEN
!
!---        Estimate gas saturation and transition to unsaturated
!           conditions if the gas saturation exceeds 0.001  ---
!
            PVAX = MAX( PX+ENPR-PSBX,0.D+0 )
            ISRX = 2
            CALL DENS_GC( T(2,N),PVAX,XMGCX,RHOGAX,ISRX,ZGX )
                Z1X = ZGX
            SGX = (XLA(2,N)-XLAX)*RHOL(2,N)/RHOGAX
            IF( SGX.LT.1.D-3 ) THEN
              PG(2,N) = PL(2,N)+ENPR-EPSLX
              NPHAZ(2,N) = 1
            ELSE
              SLX = 1.D+0-MIN( SGX,1.D-1 )
              SGTX = 0.D+0
              CALL CAP_COMP( IZN,SLX,SGTX,BTGL(2,N),PCX )
              PG(2,N) = PL(2,N)+PCX
!
! --- Calculate partial pressures
!
              XMVA = 0.
              DO IGC = 1, NGC
                XMVA = XMVA + XLC(IGC,2,N)/GCPP(1,IGC)
              ENDDO
              IF ( XMVA.EQ.0 ) XMVA = 1.D+0
              DO IGC = 1, NGC
                XMVGC(IGC,2,N) = XLC(IGC,2,N)/GCPP(1,IGC)/XMVA
                PVC(IGC,2,N) = PVA(2,N)*XMVGC(IGC,2,N)
              END DO
!
              NPHAZ(2,N) = 2
            ENDIF
!
!---      Water-vapor pressure plus CO2 bubbling pressure less than
!         gas-entry pressure, remain as saturated condition  ---
!
          ELSE
            PG(2,N) = PL(2,N)+ENPR-EPSLX
            NPHAZ(2,N) = 1
          ENDIF
!
!---    Unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
!
!---      Gas pressure falls below gas-entry pressure transition
!         to saturated condition  ---
!
          IF( PG(2,N).LE.(PL(2,N)+ENPR) ) THEN
!
!---        Trapped gas exists transition to saturated condition
!           w/ entrapped gas  ---
!
            IF( SGRMX.GT.EPSL .AND. SGT(2,N).GT.EPSL ) THEN
              PG(2,N) = PL(2,N)+ENPR-EPSLX
              CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
              BTGL(2,N) = 1.D+0
              IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &          BTGL(2,N) = SCHR(16,IZN)/SFTLX
              HDGL = MAX( BTGL(2,N)*(PG(2,N)-PL(2,N))/RHORL/GRAV,
     &          1.D-14 )
              REALX = REAL(ISM(IZN))
              HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
              SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
              R = 1.D+0/SGRMX - 1.D+0
              ASGTX = (1.D+0-ASLMIN(2,N))/(1.D+0+R*(1.D+0-ASLMIN(2,N)))
              SG(2,N) = ASGTX*(1.D+0-SLRX)
              SGT(2,N) = SGT(2,N)
              NPHAZ(2,N) = 3
!
!---        No trapped gas transition to saturated condition
!           w/o entrapped gas  ---
!
            ELSE
              PX = PL(2,N)+PATM
              CALL SOL_LS( T(2,N),XLSMX )
              XLSX = MIN( YLS(2,N),XLSMX )
              CALL SP_B( T(2,N),XLSX,PSBX )
              PVBX = PSBX
!
! --- Use Raoult's Law to calculate dry gas mole fractions
!
            XMVGC(1,2,N) = PVA(2,N)
            DO IGC = 2, NGC
              XMVGC(1,2,N) = XMVGC(1,2,N)-PVC(IGC,2,N)
              XMVGC(IGC,2,N) = PVC(IGC,2,N)
            ENDDO
            DO IGC = 1,NGC
              XMVGC(IGC,2,N) = XMVGC(IGC,2,N) / PVA(2,N)
            ENDDO
!
!---      Call gas/aqueous equilibrium routine
!
              DO IGC = 1, NGC
                XMDGCX(IGC) = XMVGC(IGC,2,N)
                XMDLCX(IGC) = XMVLC(IGC,2,N)
              ENDDO
              XMDGWX = XMVGW(2,N)
              XMDLWX = XMVLW(2,N)
              CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &          XGCX,XGWX,XCL,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &          XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &          XMDLCX,XMDLWX)
              DO IGC = 1, NGC
                XMVLC(IGC,2,N) = XMDLCX(IGC) 
              ENDDO
              XMVGW(2,N) = XMDGWX
              XMVLW(2,N) = XMDLWX
!
              XCLX = 0.D+0
              DO IGC = 1, NGC
               XLC(IGC,2,N) = XCL(IGC)
               XCLX  =  XCLX + XCL(IGC)
              ENDDO
              XLA(2,N) = XCLX
              PG(2,N) = PL(2,N)+ENPR-EPSLX
              NPHAZ(2,N) = 1
            ENDIF
!
!---      Gas pressure remains above gas-entry pressure, remain
!         as unsaturated condition  ---
!
          ELSE
            NPHAZ(2,N) = 2
          ENDIF
!
!---    Saturated system w/ entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - trapped gas saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
          CALL SFT_L( T(2,N),XLS(2,N),SFTLX )
          BTGL(2,N) = 1.D+0
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &      BTGL(2,N) = SCHR(16,IZN)/SFTLX
          HDGL = MAX( BTGL(2,N)*(PG(2,N)-PL(2,N))/RHORL/GRAV,1.D-14 )
          REALX = REAL(ISM(IZN))
          HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
          SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
          ASGTX = SG(2,N)/(1.D+0-SLRX)
          ASLMINX = MIN( ASLMIN(1,N),(1.D+0-ASGTX) )
          SGRX = ASGTX*(1.D+0-ASLMINX)/((1.D+0-ASLMINX)
     &      - ASGTX*ASLMINX)
!
!---      Trapped gas saturation disappears, transition to
!         saturated condition w/o entrapped gas  ---
!
          IF( SG(2,N).LT.EPSL ) THEN
            PX = PL(2,N)+PATM
            CALL SOL_LS( T(2,N),XLSMX )
            XLSX = MIN( YLS(2,N),XLSMX )
            CALL SP_B( T(2,N),XLSX,PSBX )
            PVBX = PSBX
!
!---      Call gas/aqueous equilibrium routine
!
            DO IGC = 1, NGC
              XMDGCX(IGC) = XMVGC(IGC,2,N)
              XMDLCX(IGC) = XMVLC(IGC,2,N)
            ENDDO
            XMDGWX = XMVGW(2,N)
            XMDLWX = XMVLW(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGCX,XGWX,XCL,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &        XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
            DO IGC = 1, NGC
              XMVLC(IGC,2,N) = XMDLCX(IGC) 
            ENDDO
            XMVGW(2,N) = XMDGWX
            XMVLW(2,N) = XMDLWX
!
            XCLX = 0.D+0
            DO IGC = 1, NGC
              XLC(IGC,2,N) = XCL(IGC)
              XCLX  =  XCLX + XCL(IGC)
            ENDDO
            XLA(2,N) = XCLX
            SGT(2,N) = 0.D+0
            SG(2,N) = 0.D+0
            PG(2,N) = PL(2,N)+ENPR-EPSLX
            NPHAZ(2,N) = 1
!
!---      Trapped gas saturation exceeds maximum trapped gas
!         saturation limit, transition to unsaturated condition  ---
!
          ELSEIF( (SGRX-SGRMX).GT.EPSL ) THEN
            SLX = 1.D+0-MIN( SG(2,N),SGRMX )
            SGTX = SGT(2,N)
            CALL CAP_COMP( IZN,SLX,SGTX,BTGL(2,N),PCX )
            PG(2,N) = PL(2,N)+PCX
            NPHAZ(2,N) = 2
!
!---      Trapped gas saturation neither disappears nor
!         exceeds maximum trapped gas saturation limit,
!         remain as saturated w/ entrapped gas condition  ---
!
          ELSE
            PG(2,N) = PL(2,N) + ENPR-EPSLX
            NPHAZ(2,N) = 3
          ENDIF
        ENDIF
 100  CONTINUE
!
!--- Compute increments  ---
!
      DO 300 N = 1,NFLD
        IF( IXP(N).LE.0 ) GOTO 300
        IZN = IZ(N)
        N_DB = N
!
!---    Assign gas-entry pressure for non Brooks-Corey;
!       Brooks-Corey; Brooks-Corey, Dual Porosity; and
!       Brooks-Corey, Entrapment  ---
!
        IF( ISCHR(IZN).EQ.2 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.102 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
        ELSE
          ENPR = 0.D+0
        ENDIF
!
!---    Initial trapped gas saturation for the van Genuchten or
!       Brooks/Corey entrapment model  ---
!
        IF( ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 ) THEN
          SGRMX = SCHR(15,IZN)
        ELSE
          SGRMX = 0.D+0
        ENDIF
!
!---    Energy (temperature), isothermal option  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          DNR(IEQT,N) = -1.D-4
        ENDIF
!
!---    Assign dissolved-salt mass fraction increments,
!       isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          CALL SOL_LS( T(2,N),XLSMX )
          XLSX = MIN( YLS(2,N),XLSMX )
          DNR(IEQS,N) = 1.D-5*XLSMX
        ENDIF
!
!---    Saturated w/o trapped gas  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          DNR(IEQW,N) = MAX( 1.D-1,1.D-6*(PL(2,N)+PATM) )
          PX = PL(2,N)+PATM
          CALL SP_B( T(2,N),XLSX,PSBX )
          PVBX = PSBX
!
! --- Calculate dry gas mole fraction
!
           XMVA = 0.
           DO IGC = 1, NGC
             XMVA = XMVA + XLC(IGC,2,N)/GCPP(1,IGC)
           ENDDO
           IF ( XMVA.EQ.0 ) XMVA = 1.D+0
           DO IGC = 1, NGC
             XMVGC(IGC,2,N) = XLC(IGC,2,N)/GCPP(1,IGC)/XMVA
           END DO
!
!---      Call gas/aqueous equilibrium routine
!
          DO IGC = 1, NGC
            XMDGCX(IGC) = XMVGC(IGC,2,N)
            XMDLCX(IGC) = XMVLC(IGC,2,N)
          ENDDO
          XMDGWX = XMVGW(2,N)
          XMDLWX = XMVLW(2,N)
          CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &      XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,
     &      XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &      XMDLCX,XMDLWX)
          DO IGC = 1, NGC
            XMVLC(IGC,2,N) = XMDLCX(IGC) 
          ENDDO
          XMVGW(2,N) = XMDGWX
          XMVLW(2,N) = XMDLWX
!
          XLAMX = 0.D+0
          DO IGC = 1, NGC
            XLAMX = XLAMX + XLCMX(IGC)
          ENDDO
          DO IGC = 1, NGC
          IF ( XLCMX(IGC).GT.EPSLX ) THEN
            IF( XLC(IGC,2,N).GT.(1.D-2*XLCMX(IGC)) ) THEN
              DNR(IEQGC(IGC),N) = SIGN( 1.D-4*XLCMX(IGC),5.D-1*
     &          XLCMX(IGC)-XLC(IGC,2,N) )
            ELSE
              DNR(IEQGC(IGC),N) = SIGN( 1.D-3*XLCMX(IGC),5.D-1*
     &          XLCMX(IGC)-XLC(IGC,2,N) )
            ENDIF
          ELSE
            DNR(IEQGC(IGC),N) = 1.0e-6
          ENDIF
          ENDDO
!
!---    Unsaturated w/ or w/o trapped gas  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          DNR(IEQW,N) = -MAX( 1.D-1,1.D-7*ABS(PG(2,N)-PL(2,N)) )
          DNR(IEQGC(1),N) = -DNR(IEQW,N)
          DO IGC = 2, NGC
            IF( PVC(IGC,2,N).LE.EPSL) THEN
              DNR(IEQGC(IGC),N) = 1.D-2*PL(2,N)
            ELSEIF( PVC(IGC,2,N).LT.1.D-2*PL(2,N) ) THEN
              DNR(IEQGC(IGC),N) =
     &          MAX(-DNR(IEQW,N),1.D-2*PL(2,N)-PVC(IGC,2,N))
            ELSE
              DNR(IEQGC(IGC),N) = -DNR(IEQW,N)
            END IF
          ENDDO
!
!---    Saturated w/ trapped gas  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
          DNR(IEQW,N) = MAX( 1.D-1,1.D-6*(PL(2,N)+PATM) )
          DO IGC = 1, NGC
            DNR(IEQGC(IGC),N) = SIGN( 1.D-6,5.D-1-SG(2,N) )
          ENDDO
        ENDIF
!
!--- Increment the primary variables  ---
!
        DO 200 M = 3,ISVC+2
          T(M,N) = T(2,N)
          PL(M,N) = PL(2,N)
          PG(M,N) = PG(2,N)
          DO IGC = 1, NGC
            PVC(IGC,M,N) = PVC(IGC,2,N)
          ENDDO
          XLA(M,N) = XLA(2,N)
          SG(M,N) = SG(2,N)
          YLS(M,N) = YLS(2,N)
          DO IGC = 1, NGC
            XLC(IGC,M,N) = XLC(IGC,2,N)
          ENDDO
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            IF( M.EQ.IEQT+2 ) T(M,N) = T(M,N) + DNR(IEQT,N)
          ENDIF
!
!---      Gas component
!
          IGC = M-(IEQGC(1)-1)-2
!
!---      Saturated w/o trapped gas  ---
!
          IF( NPHAZ(2,N).EQ.1 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
              PG(M,N) = PL(M,N) + ENPR-EPSLX
            ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
              XLC(IGC,M,N) = XLC(IGC,M,N) + DNR(IEQGC(IGC),N)
              XLA(M,N) = XLC(IGC,M,N)
            ELSEIF( M.GT.IEQGC(1)+2 ) THEN
              XLC(IGC,M,N) = XLC(IGC,M,N) + DNR(IEQGC(IGC),N)
              XLA(M,N) = XLA(M,N) + XLC(IGC,M,N)
!
!---        Isobrine option  ---
!
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Unsaturated w/ or w/o trapped gas  ---
!
          ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
              PG(M,N) = PG(M,N) + DNR(IEQGC(IGC),N)
            ELSEIF( M.GT.IEQGC(1)+2 ) THEN
              PVC(IGC,M,N) = PVC(IGC,M,N) + DNR(IEQGC(IGC),N)
!
!---        Isobrine option  ---
!
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Saturated w/ trapped gas  ---
!
          ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
              PG(M,N) = PL(M,N) + ENPR-EPSLX
            ELSEIF( M.GE.IEQGC(1)+2 ) THEN
              DNRX = 0.D+0
              SG(M,N) = SG(M,N) + DNR(IEQGC(IGC),N)
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
!---  End of INCRM_COMP group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KSP_COMP( N,IZN,M,PGX,PLX,SGX,SLX,RKLX,RKGX,
     &  ASLX,ASLMINX,ASGTX,SGRMX,SLRX,BTGLX,INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Salt Mode
!
!     Compute the aqueous saturation from the gas/aqueous capillary
!     pressure, and compute the aqueous relative permeability from the
!     aqueous saturation.
!
!     INDX = 0 : Trapped-gas saturation computed.
!     INDX = 1 : Trapped-gas saturation given.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/KSP_COMP'
!
!---  van Genuchten saturation function with the Webb extension  ---
!
      IF( ISCHR(IZN).EQ.1 .AND. ISM(IZN).EQ.2 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        HMPX = SCHR(9,IZN)
        SLRX = SCHR(4,IZN)
!
!---    Capillary head above the matching point head,
!       use Webb extension  ---
!
        IF( HDGL.GT.HMPX ) THEN
          SMPX = SCHR(8,IZN)
          HDGL = MIN( HDGL,HDOD )
          DMPX = SMPX/(LOG10(HMPX)-LOG10(HDOD))
          SLX = HDGL*DMPX
          ASLX = MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 )
!
!---    Capillary head at or below the matching point head,
!       use van Genuchten function
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
        ENDIF
        SGX = MAX( 1.D+0-SLX,0.D+0 )
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  van Genuchten saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.1 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
        ESLX = ASLX
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        SLX = ASLX*(1.D+0-SLRX) + SLRX
        SGX = MAX( 1.D+0-SLX,0.D+0 )
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Brooks and Corey saturation function with the Webb extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 .AND. ISM(IZN).EQ.2 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        HMPX = SCHR(9,IZN)
        SLRX = SCHR(4,IZN)
!
!---    Capillary head above the matching point head,
!       use Webb extension  ---
!
        IF( HDGL.GT.HMPX ) THEN
          SMPX = SCHR(8,IZN)
          HDGL = MIN( HDGL,HDOD )
          DMPX = SMPX/(LOG10(HMPX)-LOG10(HDOD))
          SLX = HDGL*DMPX
          ASLX = MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 )
!
!---    Capillary head at or below the matching point head,
!       use van Genuchten function
!
        ELSE
          CL = MAX( SCHR(3,IZN),SMALL )
          IF( HDGL.LE.SCHR(1,IZN) ) THEN
            ASLX = 1.D+0
          ELSE
            ASLX = (SCHR(1,IZN)/HDGL)**CL
          ENDIF
          SLRX = SCHR(4,IZN)
          SLX = ASLX*(1.D+0-SLRX) + SLRX
        ENDIF
        SGX = MAX( 1.D+0-SLX,0.D+0 )
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Brooks and Corey saturation function w/ or w/o
!     the Fayer-Simmons extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        CL = MAX( SCHR(3,IZN),SMALL )
        IF( HDGL.LE.SCHR(1,IZN) ) THEN
          ASLX = 1.D+0
        ELSE
          ASLX = (SCHR(1,IZN)/HDGL)**CL
        ENDIF
        ESLX = ASLX
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        SLX = ASLX*(1.D+0-SLRX) + SLRX
        SGX = MAX( 1.D+0-SLX,0.D+0 )
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Dual porosity van Genuchten saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.3 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLR_M = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        SDPM(N) = ASL_M*(1.D+0-SLR_M) + SLR_M
        CN = MAX( SCHR(6,IZN),SMALL )
        IF( SCHR(15,IZN).LE.EPSL ) THEN
          IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(15,IZN)
        ENDIF
        ASL_F = (1.D+0/(1.D+0 + (SCHR(5,IZN)*HDGL)**CN))**CM
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLR_F = MAX( (1.D+0-HSCL)*SCHR(7,IZN),ZERO )
        SDPF(N) = ASL_F*(1.D+0-SLR_F) + SLR_F
        PORD_MX = (1.D+0-POR(4,IZN))*POR(2,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
        PORD_FX = POR(4,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
        SLX = SDPF(N)*PORD_FX + SDPM(N)*PORD_MX
        SGX = MAX( 1.D+0-SLX,0.D+0 )
        ASLX = ASL_F*PORD_FX + ASL_M*PORD_MX
        ESLX = ASLX
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Dual porosity Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.4 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        CL = MAX( SCHR(3,IZN),SMALL )
        IF( HDGL.LE.SCHR(1,IZN) ) THEN
          ASL_M = 1.D+0
        ELSE
          ASL_M = (SCHR(1,IZN)/HDGL)**CL
        ENDIF
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLR_M = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        SDPM(N) = ASL_M*(1.D+0-SLR_M) + SLR_M
        CL = MAX( SCHR(6,IZN),SMALL )
        IF( HDGL.LE.SCHR(5,IZN) ) THEN
          ASL_F = 1.D+0
        ELSE
          ASL_F = (SCHR(5,IZN)/HDGL)**CL
        ENDIF
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLR_F = MAX( (1.D+0-HSCL)*SCHR(7,IZN),ZERO )
        SDPF(N) = ASL_F*(1.D+0-SLR_F) + SLR_F
        PORD_MX = (1.D+0-POR(4,IZN))*POR(2,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
        PORD_FX = POR(4,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
        SLX = SDPF(N)*PORD_FX + SDPM(N)*PORD_MX
        SGX = MAX( 1.D+0-SLX,0.D+0 )
        ASLX = ASL_F*PORD_FX + ASL_M*PORD_MX
        ESLX = ASLX
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Haverkamp saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.5 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( HDGL.LE.SCHR(1,IZN) ) THEN
          ASLX = 1.D+0
        ELSE
          ASLX = SCHR(2,IZN)/(SCHR(2,IZN)
     &      + (((HDGL-SCHR(1,IZN))/SCHR(5,IZN))**SCHR(3,IZN)))
        ENDIF
        ESLX = ASLX
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        SLX = ASLX*(1.D+0-SLRX) + SLRX
        SGX = MAX( 1.D+0-SLX,0.D+0 )
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Linear or linear-log interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.10 .OR. ISCHR(IZN).EQ.12 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( ISCHR(IZN).EQ.12 ) HDGL = LOG(HDGL)
        ITBX = 0
        IF( M.NE.2 ) ITBX = 1
        ASLX = FNTBLY( HDGL,ISLTBL(1,IZN),ISLTBL(2,IZN),ITBX )
        ESLX = ASLX
        SLX = ASLX
        SGX = MAX( 1.D+0-SLX,0.D+0 )
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Cubic-spline or cubic-spline-log interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.13 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( ISCHR(IZN).EQ.13 ) HDGL = LOG(HDGL)
        ASLX = FSPLNY( HDGL,ISLTBL(1,IZN),ISLTBL(2,IZN) )
        ESLX = ASLX
        SLX = ASLX
        SGX = MAX( 1.D+0-SLX,0.D+0 )
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  van Genuchten saturation function w/ gas entrapment  ---
!
      ELSEIF( ISCHR(IZN).EQ.101 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        ASLX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGL)**CN))**CM
!
!---    Compute trapped gas saturation, using the minimum
!       apparent aqueous saturation  ---
!
        IF( INDX.EQ.0 ) THEN
          ASLM = MIN( ASLX,ASLMINX )
          IF( SGRMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
            R = 1.D+0/SGRMX - 1.D+0
            ASGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &        (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
          ELSE
            ASGTX = 0.D+0
          ENDIF
!
!---    Trapped-gas saturation given, compute the minimum
!       apparent aqueous saturation  ---
!
        ELSEIF( INDX.EQ.1 ) THEN
          IF( SGRMX.GT.EPSL ) THEN
            R = 1.D+0/SGRMX - 1.D+0
            ASGTX = SGX/(1.D+0-SLRX)
            ASLM = (ASGTX + 2.D+0*ASGTX*R - ASGTX*R*ASLX + ASGTX*R*R
     &        - ASGTX*R*R*ASLX - ASLX)/(ASGTX*R + ASGTX*R*R 
     &        - ASGTX*R*R*ASLX - 1.D+0)
            ASLM = MIN( MAX( ASLM,0.D+0 ),ASLX )
            ASLMINX = ASLM
          ELSE
            ASGTX = 0.D+0
            ASLM = ASLX
            ASLMINX = ASLM
          ENDIF
!
!---    Unrecognized option  ---
!
        ELSE
          CHMSG = 'Unrecognized Trapped-Gas Option'
          INDX = 12
          IMSG = INDX
          CALL WRMSGS( INDX )
        ENDIF
        ESLX = ASLX - ASGTX
        SLX = ESLX*(1.D+0-SLRX) + SLRX
        SGX = MAX( 1.D+0-SLX,0.D+0 )
!
!---  Brooks and Corey saturation function w/ gas entrapment  ---
!
      ELSEIF( ISCHR(IZN).EQ.102 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        CL = MAX( SCHR(3,IZN),SMALL )
        IF( HDGL.LE.SCHR(1,IZN) ) THEN
          ASLX = 1.D+0
        ELSE
          ASLX = (SCHR(1,IZN)/HDGL)**CL
        ENDIF
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
!
!---    Compute trapped gas saturation, using the minimum
!       apparent aqueous saturation  ---
!
        IF( INDX.EQ.0 ) THEN
          ASLM = MIN( ASLX,ASLMINX )
          IF( SGRMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
            R = 1.D+0/SGRMX - 1.D+0
            ASGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &        (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
          ELSE
            ASGTX = 0.D+0
          ENDIF
!
!---    Trapped-gas saturation given, compute the minimum
!       apparent aqueous saturation  ---
!
        ELSEIF( INDX.EQ.1 ) THEN
          IF( SGRMX.GT.EPSL ) THEN
            R = 1.D+0/SGRMX - 1.D+0
            ASGTX = SGX/(1.D+0-SLRX)
            ASLM = (ASGTX + 2.D+0*ASGTX*R - ASGTX*R*ASLX + ASGTX*R*R
     &        - ASGTX*R*R*ASLX - ASLX)/(ASGTX*R + ASGTX*R*R 
     &        - ASGTX*R*R*ASLX - 1.D+0)
            ASLM = MIN( MAX( ASLM,0.D+0 ),ASLX )
            ASLMINX = ASLM
          ELSE
            ASGTX = 0.D+0
            ASLM = ASLX
            ASLMINX = ASLM
          ENDIF
!
!---    Unrecognized option  ---
!
        ELSE
          CHMSG = 'Unrecognized Trapped-Gas Option'
          INDX = 12
          IMSG = INDX
          CALL WRMSGS( INDX )
        ENDIF
        ESLX = ASLX - ASGTX
        SLX = ESLX*(1.D+0-SLRX) + SLRX
        SGX = MAX( 1.D+0-SLX,0.D+0 )
      ENDIF
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
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101  ) THEN
          RKLX(1) = SQRT(SLPX)*((1.D+0-(1.D+0-SLPX**(1.D+0/RPLC(2,IZN)))
     &      **RPLC(2,IZN))**2)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 ) THEN
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
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 ) THEN
          RKLX(1) = SQRT(ESLX)*((1.D+0-(1.D+0-ESLX**(1.D+0/RPLC(2,IZN)))
     &      **RPLC(2,IZN))**2)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 ) THEN
          RKLX(1) = ESLX**(2.5D+0 + 2.0D+0/RPLC(2,IZN))
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
          RKLX(1) = (ESLX**2)*(1.D+0-(1.D+0-ESLX**(1.D+0/RPLC(2,IZN)))
     &      **RPLC(2,IZN))
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 ) THEN
          RKLX(1) = ESLX**(3.0D+0 + 2.0D+0/RPLC(2,IZN))
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
        RKLX(1) = ESLX**4
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Fatt and Klikoff relative permeability function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.4 ) THEN
        RKLX(1) = ESLX**3
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Haverkamp relative permeability function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.5 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR. 
     &    ISCHR(IZN).EQ.5 ) THEN
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
        RKLX(1) = RPLC(1,IZN)*(ESLX**RPLC(2,IZN))
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
!---  Tabular function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.10 ) THEN
        ITBX = 0
        IF( M.NE.2 ) ITBX = 1
        RKLX(1) = FNTBLY( SLX,IRLTBL(1,IZN),IRLTBL(2,IZN),ITBX )
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
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 ) THEN
          SGP = 1.D+0-ESLX
          RKGX = SQRT(SGP)*((1.D+0-ESLX**(1.D+0/RPGC(3,IZN)))
     &      **RPGC(3,IZN))**2
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 ) THEN
          SGP = 1.D+0 - ESLX
          RKGX = SQRT(SGP)*(1.D+0-ESLX**(1.D+0+1.D+0/RPGC(3,IZN)))**2
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
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 ) THEN
          SGP = 1.D+0-ESLX
          RKGX = (SGP**2)*((1.D+0-ESLX**(1.D+0/RPGC(3,IZN)))
     &      **RPGC(3,IZN))
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 ) THEN
          SGP = 1.D+0 - ESLX
          RKGX = (SGP**2)*(1.D+0-ESLX**(1.D+0 + 2.D+0/RPGC(3,IZN)))
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
        RKGX = ((1.D+0-SLPX)**2)*(1.D+0-SLPX**2)
!
!---  Fatt and Klikoff relative permeability function  ---
!
      ELSEIF( IRPG(IZN).EQ.4 ) THEN
        RKGX = (1.D+0-ESLX)**3
!
!---  Free-Corey relative permeability function  ---
!
      ELSEIF( IRPG(IZN).EQ.7 ) THEN
        SLRX = RPGC(3,IZN)
        SGRX = RPGC(4,IZN)
        SGPX = MIN( MAX( (1.D+0-SLX-SGRX)/(1.D+0-SLRX-SGRX),0.D+0 ),
     &    1.D+0 )
        RKGX = RPGC(1,IZN)*(SGPX**(RPGC(2,IZN)))
      ENDIF
      RKLX(1) = MAX( RKLX(1),1.D-24 )
      RKLX(2) = MAX( RKLX(2),1.D-24 )
      RKLX(3) = MAX( RKLX(3),1.D-24 )
      RKGX = MAX( RKGX,0.D+0 )
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KSP_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Load the current time step values into the old time step
!     variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      USE FDVGC
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
      SUB_LOG(ISUB_LOG) = '/LDO_COMP'
!
!---  Assign old time step values  ---
!
      DO 100 N = 1,NFLD
        ASLMIN(1,N) = MIN( ASL(N),ASLMIN(2,N) )
        ASLMIN(2,N) = ASLMIN(1,N)
        BTGL(1,N) = BTGL(2,N)
        DFGW(1,N) = DFGW(2,N)
        NPHAZ(1,N) = NPHAZ(2,N)
        PG(1,N) = PG(2,N)
        PL(1,N) = PL(2,N)
        PORD(1,N) = PORD(2,N)
        PORT(1,N) = PORT(2,N)
        POSM(1,N) = POSM(2,N)
        PVA(1,N) = PVA(2,N)
        PVW(1,N) = PVW(2,N)
        RHOG(1,N) = RHOG(2,N)
        RHOL(1,N) = RHOL(2,N)
        RHOSP(1,N) = RHOSP(2,N)
        RKG(1,N) = RKG(2,N)
        RKL(1,1,N) = RKL(1,2,N)
        RKL(2,1,N) = RKL(2,2,N)
        RKL(3,1,N) = RKL(3,2,N)
        TMS(1,N) = TMS(2,N)
        SG(1,N) = SG(2,N)
        SL(1,N) = SL(2,N)
        SS(1,N) = SS(2,N)
        T(1,N) = T(2,N)
        TORG(1,N) = TORG(2,N)
        TORL(1,N) = TORL(2,N)
        VISG(1,N) = VISG(2,N)
        VISDG(1,N) = VISDG(2,N)
        VISL(1,N) = VISL(2,N)
        XGA(1,N) = XGA(2,N)
        XGW(1,N) = XGW(2,N)
        XLA(1,N) = XLA(2,N)
        XLS(1,N) = XLS(2,N)
        XLW(1,N) = XLW(2,N)
        XMGW(1,N) = XMGW(2,N)
        YLS(1,N) = YLS(2,N)
        XGO(1,N) = XGO(2,N)
        XMGO(1,N) = XMGO(2,N)
        DO IGC = 1, NGC
          PVC(IGC,1,N) = PVC(IGC,2,N)
          XGC(IGC,1,N) = XGC(IGC,2,N)
          XMGC(IGC,1,N) = XMGC(IGC,2,N)
          XLC(IGC,1,N) = XLC(IGC,2,N)
          XMLC(IGC,1,N) = XMLC(IGC,2,N)
          DFLC(IGC,1,N) = DFLC(IGC,2,N)
          DFGC(IGC,1,N) = DFGC(IGC,2,N)
          HGC(IGC,1,N) = HGC(IGC,2,N)
        ENDDO
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          DO IGC = 1, NGC
           HGC(IGC,1,N) = HGC(IGC,2,N)
          ENDDO
          HGA(1,N) = HGA(2,N)
          HGW(1,N) = HGW(2,N)
          HG(1,N) = HG(2,N)
          HGD(1,N) = HGD(2,N)
          HL(1,N) = HL(2,N)
          HSP(1,N) = HSP(2,N)
          THKL(1,N) = THKL(2,N)
          THKG(1,N) = THKG(2,N)
          UEG(1,N) = UEG(2,N)
          PSO(1,N) = PSO(2,N)
        ENDIF
        DO 90 NSL = 1,NSOLU+NEQC+NEQK
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO_COMP group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Compute hydrologic, thermodynamic and physical properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
       REAL*8 VISGC(LNGC)
       REAL*8 RHOGC(LNGC)
       REAL*8 XGCX(LNGC),XLCX(LNGC)
       REAL*8 XMGCX(LNGC),XMLCX(LNGC),XMGX(LNGC)
       REAL*8 THKGCX(LNGC),UEGCX(LNGC)
       REAL*8 XMDGCX(LNGC),HGCX(LNGC),XCLX(LNGC)
       REAL*8 XMDLCX(LNGC)
       REAL*8 DFGCWX(LNGC),DFGC2X(LNGC,LNGC)
       REAL*8 DHCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROP_COMP'
      DO 200 N = 1,NFLD
        IF( IXP(N).LE.0 ) GOTO 200
        IZN = IZ(N)
        N_DB = N
#ifdef ecke
        POR0(1,N) = POR0(1,N)
        POR0(2,N) = POR0(2,N)
#endif
!
!---    Assign gas-entry pressure for non Brooks-Corey;
!       Brooks-Corey; Brooks-Corey, Dual Porosity; and
!       Brooks-Corey, Entrapment  ---
!
        IF( ISCHR(IZN).EQ.2 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.102 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
        ELSE
          ENPR = 0.D+0
        ENDIF
!
!---    Initial trapped gas saturation for the van Genuchten or
!       Brooks/Corey entrapment model  ---
!
        IF( ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 ) THEN
          SGRMX = SCHR(15,IZN)
        ELSE
          SGRMX = 0.D+0
        ENDIF
!
!---    Minimum apparent aqueous saturation  ---
!
        ASLMINX = ASLMIN(1,N)
!
!---    Loop over increment indices  ---
!
        DO 100 M = 2,ISVC+2
          PLX = PL(M,N)+PATM
          PGX = PG(M,N)+PATM
          PX = MAX( PGX,PLX )
!
!---      Saturated system w/o entrapped gas
!         Water mass - aqueous pressure
!         CO2 mass - dissolved-CO2 aqueous mass fraction
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          IF( NPHAZ(2,N).EQ.1 ) THEN
            PCX = 0.D+0
            CALL SOL_LS( T(M,N),XLSMX )
            XLS(M,N) = MIN(YLS(M,N),XLSMX)
            CALL SP_B( T(M,N),XLS(M,N),PSBX )
            CALL DENS_B( T(M,N),PX,XLS(M,N),RHOBX )
            IF( ISLC(44).EQ.1 ) THEN
              PVBX = PSBX
            ELSE
              CALL VPL_B( T(M,N),PSBX,PCX,RHOBX,PVBX,XLS(M,N) )
            ENDIF
!
! --- Calculate dry gas mole fraction
!
           XMVA = 0.
           DO IGC = 1, NGC
             XMVA = XMVA + XLC(IGC,M,N)/GCPP(1,IGC)
           ENDDO
           IF ( XMVA.EQ.0 ) XMVA = 1.D+0
           DO IGC = 1, NGC
             XMVGC(IGC,M,N) = XLC(IGC,M,N)/GCPP(1,IGC)/XMVA
           END DO
!
!---      Call gas/aqueous equilibrium routine
!
            DO IGC = 1, NGC
              XMDGCX(IGC) = XMVGC(IGC,M,N)
              XMDLCX(IGC) = XMVLC(IGC,M,N)
            ENDDO
            XMDGWX = XMVGW(M,N)
            XMDLWX = XMVLW(M,N)
            CALL EQUIL( T(M,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGCX,XGWX,XLCX,XLS(M,N),XLWX,XMGCX,XMGWX,
     &        XMLCX,XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &        XMDLCX,XMDLWX)
            DO IGC = 1, NGC
              XMVLC(IGC,M,N) = XMDLCX(IGC) 
            ENDDO
            XMVGW(M,N) = XMDGWX
            XMVLW(M,N) = XMDLWX
             XLA(M,N) = 0.D+0
             DO IGC = 1, NGC
              XLA(M,N) = XLA(M,N) + XLC(IGC,M,N)
             ENDDO
            XLW(M,N) = 1.D+0-XLA(M,N)-XLS(M,N)
            WTMLX = 1.D+0/(XLA(M,N)/WTMGCX+XLS(M,N)/WTMS+XLW(M,N)/WTMW)
            XMLA(M,N) = 0.D+0
            DO IGC = 1, NGC
              XMLC(IGC,M,N) = XLC(IGC,M,N)*WTMLX/GCPP(1,IGC)
              XMLA(M,N) = XMLC(IGC,M,N) + XMLA(M,N)
            ENDDO
            XMLS(M,N) = XLS(M,N)*WTMLX/WTMS
            XMLW(M,N) = XLW(M,N)*WTMLX/WTMW
            DO IGC = 1, NGC
              IF(XMLCX(IGC).LE.SMALL) THEN
                XMGC(IGC,M,N)=0.D+0
              ELSE
                XMGC(IGC,M,N) = XMLC(IGC,M,N)*XMGCX(IGC)/XMLCX(IGC)
              END IF
              IF(XMLCX(IGC).LE.SMALL) XMGC(IGC,M,N)=0.D+0
            ENDDO
            XMGA(M,N) = 0.D+0
            DO IGC = 1, NGC
             XMGA(M,N) = XMGA(M,N) + XMGC(IGC,M,N)
            ENDDO
            XMGW(M,N) = 1.D+0-XMGA(M,N)
            WTMGX = XMGA(M,N)*WTMGCX + XMGW(M,N)*WTMW
            DO IGC = 1, NGC
             XGC(IGC,M,N) = XMGC(IGC,M,N)*GCPP(1,IGC)/WTMGX
            ENDDO
            XGA(M,N) = 0.D+0
            DO IGC = 1, NGC
             XGA(M,N) = XGA(M,N) + XGC(IGC,M,N)
            ENDDO
            XGW(M,N) = 1.D+0-XGA(M,N)
            ISRX = 2
            CALL DENS_W( T(M,N),PVBX,RHOLWX,RHOGWX,ISRX )
            SL(M,N) = 1.D+0
             DO IGC = 1, NGC
                 PVC(IGC,M,N) = XMGC(IGC,M,N)*PGAX
             ENDDO
             PVA(M,N) = 0.D+0
             DO IGC = 1,NGC
               PVA(M,N) = PVA(M,N) + PVC(IGC,M,N)
             ENDDO
             PVW(M,N) = PVBX
!
!---      Unsaturated system w/o entrapped gas
!         Water mass - aqueous pressure
!         CO2 mass - gas pressure
!         NaCl mass - total NaCl aqueous mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
            PCX = PG(M,N)-PL(M,N)
            CALL SOL_LS( T(M,N),XLSMX )
            XLS(M,N) = MIN(YLS(M,N),XLSMX)
            CALL SP_B( T(M,N),XLS(M,N),PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(M,N),PX,XLS(M,N),RHOBX )
            IF( ISLC(44).EQ.1 ) THEN
              PVBX = PSBX
            ELSE
              CALL VPL_B( T(M,N),PSBX,PCX,RHOBX,PVBX,XLS(M,N) )
            ENDIF
            ISRX = 2
            CALL DENS_W( T(M,N),PVBX,RHOLWX,RHOGWX,ISRX )
            PVA(M,N) = MAX( PX-PVBX,0.D+0 )
            PVW(M,N) = PVBX
            PVX = PVA(M,N) + PVW(M,N)
            PVC(1,M,N) = PVA(M,N)
            DO IGC = 2,NGC
              PVC(1,M,N) = PVC(1,M,N)-PVC(IGC,M,N)
            END DO
!
! --- Use Raoult's Law to calculate dry gas mole fractions
!
            XMVGC(1,M,N) = PVA(M,N)
            DO IGC = 2, NGC
              XMVGC(1,M,N) = XMVGC(1,M,N)-PVC(IGC,M,N)
              XMVGC(IGC,M,N) = PVC(IGC,M,N)
            ENDDO
            DO IGC = 1,NGC
              XMVGC(IGC,M,N) = XMVGC(IGC,M,N) / PVA(M,N)
            ENDDO
!
!--- Call equilibrium routine here!
!
            DO IGC = 1, NGC
              XMDGCX(IGC) = XMVGC(IGC,M,N)
              XMDLCX(IGC) = XMVLC(IGC,M,N)
            ENDDO
            XMDGWX = XMVGW(M,N)
            XMDLWX = XMVLW(M,N)
            CALL EQUIL( T(M,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGCX,XGW(M,N),XLCX,XLS(M,N),XLW(M,N),XMGCX,
     &        XMGW(M,N),XMLCX,XMLS(M,N),XMLW(M,N),ZGX,WTMGCX,XMDGCX,
     &        XMDGWX,XMDLCX,XMDLWX)
            DO IGC = 1, NGC
              XMVLC(IGC,M,N) = XMDLCX(IGC) 
            ENDDO
            XMVGW(M,N) = XMDGWX
            XMVLW(M,N) = XMDLWX
!
            XLA(M,N) = 0.D+0
            XMLA(M,N) = 0.D+0
            XGA(M,N) = 0.D+0
            XMGA(M,N) = 0.D+0
            DO IGC = 1, NGC
             XGC(IGC,M,N) = XGCX(IGC)
             XGA(M,N) = XGA(M,N) + XGC(IGC,M,N)
             XLC(IGC,M,N) = XLCX(IGC)
             XLA(M,N) = XLA(M,N) + XLC(IGC,M,N)
             XMGC(IGC,M,N) = XMGCX(IGC)
             XMGA(M,N) = XMGA(M,N) + XMGC(IGC,M,N)
             XMLC(IGC,M,N) = XMLCX(IGC)
             XMLA(M,N) = XMLA(M,N) + XMLC(IGC,M,N)
            ENDDO
      
          ENDIF
!
!---      Porous-media porosity  ---
!
          CALL PORSTY( N,PX,PCMP(N),PORD(M,N),PORT(M,N) )
          PORD(M,N) = MAX( PORD(M,N),EPSL )
          PORT(M,N) = MAX( PORT(M,N),PORD(M,N) )
!
!---      Surface tension and saturation  ---
!
          CALL SFT_L( T(M,N),XLS(M,N),SFTLX )
          BTGL(M,N) = 1.D+0
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &      BTGL(M,N) = SCHR(16,IZN)/SFTLX
          INDX = 0
          IF( NPHAZ(2,N).EQ.3 ) THEN
            SGT(M,N) = SG(M,N)
            INDX = 1
          ENDIF
          CALL KSP_COMP( N,IZN,M,PG(M,N),PL(M,N),SG(M,N),SL(M,N),
     &      RKL(1,M,N),RKG(M,N),ASLX,ASLMINX,ASGTX,SGRMX,
     &      SLRX,BTGL(M,N),INDX )
          IF( M.EQ.2 ) THEN
            ASL(N) = ASLX
            ASGT(N) = ASGTX
            ASLMIN(2,N) = ASLMINX
          ENDIF
          SGT(M,N) = ASGTX*(1.D+0-SLRX)
!
!---      Gas density and component fractions  ---
!
          ISRX = 2
          CALL DENS_GC( T(M,N),PVA(M,N),XMGCX,RHOGAX,ISRX,ZGX )
                Z1X = ZGX
          RHOG(M,N) = RHOGWX+RHOGAX
          WTMGX = 1.D+0/(XGA(M,N)/WTMGCX + XGW(M,N)/WTMW)
          RHOMG(M,N) = RHOG(M,N)/WTMGX
!
!---      Gas viscosity  ---
!
            DO IGC = 1, NGC
             XMGCX(IGC)=XMGC(IGC,M,N)
            ENDDO
            CALL VISC_GC( T(M,N),PX,ZGX,XMGCX,RHOGAX,VISGC,
     &                   VISGCX )
          CALL VISC_W( T(M,N),PVW(M,N),RHOGWX,VISGWX )
          CALL VISC_G( VISGC,VISGWX,XMGCX,XMGW(M,N),VISG(M,N) )
          VISG(M,N)=VISG(M,N)+VISGCX
          VISDG(M,N) = VISGCX
!
!---      Water-vapor diffusion coefficient  ---
!
          CALL DIFC_GW( T(M,N),PX,DFGW(M,N) )
!
!---      Aqueous component fractions and density  ---
!
          XLW(M,N) = MAX( 1.D+0-XLS(M,N)-XLA(M,N),0.D+0 )
          WTMLX = 1.D+0/(XLA(M,N)/WTMGCX + XLW(M,N)/WTMW + XLS(M,N)/
     &             WTMS)
          DO IGC = 1, NGC
            XLCX(IGC) = XLC(IGC,M,N)
          ENDDO
          XMLW(M,N) = XLW(M,N)*WTMLX/WTMW
          XMLS(M,N) = XLS(M,N)*WTMLX/WTMS
          CALL DENS_L( T(2,N),RHOBX,XLCX,RHOL(M,N) )
          RHOML(M,N) = RHOL(M,N)/WTMLX
!
!---      Aqueous viscosity  ---
!
          ISRX = 1
          CALL DENS_W( T(M,N),PX,RHOLWX,RHOX,ISRX )
          CALL VISC_W( T(M,N),PX,RHOLWX,VISLWX )
          CALL VISC_B( T(M,N),XLS(M,N),VISLWX,VISBX )
          DO IGC = 1, NGC
           XLCX(IGC) = XMLC(IGC,M,N)
          ENDDO
          CALL VISC_L( XLCX,VISBX,VISGC,VISL(M,N) )
!
!---      Gas diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            IF( NSTEP-NRST.EQ.0 ) THEN
              DO IGC = 1,NGC
                DFGC(IGC,M,N) = DFGAC
              ENDDO
            ENDIF
          ELSEIF( ISLC(2).EQ.2 ) THEN
            DO IGC = 1,NGC
              CALL BNDFGCW( T(M,N),PGX,DFGCWX(IGC),IGC )
            ENDDO
            DO IGC = 1,NGC-1
              DO JGC = IGC+1,NGC
                CALL BNDFGC2( T(M,N),PGX,DFGC2X(IGC,JGC),IGC,JGC )
              ENDDO
            ENDDO
            DO IGC = 1,NGC
              IF( (1.D+0-XMGC(IGC,M,N)).LT.EPSL ) THEN
                DFGCNX = 2.D+0*DFGCWX(IGC)
                DFGCDX = DFGCWX(IGC)
                DO JGC = 1,NGC
                  IF( JGC.NE.IGC ) THEN
                    MNGC = MIN(IGC,JGC)
                    MXGC = MAX(IGC,JGC)
                    DFGCNX = DFGCNX*DFGC2X(MNGC,MXGC)
                    DFGCDX = DFGCDX + DFGC2X(MNGC,MXGC)
                  ENDIF
                ENDDO
                DFGC(IGC,M,N) = DFGCNX/DFGCDX
              ELSE
                DFGCNX = 1.D+0-XMGC(IGC,M,N)
                DFGCDX = XMGW(M,N)/DFGCWX(IGC)
                DO JGC = 1,NGC
                  IF( JGC.NE.IGC ) THEN
                    MNGC = MIN(IGC,JGC)
                    MXGC = MAX(IGC,JGC)
                    DFGCDX = DFGCDX + XMGC(JGC,M,N)/DFGC2X(MNGC,MXGC)
                  ENDIF
                ENDDO
                DFGC(IGC,M,N) = DFGCNX/DFGCDX
              ENDIF
            ENDDO
          ENDIF
!
!---      Aqueous diffusion coefficients  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            IF( NSTEP-NRST.EQ.0 ) THEN
              DO IGC = 1,NGC
                DFLC(IGC,M,N) = DFLAC
              ENDDO
            DFLS(M,N) = DFLSC
            ENDIF
          ELSEIF( ISLC(4).EQ.2 ) THEN
            DO IGC = 1,NGC
              CALL BNDFLC( T(M,N),VISL(M,N),VISGC(IGC),
     &          DFLC(IGC,M,N),IGC )
            ENDDO
            CALL DIFC_LS( T(M,N),XLS(M,N),VISL(M,N),DFLS(M,N) )
          ENDIF
!
!---      Precipitated NaCl density, saturation, and enthalpy  ---
!
          CALL DENS_S( T(M,N),PX,RHOSP(M,N) )
          SS(M,N) = MAX(YLS(M,N)-XLS(M,N),0.D+0)*RHOL(M,N)*SL(M,N)/
     &      RHOSP(M,N)
!
!---      Permeability reduction factor  ---
!
          IF( IPRF(IZN).EQ.1 ) CALL PERM_R( SS(M,N),PERMRF(M,N),
     &      PORD(M,N),IZN )
!
!---      Change in porosity as a function of precipitated salt  ---
!
          PORD(M,N) = PORD(M,N)*MAX( (1.D+0-SS(M,N)),1.D-6 )
!
!---      Aqueous and gas tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL(M,N),SG(M,N),ZERO,
     &      PORD(M,N),TORL(M,N),TORG(M,N),TORNX )
!
!---      NaCl volumetric concentration  ---
!
          TMS(M,N) = YLS(M,N)*RHOL(M,N)*SL(M,N)*PORD(M,N)
!
!---      Nonisothermal simulation  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Gas enthalpy and internal energy  ---
!
            ISRX = 2
            WTMX = 0.D+0
            DO IGC = 1, NGC
              WTMX = WTMX + XMGC(IGC,M,N)*GCPP(1,IGC)
            ENDDO
            PVX = 0.D+0
            DO IGC = 1, NGC
             PVX = PVX + PVC(IGC,M,N)
            ENDDO
            DO IGC = 1, NGC
              XMGX(IGC) = XMGC(IGC,M,N)
            ENDDO
            UEGX = 0.D+0
            IF ( PVA(M,N).NE.0 ) THEN
              CALL ENTH_GC( T(M,N),PVA(M,N),XMGX,XMGW(M,N),
     &          WTMX,HGD(M,N),UEG(M,N),ISRX )
            END IF
            CALL ENTH_W( T(M,N),PVW(M,N),HLWX,HGW(M,N),ISRX )
            UGWX = HGW(M,N) - PVW(M,N)/RHOGWX
                UEGX = UEG(M,N)
            UEG(M,N) = UEG(M,N)+XGW(M,N)*UGWX
            HGA(M,N) = HGD(M,N)
            HG(M,N) = XGW(M,N)*HGW(M,N)
            HG(M,N) = XGA(M,N)*HGD(M,N) + HG(M,N)
!
!---        Gas thermal conductivity  ---
!
            ISRX = 2
            DO IGC = 1, NGC
              XMGX(IGC) = XMGC(IGC,M,N)
            ENDDO
              CALL THK_GC( PGX,RHOGAX,THKGAX,T(M,N),XMGX,Z1X )
            CALL THK_W( T(M,N),PGX,RHOGWX,THKGWX )
            DO IGC = 1, NGC
              XMGX(IGC) = XMGC(IGC,M,N)
            ENDDO
            CALL THK_G( T(M,N),THKGAX,THKGWX,XMGX,
     &        XMGW(M,N),THKG(M,N) )
!
!---        Aqueous enthalpy and internal energy  ---
!
            ISRX = 1
            CALL ENTH_W( T(M,N),PVW(M,N),HLWX,HGW(M,N),ISRX )
            UGWX = HGW(M,N) - PVW(M,N)/RHOGWX
            UEG(M,N) = XGW(M,N)*UGWX+XGA(M,N)*UEGX
            HG(M,N) = XGW(M,N)*HGW(M,N)+XGA(M,N)*HGA(M,N)
            HGD(M,N)=HGA(M,N)
!
!---        Gas thermal conductivity  ---
!
            ISRX = 2
            DO IGC = 1, NGC
              XMGX(IGC) = XMGC(IGC,M,N)
            ENDDO
            CALL THK_GC( PGX,RHOGAX,THKGAX,T(M,N),XMGX,Z1X )
            CALL THK_W( T(M,N),PGX,RHOGWX,THKGWX )
            DO IGC = 1, NGC
              XMGX(IGC) = XMGC(IGC,M,N)
            ENDDO
            CALL THK_G( T(M,N),THKGAX,THKGWX,XMGX,
     &        XMGW(M,N),THKG(M,N) )
!
!---        Aqueous enthalpy and internal energy  ---
!
            ISRX = 1
            CALL ENTH_W( T(M,N),PX,HLWX,HGWX,ISRX )
            CALL ENTH_B( T(M,N),XLS(M,N),HLWX,HBX )
!
!---        Change in solubility with temperature ---
!
            DTX = 1.D-6
            CALL EQUIL( T(M,N)+DTX,PX,PGAX,PGWX,PSBX,PVBX,
     &        XGCX,XGW(M,N),XLCX,XLS(M,N),XLW(M,N),XMGCX,
     &        XMGW(M,N),XMLCX,XMLS(M,N),XMLW(M,N),ZGX,WTMGCX,XMDGCX,
     &        XMDGWX,XMDLCX,XMDLWX)
            DO IGC = 1, NGC
!              IF( XMLCX(IGC).GT.EPSL ) THEN
!                DHCX(IGC) = LOG(XMLC(IGC,M,N)/XMLCX(IGC))/DTX
!              ELSE
                DHCX(IGC) = 0.D+0
!              ENDIF
            END DO
            CALL ENTH_L( T(M,N),XLS(M,N),XLCX,DHCX,
     &        HBX,HGD(M,N),HL(M,N) )
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
  100   CONTINUE
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDBC_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-NComponent-NaCl-Energy Mode (STOMP-WNSE)
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
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      CHARACTER*64 ADUM,BDUM(LUK+LSOLU+2),FDUM,SDUM
      CHARACTER*64 UNTS
      CHARACTER*32 CHTYP(44)
      CHARACTER*512 CHDUM
      REAL*8 VAR(LBTM,LBCV)
      INTEGER ITYP(LUK+LSOLU+2)
      INTEGER IBCSPX(LSPBC+1)
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
     &  'Phase Rel. Sat.','Inflow Phase Rel. Sat.','Phase Mass Frac.',
     &  'Inflow Phase Mass Frac.','Vol. Conc.','Inflow Vol. Conc.',
     &  'Phase Conc.','Inflow Phase Conc.','Dirichlet-Inflow',
     &  'Phase Mole Frac.','Inflow Phase Mole Frac.'/
!    &  'Aqu. Rel. Sat.','Inflow Aqu. Rel. Sat.','Aqu. Mass Frac.',
!    &  'Inflow Aqu. Mass Frac.','Vol. Conc.','Inflow Vol. Conc.',
!    &  'Aqu. Conc.','Inflow Aqu. Conc.','Dirichlet-Inflow',
!    &  'Inflow-Outflow Gas'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDBC_COMP'
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
        ISTART = 1
        VARB = 'Boundary Condition Orientation: '
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(/,A,$)') VARB(1:IVR)
        IF( INDEX(ADUM(1:),'west').NE.0 ) THEN
          IBCDX = -1
          WRITE(IWR,'(A)') 'X-Direction: West Surface'
        ELSEIF( INDEX(ADUM(1:),'east').NE.0 ) THEN
          IBCDX = 1
          WRITE(IWR,'(A)') 'X-Direction: East Surface'
        ELSEIF( INDEX(ADUM(1:),'south').NE.0 ) THEN
          IBCDX = -2
          WRITE(IWR,'(A)') 'Y-Direction: South Surface'
        ELSEIF( INDEX(ADUM(1:),'north').NE.0 ) THEN
          IBCDX = 2
          WRITE(IWR,'(A)') 'Y-Direction: North Surface'
        ELSEIF( INDEX(ADUM(1:),'bottom').NE.0 ) THEN
          IBCDX = -3
          WRITE(IWR,'(A)') 'Z-Direction: Bottom Surface'
        ELSEIF( INDEX(ADUM(1:),'top').NE.0 ) THEN
          IBCDX = 3
          WRITE(IWR,'(A)') 'Z-Direction: Top Surface'
        ELSEIF( INDEX(ADUM(1:),'file').NE.0 ) THEN
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
          NCH = INDEX(FDUM,'  ')-1
          OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='OLD',FORM='FORMATTED')
          WRITE(IWR,'(/,2A)') 'Boundary Condition Domain File: ',
     &      FDUM(1:NCH)
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
!---    Read boundary types  ---
!
        VARB = 'Boundary Condition Type'
!
!---    Energy, aqueous, gas, aqueous component, gas component
!---    and salt boundary condition types,
!       isothermal option  ---
!
        IF( ISLC(30).EQ.0 )
     &    CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(2))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(3))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(4))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(5))
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 )
     &    CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(6))
!
!---    Solute boundary condition types  ---
!
        IF( IEQC.GT.0 ) THEN
          DO 15 NSL = 1,NSOLU
!
!---        Allow for returns in input lines  ---
!
            CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
            IF( INDX.EQ.0 ) THEN
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              ISTART = 1
            ENDIF
            BDUM(NSL+LUK) = 'zero flux'
            IDFLT = 1
            VARB = 'Solute Boundary Condition Type: '
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(NSL+LUK))
   15     CONTINUE
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
!---      Gas species boundary condition types, 
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
          VARB = 'Gas Reactive Species Boundary Condition Type: '
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
!---        Gas species  ---
!
            DO 18 M = 1,NSPG
              IF( SPNMG(M).EQ.SDUM ) THEN
                MX = M + NSPL + NSPS
                GOTO 19
              ENDIF
   18       CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Aqueous or Gas Species Name: '
     &         // SDUM(1:NCH)
            CALL WRMSGS( INDX )
   19       CONTINUE
   20     CONTINUE
!        ENDIF
#endif
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
!
!---      Read energy boundary condition type ---
!
          IF( INDEX(BDUM(1)(1:),'dirichlet').NE.0 ) THEN
            ITYP(IEQT) = 1
          ELSEIF( INDEX(BDUM(1)(1:),'neumann').NE.0 ) THEN
            ITYP(IEQT) = 2
          ELSEIF( INDEX(BDUM(1)(1:),'zero flux').NE.0 ) THEN
            ITYP(IEQT) = 3
          ELSEIF( INDEX(BDUM(1)(1:),'outflow').NE.0 ) THEN
            ITYP(IEQT) = 7
          ELSEIF( INDEX(BDUM(1)(1:),'initial cond').NE.0 ) THEN
            ITYP(IEQT) = 12
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Energy Boundary Condition: '//BDUM(1)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---    Read aqueous boundary condition type ---
!
        IF( INDEX(BDUM(2)(1:),'dirichlet').NE.0 ) THEN
           ITYP(IEQW) = 1
        ELSEIF( INDEX(BDUM(2)(1:),'neumann').NE.0 ) THEN
           ITYP(IEQW) = 2
        ELSEIF( INDEX(BDUM(2)(1:),'zero flux').NE.0 ) THEN
           ITYP(IEQW) = 3
        ELSEIF( INDEX(BDUM(2)(1:),'saturated').NE.0 ) THEN
           ITYP(IEQW) = 4
        ELSEIF( INDEX(BDUM(2)(1:),'unit gradient').NE.0 ) THEN
           ITYP(IEQW) = 5
        ELSEIF( INDEX(BDUM(2)(1:),'hydraulic gradient').NE.0 ) THEN
           ITYP(IEQW) = 11
        ELSEIF( INDEX(BDUM(2)(1:),'initial cond').NE.0 ) THEN
           ITYP(IEQW) = 12
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Aqueous Boundary Condition: '//BDUM(2)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read gas boundary condition type ---
!
        IF( INDEX(BDUM(3)(1:),'dirichlet').NE.0 ) THEN
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
          CHMSG = 'Unrecognized Gas Boundary Condition: '//BDUM(3)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read aqueous component boundary condition type ---
!
! ZFZ - allows only concentration bc for aqueous components
!
        VARB = 'Aqueous Component Boundary Condition Type'
!       CALL RDCHR(ISTART,ICOMMA,NCH4,CHDUM,BDUM(4))
        IF( INDEX(BDUM(4)(1:),'mole frac').NE.0 ) THEN
          ITYP(4) = 43
          WRITE(IWR,'(2X,A)') 'Aqueous Component: ' //
     &      'Mole Fraction'
        ELSEIF( INDEX(BDUM(4)(1:),'mass frac').NE.0 ) THEN
          ITYP(4) = 36
          WRITE(IWR,'(2X,A)') 'Aqueous Component: ' //
     &      'Mass Fraction'
        ELSEIF( INDEX(BDUM(4)(1:),'relative sat').NE.0 ) THEN
          ITYP(4) = 34
          WRITE(IWR,'(2X,A)') 'Aqueous Component: ' //
     &      'Relative Saturation'
        ELSEIF( INDEX(BDUM(4)(1:),'outflow').NE.0 ) THEN
          ITYP(4) = 7
          WRITE(IWR,'(2X,A)') 'Aqueous Component: ' //
     &      'Outflow'
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Aqueous Component Boundary ' //
     &      'Condition: '//BDUM(4)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read gas component boundary condition type ---
!
! ZFZ 9/29/09 Gas component OUTFLOW and INFLOW BC can only work
!             when the gas phase BC is DIRICHLET
!
        VARB = 'Gas Component Boundary Condition Type'
!       CALL RDCHR(ISTART,ICOMMA,NCH5,CHDUM,BDUM(5))
        IF( INDEX(BDUM(3)(1:),'dirichlet').NE.0) THEN
          IF(INDEX(BDUM(5)(1:),'inflow').NE.0 .AND.
     &      INDEX(BDUM(5)(1:),'mass frac').NE.0 ) THEN
            ITYP(5) = 37
            WRITE(IWR,'(2X,A)') 'Gas Component: ' //
     &        'Mass Fraction'
          ELSEIF( INDEX(BDUM(5)(1:),'inflow').NE.0 .AND.
     &      INDEX(BDUM(5)(1:),'mole frac').NE.0 ) THEN
            ITYP(5) = 44
            WRITE(IWR,'(2X,A)') 'Gas Component: ' //
     &        'Mole Fraction'
          ELSEIF( INDEX(BDUM(5)(1:),'zero flux').NE.0 ) THEN
            ITYP(5) = 3
            WRITE(IWR,'(2X,A)') 'Gas Component: ' //
     &        'Zero Flux'
          ELSEIF( INDEX(BDUM(5)(1:),'outflow').NE.0 ) THEN
            ITYP(5) = 7
            WRITE(IWR,'(2X,A)') 'Gas Component: ' //
     &        'Outflow'
          ELSEIF( INDEX(BDUM(5)(1:),'dirichlet').NE.0 .AND.
     &      INDEX(BDUM(5)(1:),'mole frac').NE.0 ) THEN
            ITYP(5) = 43
            WRITE(IWR,'(2X,A)') 'Gas Component: ' //
     &        'Mole Fraction'
          ELSEIF( INDEX(BDUM(5)(1:),'dirichlet').NE.0 .AND.
     &      INDEX(BDUM(5)(1:),'mass frac').NE.0 ) THEN
            ITYP(5) = 36
            WRITE(IWR,'(2X,A)') 'Gas Component: ' //
     &        'Mass Fraction'
          ENDIF
        ELSEIF( INDEX(BDUM(3)(1:),'neumann').NE.0 .AND.
     &    INDEX(BDUM(5)(1:),'mole frac').NE.0 ) THEN
          ITYP(5) = 43
          WRITE(IWR,'(2X,A)') 'Gas Component: ' //
     &      'Mole Fraction'
        ELSEIF( INDEX(BDUM(3)(1:),'neumann').NE.0 .AND.
     &    INDEX(BDUM(5)(1:),'mass frac').NE.0 ) THEN
          ITYP(5) = 36
          WRITE(IWR,'(2X,A)') 'Gas Component: ' //
     &      'Mass Fraction'
        ELSEIF( INDEX(BDUM(3)(1:),'zero flux').NE.0 ) THEN
          ITYP(5) = 3
          WRITE(IWR,'(2X,A)') 'Gas Component: ' //
     &      'Zero Flux'
        ELSEIF( INDEX(BDUM(3)(1:),'initial cond').NE.0 ) THEN
          ITYP(5) = 12
          WRITE(IWR,'(2X,A)') 'Gas Component: Initial Condition'
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Gas Component Boundary ' //
     &      'Condition: '//BDUM(5)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
!
!---    Read salt boundary condition type ---
!
          IF( INDEX(BDUM(6)(1:),'inflow').NE.0 .AND.
     &      INDEX(BDUM(6)(1:),'vol').NE.0 .AND.
     &      INDEX(BDUM(6)(1:),'conc').NE.0 ) THEN
             ITYP(6) = 39
          ELSEIF( INDEX(BDUM(6)(1:),'inflow').NE.0 .AND.
     &      INDEX(BDUM(6)(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(6)(1:),'conc').NE.0 ) THEN
             ITYP(6) = 41
          ELSEIF( INDEX(BDUM(6)(1:),'inflow').NE.0 .AND.
     &      (INDEX(BDUM(6)(1:),'rel').NE.0 .OR.
     &      INDEX(BDUM(6)(1:),'sat').NE.0) ) THEN
            ITYP(6) = 35
          ELSEIF( INDEX(BDUM(6)(1:),'inflow').NE.0 .AND.
     &      (INDEX(BDUM(6)(1:),'mass').NE.0 .OR.
     &      INDEX(BDUM(6)(1:),'frac').NE.0) ) THEN
            ITYP(6) = 37
          ELSEIF( INDEX(BDUM(6)(1:),'vol').NE.0 .AND.
     &      INDEX(BDUM(6)(1:),'conc').NE.0 ) THEN
             ITYP(6) = 38
          ELSEIF( INDEX(BDUM(6)(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(6)(1:),'conc').NE.0 ) THEN
             ITYP(6) = 40
          ELSEIF( INDEX(BDUM(6)(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(6)(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(6)(1:),'frac').NE.0 ) THEN
             ITYP(6) = 36
          ELSEIF( INDEX(BDUM(6)(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(6)(1:),'rel').NE.0 .AND.
     &      INDEX(BDUM(6)(1:),'sat').NE.0 ) THEN
            ITYP(6) = 34
          ELSEIF( INDEX(BDUM(6)(1:),'zero flux').NE.0 ) THEN
            ITYP(6) = 3
          ELSEIF( INDEX(BDUM(6)(1:),'outflow').NE.0 ) THEN
            ITYP(6) = 7
          ELSEIF( INDEX(BDUM(6)(1:),'initial cond').NE.0 ) THEN
            ITYP(6) = 12
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Salt Boundary Condition: '//BDUM(6)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---    Read solute boundary condition type(s) ---
!
        IF( IEQC.GT.0 ) THEN
          DO 25 NSL = 1,NSOLU
            IF( INDEX(BDUM(NSL+LUK)(1:),
     &        'inflow-outflow aqu').NE.0 ) THEN
              ITYP(NSL+LUK) = 23
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),
     &        'inflow-outflow gas').NE.0 ) THEN
              ITYP(NSL+LUK) = 43
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),
     &        'inflow-outflow').NE.0 ) THEN
              ITYP(NSL+LUK) = 19
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'inflow aqu').NE.0 ) THEN
               ITYP(NSL+LUK) = 14
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'inflow gas').NE.0 ) THEN
               ITYP(NSL+LUK) = 15
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'inflow').NE.0 ) THEN
               ITYP(NSL+LUK) = 13
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'outflow').NE.0 ) THEN
               ITYP(NSL+LUK) = 7
            ELSEIF(INDEX(BDUM(NSL+LUK)(1:),
     &        'volumetric conc').NE.0 ) THEN
               ITYP(NSL+LUK) = 1
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'aqueous conc').NE.0 )THEN
               ITYP(NSL+LUK) = 8
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'gas conc').NE.0 ) THEN
               ITYP(NSL+LUK) = 9
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
#ifdef ecke
!
!---    Read aqueous reactive species boundary condition type(s) ---
!
        IF( ISLC(40).EQ.1 ) THEN
          IF( INDEX(BDUM(NSOLU+LUK+1)(1:),
     &      'inflow-outflow aqu').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 23
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),
     &      'inflow-outflow gas').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 43
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),
     &      'inflow-outflow').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 19
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'outflow').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 7
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'initial co').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 12
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'inflow aqu').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 14
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'inflow gas').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 15
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'inflow').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 13
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),
     &      'volumetric conc').NE.0 ) THEN
             ITYP(NSOLU+LUK+1) = 1
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'aqueous conc').NE.0 )THEN
            ITYP(NSOLU+LUK+1) = 8
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'gas conc').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 9
          ELSEIF( INDEX(BDUM(NSOLU+LUK+1)(1:),'zero flux').NE.0 ) THEN
            ITYP(NSOLU+LUK+1) = 3
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Reactive Species Boundary Condition: '
     &        //BDUM(NSOLU+LUK+1)
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Read gas reactive species boundary condition type(s) ---
!
          IF( INDEX(BDUM(NSOLU+LUK+2)(1:),'volumetric conc').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 1
            IF( ITYP(NSOLU+LUK+1).NE.1 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.8 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.9 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.12 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'aqueous con').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 8
            IF( ITYP(NSOLU+LUK+1).NE.1 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.8 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.9 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.12 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'gas conc').NE.0 )THEN
            ITYP(NSOLU+LUK+2) = 9
            IF( ITYP(NSOLU+LUK+1).NE.8 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.9 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.12 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'zero flux').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 3
            IF( ITYP(NSOLU+LUK+1).NE.3 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),
     &      'inflow-outflow aqu').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 23
            IF( ITYP(NSOLU+LUK+1).NE.19 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.23 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.43 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),
     &      'inflow-outflow gas').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 43
            IF( ITYP(NSOLU+LUK+1).NE.19 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.23 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.43 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),
     &      'inflow-outflow').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 19
            IF( ITYP(NSOLU+LUK+1).NE.19 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.23 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.43 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'outflow').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 7
            IF( ITYP(NSOLU+LUK+1).NE.7 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'initial co').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 12
            IF( ITYP(NSOLU+LUK+1).NE.1 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.8 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.9 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.12 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'inflow').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 13
            IF( ITYP(NSOLU+LUK+1).NE.13 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.14 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.15 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'inflow aqu').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 14
            IF( ITYP(NSOLU+LUK+1).NE.13 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.14 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.15 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
     &          'Boundary Condition: '//BDUM(NSOLU+LUK+1)//
     &          ' and '//BDUM(NSOLU+LUK+2)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),'inflow gas').NE.0 ) THEN
            ITYP(NSOLU+LUK+2) = 15
            IF( ITYP(NSOLU+LUK+1).NE.13 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.14 .AND. 
     &        ITYP(NSOLU+LUK+1).NE.15 ) THEN
              INDX = 4
              CHMSG = 'Mixed Aqueous-Gas Reactive Species '//
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
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 )
     &    WRITE(IWR,'(2X,2A)') 'Energy: ',CHTYP(ITYP(IEQT))
        WRITE(IWR,'(2X,2A)') 'Aqueous: ',CHTYP(ITYP(IEQW))
        WRITE(IWR,'(2X,2A)') 'Gas: ',CHTYP(ITYP(3))
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 )
     &    WRITE(IWR,'(2X,2A)') 'Salt: ',CHTYP(ITYP(6))
        IF( IEQC.GT.0 ) THEN
          DO 30 NSL = 1,NSOLU
            IDB = INDEX( SOLUT(NSL)(1:),'  ') - 1
            WRITE(IWR,'(2X,2A)') SOLUT(NSL)(1:IDB),CHTYP(ITYP(NSL+LUK))
   30     CONTINUE
        ENDIF
#ifdef ecke
!
!---    Write aqueous and gas species boundary condition type(s) ---
!
        IF( ISLC(40).EQ.1 ) THEN
          WRITE(IWR,'(2X,3A)') 'Aqueous: ',
     &      CHTYP(ITYP(NSOLU+LUK+1)),' Reactive Species'
          WRITE(IWR,'(2X,3A)') 'Gas: ',
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
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Read energy boundary condition variables ---
!
            IF( ITYP(IEQT).EQ.1 ) THEN
              VARB = 'Temperature'
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            ELSEIF( ITYP(IEQT).EQ.2 ) THEN
              VARB = 'Energy Flux'
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNKG = 1
              IUNS = -3
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            ELSEIF( ITYP(IEQT).EQ.3 ) THEN
              VARB = 'Zero Flux'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            ELSEIF( ITYP(IEQT).EQ.7 ) THEN
              VARB = 'Outflow'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            ELSEIF( ITYP(IEQT).EQ.12 ) THEN
              VARB = 'Initial Condition'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            ELSE
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            ENDIF
          ENDIF
!
!---      Read aqueous boundary condition variables ---
!
          IF( ITYP(IEQW).EQ.1 ) THEN
            VARB = 'Aqueous Pressure'
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
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
          ELSEIF( ITYP(IEQW).EQ.2 ) THEN
            VARB = 'Volumetric Aqueous Flux'
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
              INDX = 0
              IUNM = 1
              IUNS = -1
              CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
          ELSEIF( ITYP(IEQW).EQ.3 ) THEN
            VARB = 'Aqueous Pressure'
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
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
          ELSEIF( ITYP(IEQW).EQ.11 ) THEN
            VARB = 'Base Aqueous Pressure'
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
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
          ELSE
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          ENDIF
          IF( ITYP(3).EQ.1 ) THEN
            VARB = 'Gas Pressure'
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            ISTART = ISX
            ICOMMA = ICX
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
          ELSEIF( ITYP(3).EQ.2 ) THEN
            VARB = 'Volumetric Gas Flux'
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            ISTART = ISX
            ICOMMA = ICX
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
          ELSEIF( ITYP(3).EQ.3 ) THEN
            VARB = 'Gas Pressure'
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            ISTART = ISX
            ICOMMA = ICX
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
          ELSEIF( ITYP(3).EQ.11 ) THEN
            VARB = 'Base Gas Pressure'
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            ISTART = ISX
            ICOMMA = ICX
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
!---      Read water-vapor relative saturation ---
!
          VARB = 'Water-Vapor Relative Saturation, '
          ISX = ISTART
          ICX = ICOMMA
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
          ISTART = ISX
          ICOMMA = ICX
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',VAR(NTM,5)
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Read salt boundary condition variables ---
!
            IF( ITYP(6).EQ.38 .OR. ITYP(6).EQ.39 ) THEN
              VARB = 'Salt Volumetric Concentration, '
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),
     &          UNTS(1:NCH),': ',VAR(NTM,6)
              INDX = 0
              IUNKG = 1
              IUNM = -3
              CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
            ELSEIF( ITYP(6).EQ.40 .OR. ITYP(6).EQ.41 ) THEN
              VARB = 'Dissolved-Salt Aqueous Concentration, '
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),
     &          UNTS(1:NCH),': ',VAR(NTM,6)
              INDX = 0
              IUNKG = 1
              IUNM = -3
              CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
            ELSEIF( ITYP(6).EQ.34 .OR. ITYP(6).EQ.35 ) THEN
              VARB = 'Dissolved-Salt Aqueous Relative Saturation, '
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,6)
            ELSEIF( ITYP(6).EQ.36 .OR. ITYP(6).EQ.37 ) THEN
              VARB = 'Dissolved-Salt Aqueous Mass Fraction, '
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,6)
            ELSE
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            ENDIF
          ENDIF
!
!---      Aqueous and gas component boundary conditions  ---
!
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
!
!---      Aqueous component boundary conditions  ---
!
          DO 125 IGC = 1,NGC
            NCH = INDEX( GCNM(IGC),'  ') - 1
            VARB = 'Aqueous ' // GCNM(IGC)(1:NCH) //
     &        '  Concentration'
            IF( ITYP(4).EQ.43 .OR. ITYP(4).EQ.44 ) THEN
              NCH = INDEX( GCNM(IGC),'  ') - 1
              VARB = 'Aqueous ' // GCNM(IGC)(1:NCH) //
     &          '  Mole Fraction'
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              IDFLT = 1
              VAR(NTM,6+IGC) = 0.D+0
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6+IGC))
            ELSEIF( ITYP(4).EQ.36 .OR. ITYP(4).EQ.37 ) THEN
              NCH = INDEX( GCNM(IGC),'  ') - 1
              VARB = 'Aqueous ' // GCNM(IGC)(1:NCH) //
     &          '  Mass Fraction'
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              IDFLT = 1
              VAR(NTM,6+IGC) = 0.D+0
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6+IGC))
            ELSEIF( ITYP(4).EQ.34 .OR. ITYP(4).EQ.35 ) THEN
              NCH = INDEX( GCNM(IGC),'  ') - 1
              VARB = 'Aqueous ' // GCNM(IGC)(1:NCH) //
     &          '  Relative Saturation'
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              IDFLT = 1
              VAR(NTM,6+IGC) = 0.D+0
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6+IGC))
            ELSE
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
            ENDIF
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),': '
            WRITE(IWR,'(1PE11.4)') VAR(NTM,6+IGC)
  125     CONTINUE
!
!---      Gas component boundary conditions  ---
!
          DO 130 IGC = 1,NGC
            NCH = INDEX( GCNM(IGC),'  ') - 1
            VARB = 'Gas ' // GCNM(IGC)(1:NCH) //
     &        '  Concentration'
            IF( ITYP(5).EQ.43 .OR. ITYP(5).EQ.44 ) THEN
              NCH = INDEX( GCNM(IGC),'  ') - 1
              VARB = 'Gas ' // GCNM(IGC)(1:NCH) //
     &          '  Mole Fraction'
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              IDFLT = 1
              VAR(NTM,6+LNGC+IGC) = 0.D+0
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6+LNGC+IGC))
            ELSEIF( ITYP(5).EQ.36 .OR. ITYP(5).EQ.37 ) THEN
              NCH = INDEX( GCNM(IGC),'  ') - 1
              VARB = 'Gas ' // GCNM(IGC)(1:NCH) //
     &          '  Mass Fraction'
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              IDFLT = 1
              VAR(NTM,6+LNGC+IGC) = 0.D+0
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6+LNGC+IGC))
            ELSE
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
            ENDIF
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),': '
            WRITE(IWR,'(1PE11.4)') VAR(NTM,6+LNGC+IGC)
  130     CONTINUE
!
!---        Read solute transport boundary condition variables ---
!
          IF( IEQC.GT.0 ) THEN
            DO 50 NSL = 1,NSOLU
              IF( ITYP(NSL+LUK).EQ.1 ) THEN
                VARB = 'Volumetric Concentration'
                ISX = ISTART
                ICX = ICOMMA
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
                ISTART = ISX
                ICOMMA = ICX
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+LBCU))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL),VARB(1:IVR),
     &            ', ',UNTS(1:NCH),': ',VAR(NTM,NSL+LBCU)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VAR(NTM,NSL+LBCU),INDX)
              ELSEIF( ITYP(NSL+LUK).EQ.8 ) THEN
                VARB = 'Aqueous-Phase Volumetric Concentration'
                ISX = ISTART
                ICX = ICOMMA
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
                ISTART = ISX
                ICOMMA = ICX
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+LBCU))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL),VARB(1:IVR),
     &            ', ',UNTS(1:NCH),': ',VAR(NTM,NSL+LBCU)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VAR(NTM,NSL+LBCU),INDX)
              ELSEIF( ITYP(NSL+LUK).EQ.9 ) THEN
                VARB = 'Gas-Phase Volumetric Concentration'
                ISX = ISTART
                ICX = ICOMMA
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
                ISTART = ISX
                ICOMMA = ICX
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
#ifdef ecke
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
                  VARB = 'Aqueous-Phase Concentration, '
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
                  VARB = 'Gas-Phase Concentration, '
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
#endif
!
!---  Check for nonascending boundary condition times  ---
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
              NBCL = NBCL + 1
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
!---  Check for boundary values applied to inactive nodes  ---
!
              IF( IXP(N).EQ.0 ) THEN
                INDX = 7
                IMSG = N
                CHMSG = 'Boundary Condition Applied to an Inactive Node'
                CALL WRMSGS( INDX )
              ENDIF
!
!---  Check for boundary values applied to interior surfaces  ---
!
              IERR = 0
              IF( IBCDX.EQ.-1 .AND. I.NE.1) THEN
                IF( IXP(N-1).NE.0 ) IERR = 1
              ELSEIF( IBCDX.EQ.1 .AND. I.NE.IFLD) THEN
                IF( IXP(N+1).NE.0 ) IERR = 1
              ELSEIF( IBCDX.EQ.-2 .AND. J.NE.1) THEN
                IF( IXP(N-IFLD).NE.0 ) IERR = 1
              ELSEIF( IBCDX.EQ.2 .AND. J.NE.JFLD) THEN
                IF( IXP(N+IFLD).NE.0 ) IERR = 1
              ELSEIF( IBCDX.EQ.-3 .AND. K.NE.1) THEN
                IF( IXP(N-IJFLD).NE.0 ) IERR = 1
              ELSEIF( IBCDX.EQ.3 .AND. K.NE.KFLD) THEN
                IF( IXP(N+IJFLD).NE.0 ) IERR = 1
              ENDIF
              IF( IERR.EQ.1 ) THEN
                INDX = 7
                IMSG = NBC
                CHMSG = 'Boundary Cond. Applied to an Interior Surface'
     &            //': Boundary Number'
                CALL WRMSGS( INDX )
              ENDIF
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
              IF( ISLC(30).EQ.0 ) IBCT(IEQT,NBC) = ITYP(IEQT)
              IBCT(IEQW,NBC) = ITYP(IEQW)
              IBCT(4,NBC) = ITYP(3)
              IBCT(5,NBC) = ITYP(4)
              IBCT(6,NBC) = ITYP(5)
              IF( ISLC(32).EQ.0 ) IBCT(IEQS,NBC) = ITYP(6)
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
              IF( ITYP(3).EQ.11 ) THEN
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
!---  End of RDBC_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDGCP_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-NComponent-NaCl-Energy Mode (STOMP-WNSE)
!     Read gas component names.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      CHARACTER*64 ADUM,RDUM,UNTS
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDGCP_COMP'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Gas Components Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of gas components  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Gas Components'
      CALL RDINT( ISTART,ICOMMA,CHDUM,NGC )
      IF( NGC.GT.LNGC ) THEN
        INDX = 5
        CHMSG = 'Number of Gas Components > Parameter LNGC'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Loop over number of gas components  ---
!
      DO 500 IGC = 1,NGC
!
!---    Read new input line ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Gas Component Name'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GCNM(IGC))
        DO 100 M = 1,IGC-1
          IF( GCNM(M).EQ.GCNM(IGC) ) THEN
            INDX = 4
            CHMSG = 'Duplicate Gas Component Name: ' // GCNM(IGC)
            CALL WRMSGS( INDX )
          ENDIF
  100   CONTINUE
        WRITE (IWR,'(/,2A)') 'Gas Component Name: ',GCNM(IGC)
!
!---  Go to next gas component  ---
!
  500 CONTINUE
!
!---  End of RDGCP_COMP group ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIC_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Read input file for initial conditions information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      CHARACTER*64 ADUM,BDUM,CDUM,FDUM,FMDUM,UNTS
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
      SUB_LOG(ISUB_LOG) = '/RDIC_COMP'
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
      CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
      IF( IEO.EQ.2 ) GOTO 10
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
   10 CONTINUE
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
        ELSEIF( INDEX( ADUM(1:),'component pres' ).NE.0 ) THEN
          VARB = 'Initial Component Pressure'
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
        ELSEIF( INDEX( ADUM(1:),'gas' ).NE.0 ) THEN
          IF( INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
            VARB = 'Component Name: '
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            VARB = 'Initial Component Gas Mass Fraction'
            CDUM = BDUM(1:NCHB) // ' Gas Mass Fraction'
            NCHC = INDEX( CDUM(1:),'  ' ) - 1
          ELSEIF( INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
            VARB = 'Component Name: '
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            VARB = 'Initial Component Gas Mole Fraction'
            CDUM = BDUM(1:NCHB) // ' Gas Mole Fraction'
            NCHC = INDEX( CDUM(1:),'  ' ) - 1
          ENDIF
        ELSEIF( INDEX( ADUM(1:),'aqueous' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'specie' ).EQ.0 .AND.
     &      INDEX( ADUM(1:),'salt' ).EQ.0 ) THEN
          IF( INDEX (ADUM(1:),'mass' ).NE.0 ) THEN
            VARB = 'Component Name: '
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            VARB = 'Initial Component Aqueous Mass Fraction'
            CDUM = BDUM(1:NCHB) // ' Aqueous Mass Fraction'
            NCHC = INDEX( CDUM(1:),'  ' ) - 1
          ELSEIF( INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
            VARB = 'Component Name: '
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            VARB = 'Initial Component Aqueous Mole Fraction'
            CDUM = BDUM(1:NCHB) // ' Aqueous Mole Fraction'
            NCHC = INDEX( CDUM(1:),'  ' ) - 1
          ENDIF
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
          IVAR = 1
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
#endif
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
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( T,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( T,VAR,ADDER,IDOM,INDX )
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
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PG,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( PG,VAR,ADDER,IDOM,INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'component pres').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( PVA,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( PVA,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( PVA,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PVA,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( PVA,VAR,ADDER,IDOM,INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &    INDEX(ADUM(1:),'sat').NE.0 ) THEN
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
              CALL RDINBS( SGT,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( SGT,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( SGT,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( SGT,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( SGT,VAR,ADDER,IDOM,INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'trapped').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( SGT,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( SGT,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( SGT,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( SGT,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( SGT,VAR,ADDER,IDOM,INDX )
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
            DO 206 N = 1,NFLD
              IF( IXP(N).LE.0 ) GOTO 206
              ICBRN(N) = IVAR
  206       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( TMS,VAR(1),ADDER,IROCK,INDX )
            DO 208 N = 1,NFLD
              IF( IXP(N).LE.0 ) GOTO 208
              IF( IZ(N).EQ.IROCK ) ICBRN(N) = IVAR
  208       CONTINUE
          ELSE
            CALL RDINIS( TMS,VAR,ADDER,IDOM,INDX )
            DO 210 K = IDOM(5),IDOM(6)
            DO 210 J = IDOM(3),IDOM(4)
            DO 210 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).LE.0 ) GOTO 210
              ICBRN(N) = IVAR
  210       CONTINUE
          ENDIF
        ELSEIF( (INDEX(ADUM(1:),'gas').NE.0) .AND.
     &    (INDEX(ADUM(1:),'mass').NE.0) ) THEN
          IGCX = 0
          DO 310 IGC = 1,NGC
            IDB = INDEX(GCNM(IGC)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.GCNM(IGC)(1:IDB) ) IGCX = IGC
  310     CONTINUE
          IF( IGCX.EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Gas Component Name: '//BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNK = 1
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBGC( XGC,ADDER,UNTS,INDX,IGCX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAGC( XGC,ADDER,UNTS,INDX,IGCX )
            ELSE
              CALL RDINFGC( XGC,VAR,ADDER,UNTS,INDX,IGCX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZGC( XGC,VAR(1),ADDER,IROCK,INDX,IGCX )
          ELSE
            CALL RDINIGC( XGC,VAR,ADDER,IDOM,INDX,IGCX )
          ENDIF
        ELSEIF( (INDEX(ADUM(1:),'gas').NE.0) .AND.
     &    (INDEX(ADUM(1:),'mole').NE.0) ) THEN
          IGCX = 0
          DO 320 IGC = 1,NGC
            IDB = INDEX(GCNM(IGC)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.GCNM(IGC)(1:IDB) ) IGCX = IGC
  320     CONTINUE
          IF( IGCX.EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Gas Component Name: '//BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNK = 1
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBGC( XMGC,ADDER,UNTS,INDX,IGCX )
              XMGC = XMGC
              XMVGC = XMGC
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAGC( XMGC,ADDER,UNTS,INDX,IGCX )
              XMGC = XMGC
              XMVGC = XMGC
            ELSE
              CALL RDINFGC( XMGC,VAR,ADDER,UNTS,INDX,IGCX )
              XMGC = XMGC
              XMVGC = XMGC
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZGC( XMGC,VAR(1),ADDER,IROCK,INDX,IGCX )
            XMGC = XMGC
            XMVGC = XMGC
          ELSE
            CALL RDINIGC( XMGC,VAR,ADDER,IDOM,INDX,IGCX )
            XMGC = XMGC
            XMVGC = XMGC
          ENDIF
        ELSEIF( (INDEX(ADUM(1:),'aqueous').NE.0) .AND.
     &    (INDEX(ADUM(1:),'mass').NE.0) ) THEN
          IGCX = 0
          DO 330 IGC = 1,NGC
            IDB = INDEX(GCNM(IGC)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.GCNM(IGC)(1:IDB) ) IGCX = IGC
  330     CONTINUE
          IF( IGCX.EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Gas Component Name: '//BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNK = 1
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBGC( XLC,ADDER,UNTS,INDX,IGCX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAGC( XLC,ADDER,UNTS,INDX,IGCX )
            ELSE
              CALL RDINFGC( XLC,VAR,ADDER,UNTS,INDX,IGCX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZGC( XLC,VAR(1),ADDER,IROCK,INDX,IGCX )
          ELSE
            CALL RDINIGC( XLC,VAR,ADDER,IDOM,INDX,IGCX )
            DO K = IDOM(5),IDOM(6)
            DO J = IDOM(3),IDOM(4)
            DO I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).LE.0 ) CYCLE
              ICCOMP(IGCX,N) = 3
            ENDDO
            ENDDO
            ENDDO
          ENDIF
        ELSEIF( (INDEX(ADUM(1:),'aqueous').NE.0) .AND.
     &    (INDEX(ADUM(1:),'mole').NE.0) ) THEN
          IGCX = 0
          DO 340 IGC = 1,NGC
            IDB = INDEX(GCNM(IGC)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.GCNM(IGC)(1:IDB) ) IGCX = IGC
  340     CONTINUE
          IF( IGCX.EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Gas Component Name: '//BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNK = 1
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBGC( XMLC,ADDER,UNTS,INDX,IGCX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAGC( XMLC,ADDER,UNTS,INDX,IGCX )
            ELSE
              CALL RDINFGC( XMLC,VAR,ADDER,UNTS,INDX,IGCX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZGC( XMLC,VAR(1),ADDER,IROCK,INDX,IGCX )
          ELSE
            CALL RDINIGC( XMLC,VAR,ADDER,IDOM,INDX,IGCX )
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
#endif
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIC_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINPT_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Read input file cards.
!     Direct control to card reader subroutines.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/RDINPT_COMP'
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
        CALL RDSP_COMP
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
        GOTO 1100
      ELSE
        GOTO 900
      ENDIF
  910 CONTINUE
      INDX = 4
      CHMSG = 'Missing Aqueous Relative Permeability Card'
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
        GOTO 1300
      ELSE
        GOTO 1200
      ENDIF
 1250 CONTINUE
!
!---  Nonisothermal simulations  ---
!
      IF( ISLC(30).EQ.0 ) THEN
        INDX = 4
        CHMSG = 'Missing Thermal Properties Card'
        CALL WRMSGS( INDX )
      ELSE
        REWIND(IRD)
      ENDIF
 1300 CONTINUE
!
!---  Search input file for gas component properties card  ---
!
 1309 READ(IRD,'(A)', END=1310) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 1309
      CALL LCASE( CHDUM )
      IF( CHDUM(1:14).EQ.'~gas component' ) THEN
        CALL RDGCP_COMP
        REWIND(IRD)
        GOTO 1400
      ELSE
        GOTO 1300
      ENDIF
 1310 CONTINUE
      INDX = 18
      CHMSG = 'Missing Gas Component Properties Card'
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
        CALL RDTF_COMP
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
        CALL RDTP_COMP
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
!---  Search input file for gas species card  ---
!
 2259 READ(IRD,'(A)', END=2260) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2259
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'gas').NE.0 .AND.
     &  INDEX(CHDUM(2:),'specie').NE.0) THEN
        CALL RDGSSP
        REWIND(IRD)
        GOTO 2300
      ELSE
        GOTO 2259
      ENDIF
 2260 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Gas Species Card'
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
!---  Search input file for salt transport card --
!
 3200 CONTINUE
 3209 READ(IRD,'(A)', END=3210) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3209
      CALL LCASE( CHDUM )
      IF( CHDUM(1:15).EQ.'~salt transport' ) THEN
        CALL RDST_COMP
        REWIND(IRD)
        GOTO 3300
      ELSE
        GOTO 3200
      ENDIF
 3210 CONTINUE
      INDX = 4
      CHMSG = 'Missing Salt Transport Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for source card --
!
 3300 CONTINUE
 3309 READ(IRD,'(A)', END= 3310) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3309
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'source').NE.0 ) THEN
        CALL RDSR_COMP
        REWIND(IRD)
        GOTO 3400
      ELSE
        GOTO 3300
      ENDIF
 3310 CONTINUE
      INDX = 1
      CHMSG = 'Missing Source Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for initial conditions card --
!
 3400 CONTINUE
 3409 READ(IRD,'(A)', END=3410) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3409
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'initial').NE.0 ) THEN
        CALL RDIC_COMP
        REWIND(IRD)
        GOTO 3500
      ELSE
        GOTO 3400
      ENDIF
 3410 CONTINUE
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
 3500 CONTINUE
 3509 READ(IRD,'(A)', END=3510) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3509
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'boundary').NE.0 .AND.
     &  INDEX(CHDUM(2:),'geomech').EQ.0 ) THEN
        CALL RDBC_COMP
        REWIND(IRD)
        GOTO 3600
      ELSE
        GOTO 3500
      ENDIF
 3510 CONTINUE
      INDX = 1
      CHMSG = 'Missing Boundary Conditions Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for output control card --
!
 3600 CONTINUE
 3609 READ(IRD,'(A)', END=3610) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3609
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'output').NE.0 ) THEN
        CALL RDOU_COMP
        REWIND(IRD)
        GOTO 3700
      ELSE
        GOTO 3600
      ENDIF
 3610 CONTINUE
      INDX = 1
      CHMSG = 'Missing Output Control Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for surface flux card --
!
 3700 CONTINUE
 3709 READ(IRD,'(A)', END=3710) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3709
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'surface').NE.0 ) THEN
        CALL RDSF_COMP
        REWIND(IRD)
        GOTO 3800
      ELSE
        GOTO 3700
      ENDIF
 3710 CONTINUE
      INDX = 1
      CHMSG = 'Missing Surface Flux Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
 3800 CONTINUE
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
!---  End of RDINPT_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDOU_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Read input file for output information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE OUTPU
      USE GRID
      USE FILES
      USE FDVS
      USE FDVGC
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
      CHARACTER*64 ADUM,UNTS,SOLNM,SPNM,GCNMX
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
      SUB_LOG(ISUB_LOG) = '/RDOU_COMP'
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
!---      Gas species  ---
!
          DO 136 M = 1,NSPG
            NSP = M+NSPL+NSPS
            IF( SPNM.EQ.SPNMG(M) ) GOTO 150
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
        IREFGC(NV) = 0
        IF( INDEX(ADUM(1:),'aqueous pressure').NE.0 ) THEN
          IREF(NV) = 1
        ELSEIF( INDEX(ADUM(1:),'gas pressure').NE.0 ) THEN
          IREF(NV) = 2
        ELSEIF( INDEX(ADUM(1:),'source-well pres').NE.0 ) THEN
          IREF(NV) = 285
        ELSEIF( INDEX(ADUM(1:),'source-well temp').NE.0 ) THEN
          IREF(NV) = 286
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
        ELSEIF( INDEX(ADUM(1:),'apparent aqueous sat').NE.0 ) THEN
          IREF(NV) = 9
        ELSEIF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
          IREF(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'gas saturation').NE.0 ) THEN
          IREF(NV) = 12
        ELSEIF( INDEX(ADUM(1:),'salt saturation').NE.0 ) THEN
          IREF(NV) = 264
        ELSEIF( INDEX(ADUM(1:),'aqueous moisture cont').NE.0 ) THEN
          IREF(NV) = 15
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mole frac').NE.0 ) THEN
          CHREF(205) = 'XMLS'
          IREF(NV) = 205
        ELSEIF( INDEX(ADUM(1:),'aqueous mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous mole frac').NE.0 ) THEN
          IREF(NV) = 204
          CHREF(204) = 'XMLC'
          DO IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IREFGC(NV) = IGC
          ENDDO
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
           CHMSG = 'Unrecognized Aqueous Mole Fraction Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          END IF
        ELSEIF( INDEX(ADUM(1:),'effective trapped gas').NE.0 ) THEN
          IREF(NV) = 19
        ELSEIF( INDEX(ADUM(1:),'diffusive porosity').NE.0 ) THEN
          IREF(NV) = 20
        ELSEIF( INDEX(ADUM(1:),'h2o gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o mass frac').NE.0 ) THEN
          IREF(NV) = 21
        ELSEIF( INDEX(ADUM(1:),'gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas mass frac').NE.0 ) THEN
          IREF(NV) = 22
          CHREF(22) = 'XGC'
          DO IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )
     &        IREFGC(NV) = IGC
          ENDDO
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Gas Mass Fraction Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'h2o aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous h2o mass frac').NE.0 ) THEN
          IREF(NV) = 24
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mass').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mass').NE.0 ) THEN
          IREF(NV) = 110
        ELSEIF( INDEX(ADUM(1:),'aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous mass frac').NE.0 ) THEN
          IREF(NV) = 25
          CHREF(25) = 'XLC'
          DO IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )
     &        IREFGC(NV) = IGC
          ENDDO
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Aqu. Mass Fraction Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
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
        ELSEIF( INDEX(ADUM(1:),'total h2o mass').NE.0 ) THEN
          IREF(NV) = 37
        ELSEIF ( INDEX(ADUM(1:),'total salt mass').NE.0 ) THEN
          IREF(NV) = 50
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 38
          CHREF(38) = 'TMC'
          DO 146 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )
     &        IREFGC(NV) = IGC
  146     CONTINUE
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Total Mass Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'h2o mass source int').NE.0 ) THEN
          IREF(NV) = 40
        ELSEIF( INDEX(ADUM(1:),'mass source int').NE.0 ) THEN
          IREF(NV) = 41
          CHREF(41) = 'SRIC'
          DO 148 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )
     &        IREFGC(NV) = IGC
  148     CONTINUE
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Mass Source Integral Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF ( INDEX(ADUM(1:),'salt volumetric conc').NE.0 ) THEN
          IREF(NV) = 47
        ELSEIF ( INDEX(ADUM(1:),'salt aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt conc').NE.0 ) THEN
          IREF(NV) = 48
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
        ELSEIF( INDEX(ADUM(1:),'gas mole').NE.0 ) THEN
          IREF(NV) = 71
          CHREF(71) = 'XMGC'
          DO 155 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IREFGC(NV) = IGC
  155     CONTINUE
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Gas Mole Fraction Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'h2o gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o conc').NE.0 ) THEN
          IREF(NV) = 73
        ELSEIF( INDEX(ADUM(1:),'gas conc').NE.0 ) THEN
          IREF(NV) = 74
          CHREF(74) = 'CGC'
          DO 152 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IREFGC(NV) = IGC
  152     CONTINUE
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Gas Concentration Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'h2o aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous h2o conc').NE.0 ) THEN
          IREF(NV) = 76
        ELSEIF( INDEX(ADUM(1:),'aqueous conc').NE.0 .AND.
     &    INDEX(ADUM(1:),'specie').EQ.0 ) THEN
          IREF(NV) = 77
          CHREF(77) = 'CLC'
          DO 154 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IREFGC(NV) = IGC
  154     CONTINUE
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Aqu. Concentration Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
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
        ELSEIF( INDEX(ADUM(1:),'osmotic pressure').NE.0 ) THEN
          IREF(NV) = 101
        ELSEIF( INDEX(ADUM(1:),'osmotic eff').NE.0 ) THEN
          IREF(NV) = 102
        ELSEIF( INDEX(ADUM(1:),'water').NE.0 .AND.
     &    (INDEX(ADUM(1:),'gas').NE.0 .OR. INDEX(ADUM(1:),'vapor').NE.0)
     &    .AND. INDEX(ADUM(1:),'partial').NE.0
     &    .AND. INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IREF(NV) = 128
        ELSEIF( (INDEX(ADUM(1:),'gas').NE.0 
     &      .OR. INDEX(ADUM(1:),'vapor').NE.0)
     &     .AND. INDEX(ADUM(1:),'partial').NE.0
     &     .AND. INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IREF(NV) = 129
          CHREF(129) = 'PVC'
          DO 129 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IREFGC(NV) = IGC
  129     CONTINUE
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas-aqueous scaling').NE.0 ) THEN
          IREF(NV) = 131
        ELSEIF( INDEX(ADUM(1:),'h2o mass source rate').NE.0 ) THEN
          IREF(NV) = 140
        ELSEIF( INDEX(ADUM(1:),'mass source rate').NE.0 ) THEN
          IREF(NV) = 141
          CHREF(141) = 'SRCC'
          DO 156 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IREFGC(NV) = IGC
  156     CONTINUE
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Mass Source Rate Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'salt mass source rate').NE.0 ) THEN
          IREF(NV) = 147
        ELSEIF( INDEX(ADUM(1:),'salt mass source int').NE.0 ) THEN
          IREF(NV) = 148
        ELSEIF( INDEX(ADUM(1:),'aqueous viscosity').NE.0 ) THEN
          IREF(NV) = 176
        ELSEIF( INDEX(ADUM(1:),'gas viscosity').NE.0 ) THEN
          IREF(NV) = 289
        ELSEIF( INDEX(ADUM(1:),'integrated h2o mass').NE.0 ) THEN
          IREF(NV) = 191
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND. 
     &    INDEX(ADUM(1:),'species').EQ.0 ) THEN
          IREF(NV) = 192
          CHREF(192) = 'IMC'
          DO 162 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IREFGC(NV) = IGC
  162     CONTINUE
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Integrated Mass Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'integrated h2o aqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated aqueous h2o').NE.0 ) THEN
          IREF(NV) = 194
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
          IREF(NV) = 195
          CHREF(195) = 'IMLC'
          DO 158 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IREFGC(NV) = IGC
  158     CONTINUE
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Integrated Aqueous Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'integrated h2o gas').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated gas h2o').NE.0 ) THEN
          IREF(NV) = 197
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IREF(NV) = 198
          CHREF(198) = 'IMGC'
          DO 160 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IREFGC(NV) = IGC
  160     CONTINUE
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Integrated Aqueous Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
          IREF(NV) = 198
        ELSEIF( INDEX(ADUM(1:),'x intrinsic perm').NE.0 ) THEN
          IREF(NV) = 247
        ELSEIF( INDEX(ADUM(1:),'y intrinsic perm').NE.0 ) THEN
          IREF(NV) = 248
        ELSEIF( INDEX(ADUM(1:),'z intrinsic perm').NE.0 ) THEN
          IREF(NV) = 249
        ELSEIF( INDEX(ADUM(1:),'integrated h2o mass sour').NE.0 ) THEN
          IREF(NV) = 261
        ELSEIF( INDEX(ADUM(1:),'well-node').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IREF(NV) = 275
        ELSEIF( INDEX(ADUM(1:),'well-node').NE.0 .AND.
     &    INDEX(ADUM(1:),'temp').NE.0 ) THEN
          IREF(NV) = 276
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'node').NE.0 .AND.
     &    INDEX(ADUM(1:),'centroid').NE.0 ) THEN
          IREF(NV) = 291
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'node').NE.0 .AND.
     &    INDEX(ADUM(1:),'centroid').NE.0 ) THEN
          IREF(NV) = 292
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'node').NE.0 .AND.
     &    INDEX(ADUM(1:),'centroid').NE.0 ) THEN
          IREF(NV) = 293
        ELSEIF( INDEX(ADUM(1:),'gas viscosity').NE.0 ) THEN
          IREF(NV) = 289
        ELSEIF( INDEX(ADUM(1:),'aqueous enthalpy').NE.0 ) THEN
          IREF(NV) = 296
        ELSEIF( INDEX(ADUM(1:),'gas enthalpy').NE.0 ) THEN
          IREF(NV) = 297
        ELSEIF( INDEX(ADUM(1:),'aqueous internal energy').NE.0 ) THEN
          IREF(NV) = 341
        ELSEIF( INDEX(ADUM(1:),'gas internal energy').NE.0 ) THEN
          IREF(NV) = 342
        ELSEIF( INDEX(ADUM(1:),'aqueous thermal conduc').NE.0 ) THEN
          IREF(NV) = 344
        ELSEIF( INDEX(ADUM(1:),'gas thermal conduc').NE.0 ) THEN
          IREF(NV) = 345
        ELSEIF( INDEX(ADUM(1:),'salt aqueous diff').NE.0 ) THEN
          IREF(NV) = 349
        ELSEIF( INDEX(ADUM(1:),'aqueous diff').NE.0 ) THEN
          DO IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IREFGC(NV) = IGC
          END DO
          IF( IREFGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Aqueous Diffusion Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
          IREF(NV) = 347
        ELSEIF( INDEX(ADUM(1:),'h2o gas diff').NE.0 .OR.
     &    INDEX(ADUM(1:),'water gas diff').NE.0 ) THEN
          IREF(NV) = 348
#ifdef ecke
        ELSEIF( ( INDEX(ADUM(1:),'solute volumetric conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species volumetric conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'exchange' ).NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 12
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 12
            CHREF(INDX) = 'SPX'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute volumetric conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species volumetric conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'solid' ).NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 26
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 26
            CHREF(INDX) = 'SPS'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute volumetric conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species volumetric conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
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
        ELSEIF( ( INDEX(ADUM(1:),'solute aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species aqueous conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'exchange' ).NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 13
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 13
            CHREF(INDX) = 'SPLX'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species aqueous conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'solid' ).NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 27
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 27
            CHREF(INDX) = 'SPLS'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species aqueous conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
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
        ELSEIF( ( INDEX(ADUM(1:),'solute gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species gas conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'exchange' ).NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 14
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 14
            CHREF(INDX) = 'SPGS'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species gas conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'solid' ).NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 28
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 28
            CHREF(INDX) = 'SPGS'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species gas conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 3
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 3
            CHREF(INDX) = 'SPG'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute gas conc').NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 3
#endif
#ifdef ecke
        ELSEIF( ( INDEX(ADUM(1:),'solute nonaqueous-liquid').NE.0 .OR.
     &    INDEX(ADUM(1:),'species nonaqueous-liquid').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'exchange' ).NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 15
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 15
            CHREF(INDX) = 'SPNS'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute nonaqueous-liquid').NE.0 .OR.
     &    INDEX(ADUM(1:),'species nonaqueous-liquid').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'solid' ).NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 29
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 29
            CHREF(INDX) = 'SPNS'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute nonaqueous-liquid').NE.0 .OR.
     &    INDEX(ADUM(1:),'species nonaqueous-liquid').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 4
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 4
            CHREF(INDX) = 'SPN'
            UNREF(INDX) = 'mol/m^3'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute nonaqueous-liquid').NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 4
#endif
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
        ELSEIF( INDEX(ADUM(1:),'species gas conc').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 3
          IREF(NV) = INDX
          CHREF(INDX) = 'SPG'
          UNREF(INDX) = 'mol/m^3'
        ELSEIF( INDEX(ADUM(1:),'species source').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 11
          IREF(NV) = INDX
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
          UNREF(INDX) = ''
        ELSEIF( INDEX(ADUM(1:),'ph').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 27
          IREF(NV) = INDX
          CHREF(INDX) = 'pH'
#endif
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Reference Node Variable: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for duplicate reference node variables  ---
!
        DO 190 NX = 1,NV-1
          IF( IREF(NV).EQ.IREF(NX) .AND. IREFGC(NV).EQ.IREFGC(NX) ) THEN
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
!---      Gas species  ---
!
          DO 336 M = 1,NSPG
            NSP = M+NSPL+NSPS
            IF( SPNM.EQ.SPNMG(M) ) GOTO 350
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
        IPLOTGC(NV) = 0
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
        ELSEIF( INDEX(ADUM(1:),'trapped gas sat').NE.0 ) THEN
          IPLOT(NV) = 105
        ELSEIF( INDEX(ADUM(1:),'apparent aqueous sat').NE.0 ) THEN
          IPLOT(NV) = 9
        ELSEIF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
          IPLOT(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'gas saturation').NE.0 ) THEN
          IPLOT(NV) = 12
        ELSEIF( INDEX(ADUM(1:),'salt saturation').NE.0 ) THEN
          IPLOT(NV) = 264
        ELSEIF( INDEX(ADUM(1:),'aqueous moisture cont').NE.0 ) THEN
          IPLOT(NV) = 15
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mole frac').NE.0 ) THEN
          IPLOT(NV) = 205
        ELSEIF( INDEX(ADUM(1:),'aqueous mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous mole frac').NE.0 ) THEN
          IPLOT(NV) = 204
          CHREF(204) = 'XMLC'
          DO IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IPLOTGC(NV) = IGC
          ENDDO
          IF( IPLOTGC(NV).EQ.0 ) THEN
            INDX = 4
           CHMSG = 'Unrecognized Aqueous Mole Fraction Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          END IF
        ELSEIF( INDEX(ADUM(1:),'effective trapped gas').NE.0 ) THEN
          IPLOT(NV) = 19
        ELSEIF( INDEX(ADUM(1:),'diffusive porosity').NE.0 ) THEN
          IPLOT(NV) = 20
        ELSEIF( INDEX(ADUM(1:),'h2o gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o mass frac').NE.0 ) THEN
          IPLOT(NV) = 21
        ELSEIF( INDEX(ADUM(1:),'gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas mass frac').NE.0 ) THEN
          IPLOT(NV) = 22
          CHREF(22) = 'XGC'
          DO IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )
     &        IPLOTGC(NV) = IGC
          ENDDO
          IF( IPLOTGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Gas Mass Fraction Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'h2o aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous h2o mass frac').NE.0 ) THEN
          IPLOT(NV) = 24
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mass').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mass').NE.0 ) THEN
          IPLOT(NV) = 110
        ELSEIF( INDEX(ADUM(1:),'aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous mass frac').NE.0 ) THEN
          IPLOT(NV) = 25
          CHREF(25) = 'XLC'
          DO IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )
     &        IPLOTGC(NV) = IGC
          ENDDO
          IF( IPLOTGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Aqu. Mass Fraction Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
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
        ELSEIF ( INDEX(ADUM(1:),'total salt mass').NE.0 ) THEN
          IPLOT(NV) = 50
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IPLOT(NV) = 38
          DO 346 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IPLOTGC(NV) = IGC
  346     CONTINUE
          IF( IPLOTGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Total Mass Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'h2o mass source int').NE.0 ) THEN
          IPLOT(NV) = 40
        ELSEIF( INDEX(ADUM(1:),'mass source int').NE.0 ) THEN
          IPLOT(NV) = 41
          DO 348 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IPLOTGC(NV) = IGC
  348     CONTINUE
          IF( IPLOTGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Mass Source Integral Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF ( INDEX(ADUM(1:),'salt conc').NE.0 ) THEN
          IPLOT(NV) = 47
        ELSEIF ( INDEX(ADUM(1:),'salt aqueous conc').NE.0 ) THEN
          IPLOT(NV) = 48
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
        ELSEIF( INDEX(ADUM(1:),'gas mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas mole').NE.0 ) THEN
          IPLOT(NV) = 71
          CHREF(71) = 'XMGC'
          DO IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IPLOTGC(NV) = IGC
          ENDDO
          IF( IPLOTGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Gas Mole Fraction Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'h2o gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o conc').NE.0 ) THEN
          IPLOT(NV) = 73
        ELSEIF( INDEX(ADUM(1:),'gas conc').NE.0 ) THEN
          IPLOT(NV) = 74
          DO 352 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IPLOTGC(NV) = IGC
  352     CONTINUE
          IF( IPLOTGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Gas Concentration Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'h2o aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous h2o conc').NE.0 ) THEN
          IPLOT(NV) = 76
        ELSEIF( INDEX(ADUM(1:),'aqueous conc').NE.0 .AND.
     &    INDEX(ADUM(1:),'specie').EQ.0) THEN
          IPLOT(NV) = 77
          DO 354 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IPLOTGC(NV) = IGC
  354     CONTINUE
          IF( IPLOTGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Aqu. Concentration Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
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
        ELSEIF( INDEX(ADUM(1:),'osmotic press').NE.0 ) THEN
          IPLOT(NV) = 101
        ELSEIF( INDEX(ADUM(1:),'osmotic eff').NE.0 ) THEN
          IPLOT(NV) = 102
        ELSEIF( INDEX(ADUM(1:),'water').NE.0 .AND.
     &    (INDEX(ADUM(1:),'gas').NE.0 .OR. INDEX(ADUM(1:),'vapor').NE.0)
     &    .AND. INDEX(ADUM(1:),'partial').NE.0
     &    .AND. INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IPLOT(NV) = 128
        ELSEIF( (INDEX(ADUM(1:),'gas').NE.0 
     &      .OR. INDEX(ADUM(1:),'vapor').NE.0)
     &     .AND. INDEX(ADUM(1:),'partial').NE.0
     &     .AND. INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IPLOT(NV) = 129
          DO 229 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IPLOTGC(NV) = IGC
  229     CONTINUE
          IF( IPLOTGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas-aqueous scaling').NE.0 ) THEN
          IPLOT(NV) = 131
        ELSEIF( INDEX(ADUM(1:),'h2o mass source rate').NE.0 ) THEN
          IPLOT(NV) = 140
        ELSEIF( INDEX(ADUM(1:),'mass source rate').NE.0 ) THEN
          IPLOT(NV) = 141
          DO 157 IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IPLOTGC(NV) = IGC
  157     CONTINUE
          IF( IPLOTGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Mass Source Rate Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'salt mass source rate').NE.0 ) THEN
          IPLOT(NV) = 147
        ELSEIF( INDEX(ADUM(1:),'salt mass source int').NE.0 ) THEN
          IPLOT(NV) = 148
        ELSEIF( INDEX(ADUM(1:),'aqueous viscosity').NE.0 ) THEN
          IPLOT(NV) = 176
        ELSEIF( INDEX(ADUM(1:),'gas viscosity').NE.0 ) THEN
          IPLOT(NV) = 289
        ELSEIF( INDEX(ADUM(1:),'x intrinsic perm').NE.0 ) THEN
          IPLOT(NV) = 247
        ELSEIF( INDEX(ADUM(1:),'y intrinsic perm').NE.0 ) THEN
          IPLOT(NV) = 248
        ELSEIF( INDEX(ADUM(1:),'z intrinsic perm').NE.0 ) THEN
          IPLOT(NV) = 249
        ELSEIF( INDEX(ADUM(1:),'aqueous enthalpy').NE.0 ) THEN
          IPLOT(NV) = 296
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'node').NE.0 .AND.
     &    INDEX(ADUM(1:),'centroid').NE.0 ) THEN
          IPLOT(NV) = 291
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'node').NE.0 .AND.
     &    INDEX(ADUM(1:),'centroid').NE.0 ) THEN
          IPLOT(NV) = 292
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'node').NE.0 .AND.
     &    INDEX(ADUM(1:),'centroid').NE.0 ) THEN
          IPLOT(NV) = 293
        ELSEIF( INDEX(ADUM(1:),'gas enthalpy').NE.0 ) THEN
          IPLOT(NV) = 297
        ELSEIF( INDEX(ADUM(1:),'aqueous internal energy').NE.0 ) THEN
          IPLOT(NV) = 341
        ELSEIF( INDEX(ADUM(1:),'gas internal energy').NE.0 ) THEN
          IPLOT(NV) = 342
        ELSEIF( INDEX(ADUM(1:),'aqueous thermal conduc').NE.0 ) THEN
          IPLOT(NV) = 344
        ELSEIF( INDEX(ADUM(1:),'gas thermal conduc').NE.0 ) THEN
          IPLOT(NV) = 345
        ELSEIF( INDEX(ADUM(1:),'salt aqueous diff').NE.0 ) THEN
          IPLOT(NV) = 349
        ELSEIF( INDEX(ADUM(1:),'aqueous diff').NE.0 ) THEN
          DO IGC = 1,NGC
            NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
            GCNMX = GCNM(IGC)
            IF( INDEX(ADUM(1:),GCNMX(1:NCHX)).NE.0 ) IPLOTGC(NV) = IGC
          END DO
          IF( IPLOTGC(NV).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Aqueous Diffusion Component Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          END IF
          IPLOT(NV) = 347
        ELSEIF( INDEX(ADUM(1:),'h2o gas diff').NE.0 .OR.
     &    INDEX(ADUM(1:),'water gas diff').NE.0 ) THEN
          IPLOT(NV) = 348
#ifdef ecke
        ELSEIF( ( INDEX(ADUM(1:),'solute volumetric conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species volumetric conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'exchange' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 12
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 12
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute volumetric conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species volumetric conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'solid' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 26
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 26
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute volumetric conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species volumetric conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 1
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 1
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 1
#endif
#ifdef ecke
        ELSEIF( ( INDEX(ADUM(1:),'solute aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species aqueous conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'exchange' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 13
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 13
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species aqueous conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'solid' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 27
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 27
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species aqueous conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 2
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 2
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 2
#endif
#ifdef ecke
        ELSEIF( ( INDEX(ADUM(1:),'solute gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species gas conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'exchange' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 14
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 14
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species gas conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 28
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 28
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species gas conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 3
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 3
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
#else

        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 3
#endif
#ifdef ecke
        ELSEIF( ( INDEX(ADUM(1:),'solute nonaqueous-liquid').NE.0 .OR.
     &    INDEX(ADUM(1:),'species nonaqueous-liquid').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'exchange' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 15
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 15
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute nonaqueous-liquid').NE.0 .OR.
     &    INDEX(ADUM(1:),'species nonaqueous-liquid').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 29
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 29
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute nonaqueous-liquid').NE.0 .OR.
     &    INDEX(ADUM(1:),'species nonaqueous-liquid').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 4
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 4
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
#else

        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous-liquid').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 4
#endif
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
        ELSEIF( INDEX(ADUM(1:),'species gas conc').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 3
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
#endif
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Plot File Variable: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for duplicate plot file variables  ---
!
        DO 390 NX = 1,NV-1
          IF( IPLOT(NV).EQ.IPLOT(NX) .AND.
     &        IPLOTGC(NV).EQ.IPLOTGC(NX) ) THEN
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
!---  End of RDOU_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSF_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Read input file surface flux information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
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
      SUB_LOG(ISUB_LOG) = '/RDSF_COMP'
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
        ELSEIF( INDEX(ADUM(1:),'total').NE.0) THEN
          IF( INDEX(ADUM(1:),'advective').NE.0) THEN
            DO 10 IGC = 1,NGC
              NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
              IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )  THEN
                ISFGC(NS) = IGC
                GOTO 12
              ENDIF
   10       CONTINUE
            IF( ISFGC(NS).EQ.0 ) THEN
              INDX = 4
              CHMSG = 'Unrecognized Total-Advective Mass '
     &          // 'Flux Surface Component Name: ' // ADUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
   12       CONTINUE
            ISFT(NS) = 50
            WRITE(IWR,'(A)') 'Total-Advective ' //
     &        GCNM(IGC)(1:NCHX) // '-Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'diffusive').NE.0) THEN
            DO 20 IGC = 1,NGC
              NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
              IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )  THEN
                ISFGC(NS) = IGC
                GOTO 22
              ENDIF
   20       CONTINUE
            IF( ISFGC(NS).EQ.0 ) THEN
              INDX = 4
              CHMSG = 'Unrecognized Total-Diffusive Mass '
     &          // 'Flux Surface Component Name: ' // ADUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
   22       CONTINUE
            ISFT(NS) = 51
            WRITE(IWR,'(A)') 'Total-Diffusive ' //
     &        GCNM(IGC)(1:NCHX) // '-Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            DO 25 IGC = 1,NGC
              NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
              IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )  THEN
                ISFGC(NS) = IGC
                GOTO 27
              ENDIF
   25       CONTINUE
            IF( ISFGC(NS).EQ.0 ) THEN
              INDX = 4
              CHMSG = 'Unrecognized Total-Advective Mass '
     &          // 'Flux Surface Component Name: ' // ADUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
   27       CONTINUE
            ISFT(NS) = 54
            WRITE(IWR,'(A)') 'Total-Advective ' //
     &        GCNM(IGC)(1:NCHX) // '-Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0) THEN
          IF( INDEX(ADUM(1:),'volum').NE.0) THEN
            ISFT(NS) = 2
            WRITE(IWR,'(A)') 'Aqueous-Phase Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ELSEIF( INDEX(ADUM(1:),'advective').NE.0) THEN
            DO 30 IGC = 1,NGC
              NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
              IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )  THEN
                ISFGC(NS) = IGC
                GOTO 32
              ENDIF
   30       CONTINUE
            IF( ISFGC(NS).EQ.0 ) THEN
              INDX = 4
              CHMSG = 'Unrecognized Aqueous-Advective Mass '
     &          // 'Flux Surface Component Name: ' // ADUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
   32       CONTINUE
            ISFT(NS) = 46
            WRITE(IWR,'(A)') 'Aqueous-Advective ' //
     &        GCNM(IGC)(1:NCHX) // '-Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'diffusive').NE.0) THEN
            DO 40 IGC = 1,NGC
              NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
              IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )  THEN
                ISFGC(NS) = IGC
                GOTO 42
              ENDIF
   40       CONTINUE
            IF( ISFGC(NS).EQ.0 ) THEN
              INDX = 4
              CHMSG = 'Unrecognized Aqueous-Diffusive Mass '
     &          // 'Flux Surface Component Name: ' // ADUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
   42       CONTINUE
            ISFT(NS) = 47
            WRITE(IWR,'(A)') 'Aqueous-Diffusive ' //
     &        GCNM(IGC)(1:NCHX) // '-Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            ISFT(NS) = 5
            WRITE(IWR,'(A)') 'Aqueous-Phase Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
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
            ELSE
              DO 50 IGC = 1,NGC
                NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
                IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )  THEN
                  ISFGC(NS) = IGC
                  GOTO 52
                ENDIF
   50         CONTINUE
              IF( ISFGC(NS).EQ.0 ) THEN
                INDX = 4
                CHMSG = 'Unrecognized Gas-Advective Mass '
     &            // 'Flux Surface Component Name: ' // ADUM(1:NCH)
                CALL WRMSGS( INDX )
              ENDIF
   52         CONTINUE
              ISFT(NS) = 48
              WRITE(IWR,'(A)') 'Gas-Advective ' //
     &          GCNM(IGC)(1:NCHX) // '-Mass Flux Surface'
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
            ELSE
              DO 60 IGC = 1,NGC
                NCHX = INDEX(GCNM(IGC)(1:),'  ') - 1
                IF( INDEX(ADUM(1:),GCNM(IGC)(1:NCHX)).NE.0 )  THEN
                  ISFGC(NS) = IGC
                  GOTO 62
                ENDIF
   60         CONTINUE
              IF( ISFGC(NS).EQ.0 ) THEN
                INDX = 4
                CHMSG = 'Unrecognized Gas-Diffusive Mass ' //
     &            'Flux Surface Component Name: ' // ADUM(1:NCH)
                CALL WRMSGS( INDX )
              ENDIF
   62         CONTINUE
              ISFT(NS) = 49
              WRITE(IWR,'(A)') 'Gas-Diffusive ' //
     &          GCNM(IGC)(1:NCHX) // '-Mass Flux Surface'
              UNSF(1,NS) = 'kg/s'
              UNSF(2,NS) = 'kg'
            ENDIF
          ELSE
            ISFT(NS) = 6
            WRITE(IWR,'(A)') 'Gas-Phase Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'salt').NE.0) THEN
          ISFT(NS) = 8
          WRITE(IWR,'(A)') 'Salt Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          DO 110 NSL = 1,NSOLU
            IDB = INDEX(SOLUT(NSL)(1:),'  ')
            IF( INDEX(BDUM(1:),SOLUT(NSL)(1:IDB)).NE.0 ) THEN
              ISFT(NS) = NSL+100
              WRITE(IWR,'(2X,2A)') SOLUT(NSL),' Flux Surface'
              UNSF(1,NS) = 'sol/s'
              UNSF(2,NS) = 'sol'
              GOTO 120
            ENDIF
  110     CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Solute Name: '//BDUM
            CALL WRMSGS( INDX )
  120     CONTINUE
#ifdef ecke
!
!---    Conservation-component species surface flux input  ---
!
        ELSEIF( INDEX(ADUM(1:),'conservation').NE.0 .AND.
     &    INDEX(ADUM(1:),'component').NE.0 ) THEN
          VARB = 'Conservation-Component Species Name: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          DO 140NSL = NSOLU+1,NSOLU+NEQC
            IDB = INDEX(SOLUT(NSL)(1:),'  ')
            IF( BDUM(1:IDB).EQ.SOLUT(NSL)(1:IDB) ) THEN
              ISFT(NS) = NSL+100
              WRITE(IWR,'(2X,2A)') SOLUT(NSL),' Flux Surface'
              UNSF(1,NS) = 'mol/s'
              UNSF(2,NS) = 'mol'
              GOTO 160
            ENDIF
  140     CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Conservation-Component ' //
     &       'Species Name: '//BDUM
            CALL WRMSGS( INDX )
  160     CONTINUE
!
!---    Kinetic-component species surface flux input  ---
!
        ELSEIF( INDEX(ADUM(1:),'kinetic').NE.0 .AND.
     &    INDEX(ADUM(1:),'component').NE.0 ) THEN
          VARB = 'Kinetic-Component Species Name: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          DO 180 NSL = NSOLU+NEQC+1,NSOLU+NEQC+NEQK
            IDB = INDEX(SOLUT(NSL)(1:),'  ')
            IF( BDUM(1:IDB).EQ.SOLUT(NSL)(1:IDB) ) THEN
              ISFT(NS) = NSL+100
              WRITE(IWR,'(2X,2A)') SOLUT(NSL),' Flux Surface'
              UNSF(1,NS) = 'mol/s'
              UNSF(2,NS) = 'mol'
              GOTO 200
            ENDIF
  180     CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Kinetic-Component ' //
     &       'Species Name: '//BDUM
            CALL WRMSGS( INDX )
  200     CONTINUE
#endif

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
!---  End of RDSF_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSP_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Read input file for rock/soil saturation function information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TABL
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
      CHARACTER*64 ADUM,RDUM,UNTS
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSP_COMP'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Rock/Soil Saturation Function Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Loop over the rock/soil saturation information lines  ---
!
      N = 0
      IJK = 0
   10 CONTINUE
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
!---  Read saturation/capillary pressure function  ---
!
        VARB = 'Saturation Function Type: '
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        ISMX = 0
        IF( INDEX(ADUM(1:),'extended').NE.0 ) ISMX = 1
        IF( INDEX(ADUM(1:),'webb').NE.0 ) ISMX = 2
        IF( INDEX(ADUM(1:),'entrap').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            ISCHRX = 101
          ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &      INDEX(ADUM(1:),'corey').NE.0 ) THEN
            ISCHRX = 102
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Saturation Function: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'webb').NE.0 .OR.
     &    INDEX(ADUM(1:),'extension').NE.0 ) THEN
          ISMX = 2
          IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            ISCHRX = 1
          ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &      INDEX(ADUM(1:),'corey').NE.0 ) THEN
            ISCHRX = 2
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Saturation Function: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
          ISCHRX = 1
          IF( INDEX(RDUM(1:),'fractured').NE.0 .OR.
     &      INDEX(RDUM(1:),'dp').NE.0 .OR.
     &      INDEX(RDUM(1:),'dual').NE.0 ) ISCHRX = 3
        ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &      INDEX(ADUM(1:),'corey').NE.0 ) THEN
          ISCHRX = 2
          IF( INDEX(RDUM(1:),'fractured').NE.0 .OR.
     &      INDEX(RDUM(1:),'dp').NE.0 .OR.
     &      INDEX(RDUM(1:),'dual').NE.0 ) ISCHRX = 4
        ELSEIF( INDEX(ADUM(1:),'haverkamp').NE.0 ) THEN
          ISCHRX = 5
        ELSEIF( INDEX(ADUM(1:),'russo').NE.0 ) THEN
          ISCHRX = 9
        ELSEIF( INDEX(ADUM(1:),'tabular').NE.0 ) THEN
          IF( INDEX( ADUM(1:),'spline' ).NE.0 ) THEN
            IF( INDEX( ADUM(1:),'log' ).NE.0 ) THEN
              ISCHRX = 13
              WRITE(IWR,'(A)') 'Cubic-Spline-Log Interpolation'
            ELSE
              ISCHRX = 11
              WRITE(IWR,'(A)') 'Cubic-Spline Interpolation'
            ENDIF
          ELSE
            IF( INDEX( ADUM(1:),'log' ).NE.0 ) THEN
              ISCHRX = 12
              WRITE(IWR,'(A)') 'Linear-Log Interpolation'
            ELSE
              ISCHRX = 10
              WRITE(IWR,'(A)') 'Linear Interpolation'
            ENDIF
          ENDIF
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Saturation Function: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
        IF( IJK.GT.0 ) THEN
          DO 230 N = 1,NFLD
            ISCHR(IZ(N)) = ISCHRX
            IF( ISCHRX.EQ.7 ) IRPL(IZ(N)) = 7
            ISM(IZ(N)) = ISMX
  230     CONTINUE
        ELSE
          ISCHR(IROCK) = ISCHRX
          IF( ISCHRX.EQ.7 ) IRPL(IROCK) = 7
          ISM(IROCK) = ISMX
        ENDIF
!
!---  van Genuchten Function  ---
!
        IF( ISCHRX.EQ.1 .OR. ISCHRX.EQ.101 ) THEN
          IF( ISCHRX.EQ.1 ) THEN
            WRITE(IWR,'(A)') 'van Genuchten Function'
          ELSEIF( ISCHRX.EQ.101 ) THEN
            WRITE(IWR,'(A)') 'van Genuchten Function w/ Gas Entrapment'
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
          IF( IJK.GT.0 ) THEN
            DO 240 N = 1,NFLD
              SCHR(9,IZ(N)) = 1.D+20
  240       CONTINUE
          ELSE
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
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(14,IROCK)
          ENDIF
          IF( ISCHRX.EQ.101 ) THEN
            VARB = 'van Genuchten (Effective Gas Residual Saturation)'
            IF( IJK.GT.0 ) THEN
              INDX = 15
              LNDX = LSCHR
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(15,IROCK))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          SCHR(15,IROCK)
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
          ENDIF
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
            VARB = 'van Genuchten (Reference Interfacial Tension)'
            TX = 20.D+0
            XLSX = 0.D+0
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            IF( IJK.GT.0 ) THEN
              INDX = 16
              LNDX = LSCHR
              UNTS = 'n/m'
              IUNKG = 1
              IUNS = -2
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              IDFLT = 1
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(16,IROCK))
              UNTS = 'n/m'
              IDFLT = 1
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',SCHR(16,IROCK)
              INDX = 0
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,SCHR(16,IROCK),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(16,IROCK),', N/m)'
            ENDIF
          ENDIF
!
!---  Brooks and Corey Function  ---
!
        ELSEIF( ISCHRX.EQ.2 .OR. ISCHRX.EQ.102 ) THEN
          IF( ISCHRX.EQ.2 ) THEN
            WRITE(IWR,'(A)') 'Brooks and Corey Function'
          ELSEIF( ISCHRX.EQ.102 ) THEN
           WRITE(IWR,'(A)')'Brooks and Corey Function w/ Gas Entrapment'
          ENDIF
          VARB = 'Brooks and Corey (psi): '
          IF( IJK.GT.0 ) THEN
            INDX = 1
            LNDX = LSCHR
            UNTS = 'm'
            IUNM = 1
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(1,IROCK))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHR(1,IROCK)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,SCHR(1,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(1,IROCK),', m)'
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
          IF( IJK.GT.0 ) THEN
            DO 250 N = 1,NFLD
              SCHR(9,IZ(N)) = 1.D+20
  250       CONTINUE
          ELSE
            SCHR(9,IROCK) = 1.D+20
          ENDIF
          IF( ISCHRX.EQ.102 ) THEN
            VARB = 'Brooks and Corey (Effective Gas Residual ' //
     &        'Saturation)'
            IF( IJK.GT.0 ) THEN
              INDX = 15
              LNDX = LSCHR
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(15,IROCK))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          SCHR(15,IROCK)
            ENDIF
            VARB = 'Brooks and Corey (Critical Trapping Number)'
            IF( IJK.GT.0 ) THEN
              INDX = 9
              LNDX = LSCHR
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              IDFLT = 1
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(9,IROCK))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          SCHR(9,IROCK)
            ENDIF
          ENDIF
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
            VARB = 'Brooks and Corey (Reference Interfacial Tension)'
            TX = 20.D+0
            XLSX = 0.D+0
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            IF( IJK.GT.0 ) THEN
              INDX = 16
              LNDX = LSCHR
              UNTS = 'n/m'
              IUNKG = 1
              IUNS = -2
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              IDFLT = 1
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(16,IROCK))
              UNTS = 'n/m'
              IDFLT = 1
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',SCHR(16,IROCK)
              INDX = 0
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,SCHR(16,IROCK),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(16,IROCK),', N/m)'
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
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(1,IROCK))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHR(1,IROCK)
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,SCHR(1,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(1,IROCK),', 1/m)'
          ENDIF
          VARB = 'Matrix van Genuchten (n): '
          IF( IJK.GT.0 ) THEN
            INDX = 3
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(3,IROCK))
            WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHR(3,IROCK)
          ENDIF
          VARB = 'Matrix van Genuchten (residual saturation): '
          IF( IJK.GT.0 ) THEN
            INDX = 4
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(4,IROCK))
            WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHR(4,IROCK)
          ENDIF
          VARB = 'Fracture van Genuchten (alpha), '
          IF( IJK.GT.0 ) THEN
            INDX = 5
            LNDX = LSCHR
            UNTS = '1/m'
            IUNM = -1
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(5,IROCK))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHR(5,IROCK)
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,SCHR(5,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(5,IROCK),', 1/m)'
          ENDIF
          VARB = 'Fracture van Genuchten (n): '
          IF( IJK.GT.0 ) THEN
            INDX = 6
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(6,IROCK))
            WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHR(6,IROCK)
          ENDIF
          VARB = 'Fracture van Genuchten (residual saturation): '
          IF( IJK.GT.0 ) THEN
            INDX = 7
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(7,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(7,IROCK)
          ENDIF
          VARB = 'Matrix van Genuchten (m)'
          IF( IJK.GT.0 ) THEN
            INDX = 14
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(14,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(14,IROCK)
          ENDIF
          VARB = 'Fracture van Genuchten (m)'
          IF( IJK.GT.0 ) THEN
            INDX = 15
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(15,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(15,IROCK)
          ENDIF
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
            VARB = 'Dual Porosity van Genuchten ' //
     &        '(Reference Interfacial Tension)'
            TX = 20.D+0
            XLSX = 0.D+0
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            IF( IJK.GT.0 ) THEN
              INDX = 16
              LNDX = LSCHR
              UNTS = 'n/m'
              IUNKG = 1
              IUNS = -2
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              IDFLT = 1
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(16,IROCK))
              UNTS = 'n/m'
              IDFLT = 1
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',SCHR(16,IROCK)
              INDX = 0
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,SCHR(16,IROCK),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(16,IROCK),', N/m)'
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
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(1,IROCK))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHR(1,IROCK)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,SCHR(1,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(1,IROCK),', m)'
          ENDIF
          VARB = 'Matrix Brooks and Corey (lambda): '
          IF( IJK.GT.0 ) THEN
            INDX = 3
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(3,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(3,IROCK)
          ENDIF
          VARB = 'Matrix Brooks and Corey (residual saturation): '
          IF( IJK.GT.0 ) THEN
            INDX = 4
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(4,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(4,IROCK)
          ENDIF
          VARB = 'Fracture Brooks and Corey (psi)'
          IF( IJK.GT.0 ) THEN
            INDX = 5
            LNDX = LSCHR
            UNTS = 'm'
            IUNM = 1
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(5,IROCK))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHR(5,IROCK)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,SCHR(5,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(5,IROCK),', m)'
          ENDIF
          VARB = 'Fracture Brooks and Corey (lambda): '
          IF( IJK.GT.0 ) THEN
            INDX = 6
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(6,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(6,IROCK)
          ENDIF
          VARB = 'Fracture Brooks and Corey (residual saturation): '
          IF( IJK.GT.0 ) THEN
            INDX = 7
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(7,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(7,IROCK)
          ENDIF
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
            VARB = 'Dual Porosity Brooks and Corey ' //
     &        '(Reference Interfacial Tension)'
            TX = 20.D+0
            XLSX = 0.D+0
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            IF( IJK.GT.0 ) THEN
              INDX = 16
              LNDX = LSCHR
              UNTS = 'n/m'
              IUNKG = 1
              IUNS = -2
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              IDFLT = 1
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(16,IROCK))
              UNTS = 'n/m'
              IDFLT = 1
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',SCHR(16,IROCK)
              INDX = 0
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,SCHR(16,IROCK),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(16,IROCK),', N/m)'
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
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(1,IROCK))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHR(1,IROCK)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,SCHR(1,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(1,IROCK),', m)'
          ENDIF
          INDX = 0
          IUNM = 1
          VAR = 1.D+0
          CALL RDUNIT(UNTS,VAR,INDX)
          IF( IJK.GT.0 ) THEN
            DO 260 N = 1,NFLD
              SCHR(5,IZ(N)) = VAR
  260       CONTINUE
          ELSE
            SCHR(5,IROCK) = VAR
          ENDIF
          VARB = 'Haverkamp (alpha): '
          IF( IJK.GT.0 ) THEN
            INDX = 2
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(2,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(2,IROCK)
          ENDIF
          VARB = 'Haverkamp (beta): '
          IF( IJK.GT.0 ) THEN
            INDX = 3
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(3,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(3,IROCK)
          ENDIF
          VARB = 'Haverkamp (residual saturation): '
          IF( IJK.GT.0 ) THEN
            INDX = 4
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(4,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(4,IROCK)
          ENDIF
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
            VARB = 'Haverkamp (Reference Interfacial Tension)'
            TX = 20.D+0
            XLSX = 0.D+0
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            IF( IJK.GT.0 ) THEN
              INDX = 16
              LNDX = LSCHR
              UNTS = 'n/m'
              IUNKG = 1
              IUNS = -2
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              IDFLT = 1
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(16,IROCK))
              UNTS = 'n/m'
              IDFLT = 1
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',SCHR(16,IROCK)
              INDX = 0
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,SCHR(16,IROCK),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(16,IROCK),', N/m)'
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
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(1,IROCK))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHR(1,IROCK)
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,SCHR(1,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(1,IROCK),', 1/m)'
          ENDIF
          VARB = 'Russo (n): '
          IF( IJK.GT.0 ) THEN
            INDX = 3
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(3,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(3,IROCK)
          ENDIF
          VARB = 'Russo (residual saturation): '
          IF( IJK.GT.0 ) THEN
            INDX = 4
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(4,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(4,IROCK)
          ENDIF
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
            VARB = 'Russo (Reference Interfacial Tension)'
            TX = 20.D+0
            XLSX = 0.D+0
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            IF( IJK.GT.0 ) THEN
              INDX = 16
              LNDX = LSCHR
              UNTS = 'n/m'
              IUNKG = 1
              IUNS = -2
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              IDFLT = 1
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(16,IROCK))
              UNTS = 'n/m'
              IDFLT = 1
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',SCHR(16,IROCK)
              INDX = 0
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,SCHR(16,IROCK),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(16,IROCK),', N/m)'
            ENDIF
          ENDIF
!
!---    Tabular  ---
!
        ELSEIF( ISCHRX.EQ.10 .OR. ISCHRX.EQ.11 ) THEN
          IF( INDEX( ADUM(1:),'content' ).NE. 0 ) THEN
            WRITE(IWR,'(A)') 'Tabular Moisture Content versus ' //
     &        'Gas-Water Capillary Head'
          ELSE
            WRITE(IWR,'(A)') 'Tabular Aqueous Saturation versus ' //
     &        'Gas-Water Capillary Head'
          ENDIF
          VARB = 'Number of Table Entries'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
          WRITE(IWR,'(2A,I6)') VARB(1:IVR),': ',NLIN
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
            VARB = 'Tabular (Reference Interfacial Tension)'
            TX = 20.D+0
            XLSX = 0.D+0
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            IF( IJK.GT.0 ) THEN
              INDX = 16
              LNDX = LSCHR
              UNTS = 'n/m'
              IUNKG = 1
              IUNS = -2
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
            ELSE
              IDFLT = 1
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(16,IROCK))
              UNTS = 'n/m'
              IDFLT = 1
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',SCHR(16,IROCK)
              INDX = 0
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,SCHR(16,IROCK),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(16,IROCK),', N/m)'
            ENDIF
          ENDIF
          IF( NLIN.LT.2 ) THEN
            INDX = 4
            CHMSG = 'Saturation Invalid Table'
            CALL WRMSGS( INDX )
          ENDIF
          IF( IJK.GT.0 ) THEN
            VARB = 'Capillary Head'
            UNTS = 'm'
            IUNM = 1
            NTBLX = NTBL
            ILOG = 0
            CALL RDIJKT( ISTART,IJK,CHDUM,UNTS,TBLX,ISLTBL,NLIN,
     &        NTBLX,ILOG )
            IF( INDEX( ADUM(1:),'content' ).NE. 0 ) THEN
              VARB = 'Water Content'
            ELSE
              VARB = 'Saturation'
            ENDIF
            UNTS = 'null'
            ILOG = 0
            CALL RDIJKT( ISTART,IJK,CHDUM,UNTS,TBLY,ISLTBL,NLIN,
     &        NTBLX,ILOG )
            IF( INDEX( ADUM(1:),'content' ).NE. 0 ) THEN
              DO 280 N = 1,NFLD
                DO 280 M = ISLTBL(1,N),ISLTBL(2,N)
                  TBLY(M) = TBLY(M)/POR(2,N)
  280         CONTINUE
            ENDIF
          ELSE
            ISLTBL(1,IROCK) = NTBL + 1
            DO 400 NL = 1,NLIN
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
              IF( INDEX( ADUM(1:),'content' ).NE. 0 ) THEN
                VARB = 'Water Content'
              ELSE
                VARB = 'Saturation'
              ENDIF
              CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLY(NTBL))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',TBLY(NTBL)
              IF( INDEX( ADUM(1:),'content' ).NE. 0 ) THEN
                TBLY(NTBL) = TBLY(NTBL)/POR(2,IROCK)
              ENDIF
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
  400       CONTINUE
            ISLTBL(2,IROCK) = NTBL
            IF( ISCHRX.EQ.11 ) THEN
              CALL SPLINY( ISLTBL(1,IROCK),ISLTBL(2,IROCK) )
              CALL SPLINX( ISLTBL(1,IROCK),ISLTBL(2,IROCK) )
            ENDIF
          ENDIF
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
              IF( ISCHR(M).EQ.7 ) IRPL(M) = 7
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
        GOTO 10
 500  CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSP33 group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSR_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Read input file for source information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      CHARACTER*64 ADUM,BDUM,CDUM,FDUM,FMDUM,UNTS
      CHARACTER*512 CHDUM,CHDUMX
      LOGICAL FCHK
      REAL*8 VAR(LSTM,8+LSOLU+LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSR_COMP'
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
        IF( INDEX(ADUM(1:),'power').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'density').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Power Density Source'
            ISRTX = 2
          ELSEIF( INDEX(ADUM(1:),'power').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Power Source'
            ISRTX = 1
          ENDIF
! component sources
       ELSEIF( INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &    INDEX(ADUM(1:),'vol').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Volumetric Source'
          ISRTX = 3
          VARB = 'Dissolved-Salt Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(1:),'conc').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved-Salt Aqu. Concentration'
            ISRTX = ISRTX + 100
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved-Salt Aqu. Mass Fraction'
            ISRTX = ISRTX + 300
          ELSEIF( INDEX(BDUM(1:),'rel').NE.0 .AND.
     &      INDEX(BDUM(1:),'sat').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved-Salt Aqu. Rel. Saturation'
             ISRTX = ISRTX + 200
          ELSE
            WRITE(IWR,'(2X,A)') 'No Dissolved Salt'
          ENDIF
          VARB = 'Dissolved-Gas Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHC,CHDUM,CDUM)
          IF( INDEX(CDUM(1:),'aqu').NE.0 .AND.
     &      INDEX(CDUM(1:),'conc').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved-Gas Aqu. Concentration'
            ISRTX = ISRTX + 1000
          ELSEIF( INDEX(CDUM(1:),'mass').NE.0 .AND.
     &      INDEX(CDUM(1:),'frac').NE.0 ) THEN
            ISRTX = ISRTX + 3000
          ELSEIF( INDEX(CDUM(1:),'rel').NE.0 .AND.
     &      INDEX(CDUM(1:),'sat').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved-Gas Aqu. Rel. Saturation'
            ISRTX = ISRTX + 2000
          ELSE
            WRITE(IWR,'(2X,A)') 'No Dissolved Gas'
          ENDIF
       ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &   INDEX(ADUM(1:),'vol').NE.0 ) THEN
         WRITE(IWR,'(2X,A)') 'Gas Volumetric Source'
         VARB = 'Water-Vapor Source Option'
         CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
         IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
           WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
           ISRTX = 4
           WRITE(IWR,'(2X,A)') 'Dry-Gas Component Mole Fraction'
         ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &     INDEX(BDUM(1:),'frac').NE.0 ) THEN
           WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
           ISRTX = 5
           WRITE(IWR,'(2X,A)') 'Dry-Gas Component Mass Fraction'
         ELSE
           ISRTX = 6
           WRITE(IWR,'(2X,A)') 'No Water-Vapor'
           WRITE(IWR,'(2X,A)') 'Dry-Gas Component Mass Fraction'
         ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Mass Source'
          ISRTX = 7
          VARB = 'Dissolved-Salt Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(1:),'conc').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved-Salt Aqu. Concentration'
            ISRTX = ISRTX + 100
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved-Salt Aqu. Mass Fraction'
            ISRTX = ISRTX + 300
          ELSEIF( INDEX(BDUM(1:),'rel').NE.0 .AND.
     &      INDEX(BDUM(1:),'sat').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved-Salt Aqu. Rel. Saturation'
            ISRTX = ISRTX + 200
          ELSE
            WRITE(IWR,'(2X,A)') 'No Dissolved Salt'
          ENDIF
          VARB = 'Dissolved-Gas Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHC,CHDUM,CDUM)
          IF( INDEX(CDUM(1:),'aqu').NE.0 .AND.
     &      INDEX(CDUM(1:),'conc').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved-Gas Aqu. Concentration'
            ISRTX = ISRTX + 1000
          ELSEIF( INDEX(CDUM(1:),'mass').NE.0 .AND.
     &      INDEX(CDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved-Gas Aqu. Mass Fraction'
            ISRTX = ISRTX + 3000
          ELSEIF( INDEX(CDUM(1:),'rel').NE.0 .AND.
     &      INDEX(CDUM(1:),'sat').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved-Gas Aqu. Rel. Saturation'
            ISRTX = ISRTX + 2000
          ELSE
            WRITE(IWR,'(2X,A)') 'No Dissolved Gas'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Gas Mass Source'
          VARB = 'Water-Vapor Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
            ISRTX = 8
            WRITE(IWR,'(2X,A)') 'Dry-Gas Component Mole Fraction'
           ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &       INDEX(BDUM(1:),'frac').NE.0 ) THEN
             WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
             ISRTX = 9
             WRITE(IWR,'(2X,A)') 'Dry-Gas Component Mass Fraction'
           ELSE
             ISRTX = 10
             WRITE(IWR,'(2X,A)') 'No Water-Vapor'
             WRITE(IWR,'(2X,A)') 'Dry-Gas Component Mass Fraction'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'salt').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'density').NE.0 ) THEN
            ISRTX = 11
            WRITE(IWR,'(2X,A)') 'Salt Density Source'
          ELSE
            ISRTX = 12
            WRITE(IWR,'(2X,A)') 'Salt Source'
          ENDIF
!
! Pressure limited well model
!
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
!
! Mass Rate limited well model
!
        ELSEIF( INDEX(ADUM(1:),'z-dir').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Gas Mass Rate ' //
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
!---      Read source domain indices  ---
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
!---      Check for ill-defined source domains  ---
!
          IF( (ISRTX.GE.13.AND.ISRTX.LE.15) .OR. 
     &        (ISRTX.GE.23.AND.ISRTX.LE.25)  ) THEN
            IF( ((I2X-I1X).GT.1) .OR. ((J2X-J1X).GT.1) ) THEN
              INDX = 4
              CHMSG = 'Invalid Well Source Domain'
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
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
!---      Check for sources applied to inactive nodes  ---
!
          DO 50 K = K1X,K2X
          DO 50 J = J1X,J2X
          DO 50 I = I1X,I2X
            IF( IXP(ND(I,J,K)).LE.0 ) THEN
              INDX = 7
              IMSG = ND(I,J,K)
              CHMSG = 'Source Applied to an Inactive Node: '
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
!
!---    Loop over number of source times  ---
!
        DO 100 NTM = 1,ISRM(NS)
          DO 70 M = 1,8+NSOLU+NGC
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
!
!---      Power-Density Source  ---
!
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
!
!---      Aqueous Volumetric Source  ---
!
          ELSEIF( MOD(ISRTX,100).EQ.3 ) THEN
!
!---        Nonisothermal simulations  ---
!
            IF( ISLC(30).EQ.0 ) THEN
              VARB = 'Source Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            ENDIF
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
            IF( MOD(ISRTX,1000)/100.EQ.1 ) THEN
              VARB = 'Source Dissolved-Salt Aqueous Concentration'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &          ': ',VAR(NTM,5)
              INDX = 0
              IUNM = -3
              IUNKG = 1
              CALL RDUNIT(UNTS,VAR(NTM,5),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,5),', kg/m^3)'
            ELSEIF( MOD(ISRTX,1000)/100.EQ.2 ) THEN
              VARB = 'Source Dissolved-Salt Aqueous Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( MOD(ISRTX,1000)/100.EQ.3 ) THEN
              VARB = 'Source Dissolved-Salt Aqueous Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,5)
            ENDIF
            IF( (ISRTX/1000).EQ.1 ) THEN
             DO IGC = 1, NGC
              NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
              VARB = 'Source Dissolved-' // GCNM(IGC)(1:NCH) //
     &          '  Aqueous Concentration'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &          ': ',VAR(NTM,5+IGC)
              INDX = 0
              IUNM = -3
              IUNKG = 1
              CALL RDUNIT(UNTS,VAR(NTM,5+IGC),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,5+IGC),', kg/m^3)'
             ENDDO
            ELSEIF( (ISRTX/1000).EQ.2 ) THEN
             DO IGC = 1, NGC
              NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
              VARB = 'Source Dissolved-' // GCNM(IGC)(1:NCH) //
     &          '  Aqueous Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
              WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',
     &               VAR(NTM,5+IGC)
             ENDDO
            ELSEIF( (ISRTX/1000).EQ.3 ) THEN
             DO IGC = 1, NGC
              NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
              VARB = 'Source Dissolved-' // GCNM(IGC)(1:NCH) //
     &          '  Aqueous Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
              WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ', 
     &          VAR(NTM,5+IGC)
             ENDDO
            ENDIF
!
!---      Gas Volumetric Source  ---
!
          ELSEIF( ISRTX.EQ.4 .OR. ISRTX.EQ.5 .OR. ISRTX.EQ.6 ) THEN
!
!---        Nonisothermal simulations  ---
!
            IF( ISLC(30).EQ.0 ) THEN
              VARB = 'Source Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            ENDIF
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
            IF( ISRTX.EQ.4 ) THEN
              VARB = 'Source Water Vapor Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mole Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
                WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,5+IGC)
              ENDDO
            ELSEIF( ISRTX.EQ.5 ) THEN
              VARB = 'Source Water Vapor Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
                WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,5+IGC)
              ENDDO
            ELSEIF( ISRTX.EQ.6 ) THEN
              ISRTX = 5
              VAR(NTM,5) = 0.D+0
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
                WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,5+IGC)
              ENDDO
            ENDIF
!
!---      Aqueous Mass Source  ---
!
          ELSEIF( MOD(ISRTX,10).EQ.7 ) THEN
!
!---        Nonisothermal simulations  ---
!
            IF( ISLC(30).EQ.0 ) THEN
              VARB = 'Source Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            ENDIF
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
            IF( MOD(ISRTX,100)/10.EQ.1 ) THEN
              VARB = 'Source Dissolved-Salt Aqueous Concentration'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &          ': ',VAR(NTM,5)
              INDX = 0
              IUNM = -3
              IUNKG = 1
              CALL RDUNIT(UNTS,VAR(NTM,5),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,5),', kg/m^3)'
            ELSEIF( MOD(ISRTX,1000)/100.EQ.2 ) THEN
              VARB = 'Source Dissolved-Salt Aqueous Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( MOD(ISRTX,1000)/100.EQ.3 ) THEN
              VARB = 'Source Dissolved-Salt Aqueous Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,5)
            ENDIF
            IF( (ISRTX/1000).EQ.1 ) THEN
             DO IGC = 1, NGC
              VARB = 'Source Dissolved-' // GCNM(IGC)(1:NCH) //
     &           '  Aqueous Concentration'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &          ': ',VAR(NTM,5+IGC)
              INDX = 0
              IUNM = -3
              IUNKG = 1
              CALL RDUNIT(UNTS,VAR(NTM,5+IGC),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,5+IGC),', kg/m^3)'
             ENDDO
            ELSEIF( (ISRTX/1000).EQ.2 ) THEN
              DO IGC = 1, NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source Dissolved-' // GCNM(IGC)(1:NCH) //
     &           '  Aqueous Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
                WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',
     &             VAR(NTM,5+IGC)
              ENDDO
            ELSEIF( (ISRTX/1000).EQ.3 ) THEN
              DO IGC = 1, NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source Dissolved-' // GCNM(IGC)(1:NCH) //
     &            '  Aqueous Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
                WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ', 
     &             VAR(NTM,5+IGC)
              ENDDO
            ENDIF
!
!---      Gas Mass Source  ---
!
          ELSEIF( ISRTX.EQ.8 .OR. ISRTX.EQ.9 .OR. ISRTX.EQ.10 ) THEN
!
!---        Nonisothermal simulations  ---
!
            IF( ISLC(30).EQ.0 ) THEN
              VARB = 'Source Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            ENDIF
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
            IF( ISRTX.EQ.8 ) THEN
              VARB = 'Source Water Vapor Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mole Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
                WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,5+IGC)
              ENDDO
            ELSEIF( ISRTX.EQ.9 ) THEN
              VARB = 'Source Water Vapor Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
                WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,5+IGC)
              ENDDO
            ELSEIF( ISRTX.EQ.10 ) THEN
              ISRTX = 9
              VAR(NTM,5) = 0.D+0
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5+IGC))
                WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,5+IGC)
              ENDDO
            ENDIF
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
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mole Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7+IGC))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,7+IGC)
              ENDDO
            ELSEIF( ISRTX.EQ.14.OR.ISRTX.EQ.17.OR.ISRTX.EQ.20 ) THEN
              VARB = 'Source Water Vapor Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7+IGC))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,7+IGC)
              ENDDO
            ELSEIF( ISRTX.EQ.15.OR.ISRTX.EQ.18.OR.ISRTX.EQ.21 ) THEN
              ISRTX = ISRTX-1
              VAR(NTM,5) = 0.D+0
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7+IGC))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,7+IGC)
              ENDDO
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
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,7)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,7),INDX)
!
!---      Gas Mass Rate Well Source  ---
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
            VARB = 'Source Gas Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
              INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
!
!---        Well water vapor  ---
!
            IF( ISRTX.EQ.23 ) THEN
              VARB = 'Source Water Vapor Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mole Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7+IGC))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,7+IGC)
              ENDDO
            ELSEIF( ISRTX.EQ.24 ) THEN
              VARB = 'Source Water Vapor Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7+IGC))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,7+IGC)
              ENDDO
            ELSEIF( ISRTX.EQ.25 ) THEN
              ISRTX = ISRTX-1
              VAR(NTM,5) = 0.D+0
              DO IGC = 1,NGC
                NCH = INDEX( GCNM(IGC)(1:),'  ' )-1
                VARB = 'Source ' // GCNM(IGC)(1:NCH) //
     &            ' Dry-Gas Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7+IGC))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            VAR(NTM,7+IGC)
              ENDDO
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
!
!---        Solute well source: gas mass concentration  ---
!
            IF( ISRTX.EQ.125 ) THEN
              DO 80 NSL = 1,NSOLSR(NS)
                VARB = 'Solute Gas Mass Concentration'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8+NSL))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',VAR(NTM,8+NSL)
                INDX = 0
                IUNM = 1
                CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
                WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,8+NSL),', m)'
   80         CONTINUE
            ENDIF
!
!---      Salt Density Source  ---
!
          ELSEIF( ISRTX.EQ.11 ) THEN
            VARB = 'Source Salt Density Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &        VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            IUNM = -3
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
!
!---      Salt Source  ---
!
          ELSEIF( ISRTX.EQ.12 ) THEN
            VARB = 'Source Salt Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &        VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
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
          DO 120 M = 1,8+NSOLU+NGC
            SRC(M,NTM,NSR) = VAR(NTM,M)
  120     CONTINUE
  130   CONTINUE
  140 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSR_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDST_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Reads the salt transport card dispersivities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/RDST_COMP'
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
        IF( N.GE.NROCK ) GOTO 600
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
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDST_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTF_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Reads solute/fluid interaction card for diffusion and partition
!     coefficients, and internodal diffusion term averaging scheme for
!     single phase (aqueous) solute transport equation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/RDTF_COMP'
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
!---  End of RDTF_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTP_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Reads the solute/porous media interaction card for the
!     dispersivities, half-lives, and partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/RDTP_COMP'
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
        IF( N.GE.NROCK ) GOTO 600
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
!---  End of RDTP_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Compute the maximum relative residuals
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      USE FDVT
      USE FDVS
      USE FDVP
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
      INTEGER MPC(LNGC)
      REAL*8 XGCX(LNGC),XLCX(LNGC), XLCMX(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),XLAX(LSV,LFD)
      REAL*8 XMDGCX(LNGC)
      REAL*8 XMDLCX(LNGC)
      CHARACTER*64 PH_CND(3)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Saturated w/ Dissolved CO2',
     &             'Unsaturated', 'Saturated w/ Trapped Gas'/
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL_COMP'
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
        IF( IXP(N).LE.0 ) GOTO 200
        IZN = IZ(N)
        N_DB = N
        NMD = IXP(N)
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) MPT = IM(IEQT,NMD)
        MPL = IM(IEQW,NMD)
        DO IGC = 1, NGC
         MPC(IGC) = IM(IEQGC(IGC),NMD)
        ENDDO
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) MPS = IM(IEQS,NMD)
!
!---    Skip selected nodes in the residual calculation  ---
!
        IF( ISKP(IZN).EQ.1 ) GOTO 200
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
        ENDIF
!
!---    Saturated system w/o entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - dissolved-CO2 aqueous mass fraction
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPL))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
!
!---      Ignore residual for small dissolved-CO2  ---
!
          CALL SP_B( T(2,N),XLS(2,N),PSBX )
          PVBX = PSBX
          PX = PL(2,N)+PATM
!
! --- Calculate dry gas mole fraction
!
          XMVA = 0.
          DO IGC = 1, NGC
            XMVA = XMVA + XLC(IGC,2,N)/GCPP(1,IGC)
          ENDDO
          IF ( XMVA.EQ.0 ) XMVA = 1.D+0
          DO IGC = 1, NGC
            XMVGC(IGC,2,N) = XLC(IGC,2,N)/GCPP(1,IGC)/XMVA
          END DO
          DO IGC = 1, NGC
            XMDGCX(IGC) = XMVGC(IGC,2,N)
            XMDLCX(IGC) = XMVLC(IGC,2,N)
          ENDDO
          XMDGWX = XMVGW(2,N)
          XMDLWX = XMVLW(2,N)
          CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &      XGCX,XGWX,XLCX,XLS(2,N),XLWX,XMGCX,XMGWX,XMLCX,
     &      XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &      XMDLCX,XMDLWX)
          DO IGC = 1, NGC
            XMVLC(IGC,2,N) = XMDLCX(IGC) 
          ENDDO
          XMVGW(2,N) = XMDGWX
          XMVLW(2,N) = XMDLWX
          XLAX(2, N) = 0.D+0
          DO IGC = 1, NGC
            XLAX(2,N) = XLAX(2,N) + XLCX(IGC)
          ENDDO
          IF( XLA(2,N).GT.(1.D-6*XLAX(2,N)) ) THEN
            DO IGC = 1, NGC
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGC(IGC,2,N) +
     &        RHOL(2,N)*SL(2,N)*XLC(IGC,2,N))*DTI*VOL(N)
            RSDX = MIN( ABS(BLU(MPC(IGC)))/MAX(PVC(IGC,2,N),PATM),
     &        ABS(RSDL(IEQGC(IGC),N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQGC(IGC)) ) THEN
              RSD(IEQGC(IGC)) = RSDX
              NSD(IEQGC(IGC)) = N
            ENDIF
           ENDDO
          ENDIF
!
!---    Unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPL))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
!ignore zero mole fractions
         DO IGC = 1, NGC
          IF( XMGC(IGC,2,N).GT.0.D+0 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGC(IGC,2,N) +
     &        RHOL(2,N)*SL(2,N)*XLC(IGC,2,N))*DTI*VOL(N)
            RSDX = MIN( ABS(BLU(MPC(IGC)))/(ABS(PVC(IGC,2,N))+PATM),
     &        ABS(RSDL(IEQGC(IGC),N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQGC(IGC)) ) THEN
              RSD(IEQGC(IGC)) = RSDX
              NSD(IEQGC(IGC)) = N
            ENDIF
          ENDIF
         ENDDO
!
!---    Saturated system w/ entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - trapped gas saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPL))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
         DO IGC = 1, NGC
          ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGC(IGC,2,N) +
     &      RHOL(2,N)*SL(2,N)*XLC(IGC,2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPC(IGC))),
     &      ABS(RSDL(IEQGC(IGC),N)/(ACP+SMALL)) )
           IF( RSDX.GT.RSD(IEQGC(IGC)) ) THEN
            RSD(IEQGC(IGC)) = RSDX
            NSD(IEQGC(IGC)) = N
          ENDIF
         ENDDO
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
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          NT = NSD(IEQT)
          IF( NT.GT.0 ) THEN
            NPT = NPHAZ(2,NT)
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
        IF( NW.GT.0 ) THEN
          NPW = NPHAZ(2,NW)
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
        DO IGC = 1, NGC
         NA = NSD(IEQGC(IGC))
         IF( NA.GT.0 ) THEN
           NPA = NPHAZ(2,NA)
           NCHA = INDEX( PH_CND(NPA),'  ') - 1
           NCHX = INDEX(GCNM(IGC)(1:),'  ')-1
           WRITE(ISC,'(4X,A,A,1PE11.4,A,I6,2A)') GCNM(IGC)(1:NCHX),
     &       ' Equation Maximum Residual = ',RSD(IEQGC(IGC)),
     &       ': Node = ',NA,
     &       ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
           WRITE(IWR,'(4X,A,A,1PE11.4,A,I6,2A)') GCNM(IGC)(1:NCHX),
     &       ' Equation Maximum Residual = ',RSD(IEQGC(IGC)),
     &       ': Node = ',NA,
     &       ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
           WRITE(ISC,*) PVC(IGC,2,NA),XMGC(IGC,2,NA)
           WRITE(IWR,*) PVC(IGC,2,NA),XMGC(IGC,2,NA)
         ENDIF
        ENDDO
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          NS = NSD(IEQS)
          IF( NS.GT.0 ) THEN
            NPS = NPHAZ(2,NS)
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
            XLA(2,N) = XLA(1,N)
            DO IGC = 1, NGC
              XLC(IGC,2,N) = XLC(IGC,1,N)
              PVC(IGC,2,N) = PVC(IGC,1,N)
            ENDDO
            SG(2,N) = SG(1,N)
            YLS(2,N) = YLS(1,N)
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
!---  End of RSDL_COMP group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SBND_COMP( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Modify the Jacobian matrix for the solute transport equation
!     to incorporate boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 BCX(LSPBC+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SBND_COMP'
!
!---  Loop over number of specified boundary conditions  ---
!
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
!---          Gas species ---
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
     &            TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
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
!---              Gas species ---
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
        MF = 1
        IZN = IZ(N)
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
        XVLP = SL(2,N)*PORD(2,N)
        FCLP = 0.D+0
        IF( XVLP.GT.SMALL ) FCLP = YL(N,NSL)/XVLP
        XVGP = SG(2,N)*PORD(2,N)
        FCGP = 0.D+0
        IF( XVGP.GT.SMALL ) FCGP = YG(N,NSL)/XVGP
!
!---    Solute transport only, skip calculations for reactive
!       species transport  ---
!
        XVLB = SLB(2,NB)*PORDB(2,NB)
        XVGB = SGB(2,NB)*PORDB(2,NB)
#ifdef ecke
        IF( NSL.LE.NSOLU ) THEN
#endif
!
!---      Phase fraction factors at boundary  ---
!
          IF( IPCL(NSL).EQ.2 ) THEN
            XVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))
     &        *SLB(2,NB)
          ELSE
            XVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))
          ENDIF
!
!---      Constant gas-aqueous partition coefficient  ---
!
          IF( IPCGL(NSL).EQ.0 ) THEN
            PCGLX = PCGL(1,NSL)
!
!---      Temperature dependent gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.1 ) THEN
            TK = TB(2,NB)+TABS
            PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &        + PCGL(3,NSL)*LOG(TK)
     &        + PCGL(4,NSL)*TK + PCGL(5,NSL)*TK**2 )
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
          IF( IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.14 .OR.
     &      IBCT(NSL+LUK,NB).EQ.23 ) THEN
            BCX(1) = BCX(1)/(FCL+SMALL)
          ELSEIF( IBCT(NSL+LUK,NB).EQ.9 .OR.
     &      IBCT(NSL+LUK,NB).EQ.15 .OR.
     &      IBCT(NSL+LUK,NB).EQ.43 ) THEN
            BCX(1) = BCX(1)/(FCG+SMALL)
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
            YSPGX = 0.D+0
            DO 130 NSP = 1,IEQ_C(1,NEQ)
              DO 120 NSPX = 1,IBCSP(1,NB)
                IF( IBCSP(NSPX+1,NB).EQ.IEQ_C(NSP+1,NEQ) ) THEN
!
!---              Aqueous species ---
!
                  IF( IEQ_C(NSP+1,NEQ).LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+1,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                    YSPLX = YSPLX + EQ_C(NSP,NEQ)*BCX(NSPX+1)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+2,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                    YSPGX = YSPGX + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                  ENDIF                    
                  BCX(1) = BCX(1) + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
  120         CONTINUE
  130       CONTINUE
!
!---        Linked aqueous component   ---
!
            DO IGC = 1,NGC
            IF( ISPLK(14+NSPLK+IGC).EQ.NSL ) BCX(1) = 1.D+3*XLAB(2,NB)*
     &          RHOLB(2,N)*SLB(2,NB)*PORDB(2,NB)/WTMGCX
            END DO
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
     &                IBCT(NSOLU+LUK+1,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+1,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                    YSPLX = YSPLX + EQ_K(NSP,NEQ)*BCX(NSPX+1)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+2,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                  ENDIF                    
                  BCX(1) = BCX(1) + EQ_K(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
  140         CONTINUE
  150       CONTINUE
!
!---        Linked aqueous component   ---
!
            DO IGC = 1,NGC
            IF( ISPLK(14+NSPLK+IGC).EQ.NSL ) BCX(1) = 1.D+3*XLAB(2,NB)*
     &          RHOLB(2,N)*SLB(2,NB)*PORDB(2,NB)/WTMGCX
            END DO
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
          YGBX = 1.D+0-YLBX
          FCL = 0.D+0
          IF( XVLB/EPSL.GT.EPSL ) FCL = YLBX/XVLB
          FCG = 0.D+0
          IF( XVGB/EPSL.GT.EPSL ) FCG = YGBX/XVGB
!
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = YLBX
          YGB(NB,NSL) = YGBX
          CB(NB,NSL) = BCX(1)
        ENDIF
#endif
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
          NPZ = NSZ(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULB,VLB,WLB,N,MF )
            CALL SHDP( WLB,ULB,VLB,DISPL(IZN),DISPT(IZN),DPLB )
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGB,VGB,WGB,N,MF )
            CALL SHDP( WGB,UGB,VGB,DISPL(IZN),DISPT(IZN),DPGB )
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
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
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
!---        TVD Transport for the boundary surface  ---
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
              IF( FGB.GE.ZERO ) THEN
                WCGZ = BCX(1)*FCG*FGB
              ELSEIF( FGB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCGT = YG(NBT,NSL)/(SG(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBT,NSL)*FCGT)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                WCGZ = BCX(1)*FGB*THETA*FCG
     &            + C(N,NSL)*FGB*(1.D+0-THETA)*FCGP
              ELSEIF( FGB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCGZ = C(N,NSL)*FGB*FCGP
              ENDIF
              AB = DLZ*FCL + DGZ*FCG
              AP = DGZ*FCGP + DLZ*FCLP
              WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ+WCGZ)/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCGZ
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVLX = SL(2,NBT)*PORD(2,NBT)
                FCLT = YL(NBT,NSL)/(XVLX+SMALL)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(XVLX+SMALL)
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
              FGT = AFZ(NQZ)*WG(1,NQZ)
              IF( FGT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVGX = SG(2,NBT)*PORD(2,NBT)
                FCGT = YG(NBT,NSL)/(XVGX+SMALL)
                CRGT = ABS( WG(1,NQZ) )*DT/DZGP(NQZ)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AB*BCX(1)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLB = MIN( FLB,0.D+0 )
            FGB = MIN( FGB,0.D+0 )
!
!---      TVD Transport for the boundary surface  ---
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
              WCGZ = 0.D+0
              IF( FGB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCGT = YG(NBT,NSL)/(SG(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBT,NSL)*FCGT)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                WCGZ = BCX(1)*FGB*THETA*FCG
     &            + C(N,NSL)*FGB*(1.D+0-THETA)*FCGP
              ELSEIF( FGB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCGZ = C(N,NSL)*FGB*FCGP
              ENDIF
              AB = 0.D+0
              AP = 0.D+0
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCLZ/AFZ(NPZ) + WCGZ/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCGZ
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVLX = SL(2,NBT)*PORD(2,NBT)
                FCLT = YL(NBT,NSL)/(XVLX+SMALL)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(XVLX+SMALL)
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
              FGT = AFZ(NQZ)*WG(1,NQZ)
              IF( FGT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVGX = SG(2,NBT)*PORD(2,NBT)
                FCGT = YG(NBT,NSL)/(XVGX+SMALL)
                CRGT = ABS( WG(1,NQZ) )*DT/DZGP(NQZ)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AB*BCX(1)
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLB = MAX( FLB,0.D+0 )
            FGB = MAX( FGB,0.D+0 )
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              WCLZ = 0.D+0
              IF( FLB.GE.ZERO ) WCLZ = BCX(1)*FCL*FLB
              WCGZ = 0.D+0
              IF( FGB.GE.ZERO ) WCGZ = BCX(1)*FCG*FGB
              AB = 0.D+0
              AP = 0.D+0
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCLZ/AFZ(NPZ) + WCGZ/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCGZ
!
!---          TVD Transport for interior surface 
!             adjacent to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVLX = SL(2,NBT)*PORD(2,NBT)
                FCLT = YL(NBT,NSL)/(XVLX+SMALL)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(XVLX+SMALL)
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
              FGT = AFZ(NQZ)*WG(1,NQZ)
              IF( FGT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVGX = SG(2,NBT)*PORD(2,NBT)
                FCGT = YG(NBT,NSL)/(XVGX+SMALL)
                CRGT = ABS( WG(1,NQZ) )*DT/DZGP(NQZ)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
     &        UL,VL,WL,ULS,VLS,WLS,N,MF )
            CALL SHDP( VLS,WLS,ULS,DISPL(IZN),DISPT(IZN),DPLS )
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGS,VGS,WGS,N,MF )
            CALL SHDP( VGS,WGS,UGS,DISPL(IZN),DISPT(IZN),DPGS )
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
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLY = DIFMN(DLB,DLP,DYGF(N),DYGF(N),VL(1,NPY),INDX)
            DLY = AFY(NPY)*(DLY+DPLB)/RP(I)/(5.D-1*DYGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGY = DIFMN(DGB,DGP,DYGF(N),DYGF(N),VG(1,NPY),INDX)
            DGY = AFY(NPY)*(DGY+DPGB)/RP(I)/(5.D-1*DYGF(N))
!
!---        TVD Transport for the boundary surface  ---
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
              IF( FGS.GE.ZERO ) THEN
                VCGY = BCX(1)*FCG*FGS
              ELSEIF( FGS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCGN = YG(NBN,NSL)/(SG(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBN,NSL)*FCGN)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                VCGY = BCX(1)*FGS*THETA*FCG
     &            + C(N,NSL)*FGS*(1.D+0-THETA)*FCGP
              ELSEIF( FGS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCGY = C(N,NSL)*FGS*FCGP
              ENDIF
              AS = DLY*FCL + DGY*FCG
              AP = DGY*FCGP + DLY*FCLP
              VC(NPY,NSL) = VC(NPY,NSL) + (VCLY+VCGY)/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCGY
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NQY = NPY+IFLD
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVLX = SL(2,NBN)*PORD(2,NBN)
                FCLN = YL(NBN,NSL)/(XVLX+SMALL)
                CRLN = ABS( VL(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVLX+SMALL)
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
              FGN = AFY(NQY)*VG(1,NQY)
              IF( FGN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVGX = SG(2,NBN)*PORD(2,NBN)
                FCGN = YG(NBN,NSL)/(XVGX+SMALL)
                CRGN = ABS( VG(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AS*BCX(1)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLS = MIN( FLS,0.D+0 )
            FGS = MIN( FGS,0.D+0 )
!
!---        TVD Transport for the boundary surface  ---
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
              VCGY = 0.D+0
              IF( FGS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCGN = YG(NBN,NSL)/(SG(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBN,NSL)*FCGN)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                VCGY = BCX(1)*FGS*THETA*FCG
     &            + C(N,NSL)*FGS*(1.D+0-THETA)*FCGP
              ELSEIF( FGS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCGY = C(N,NSL)*FGS*FCGP
              ENDIF
              AS = 0.D+0
              AP = 0.D+0
              VC(NPY,NSL) = VC(NPY,NSL) + VCLY/AFY(NPY) + VCGY/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCGY
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NQY = NPY+IFLD
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVLX = SL(2,NBN)*PORD(2,NBN)
                FCLN = YL(NBN,NSL)/(XVLX+SMALL)
                CRLN = ABS( VL(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVLX+SMALL)
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
              FGN = AFY(NQY)*VG(1,NQY)
              IF( FGN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVGX = SG(2,NBN)*PORD(2,NBN)
                FCGN = YG(NBN,NSL)/(XVGX+SMALL)
                CRGN = ABS( VG(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AS*BCX(1)
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLS = MAX( FLS,0.D+0 )
            FGS = MAX( FGS,0.D+0 )
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              VCLY = 0.D+0
              IF( FLS.GE.ZERO ) VCLY = BCX(1)*FCL*FLS
              VCGY = 0.D+0
              IF( FGS.GE.ZERO ) VCGY = BCX(1)*FCG*FGS
              AS = 0.D+0
              AP = 0.D+0
              VC(NPY,NSL) = VC(NPY,NSL) + VCLY/AFY(NPY) + VCGY/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCGY
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NQY = NPY+IFLD
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVLX = SL(2,NBN)*PORD(2,NBN)
                FCLN = YL(NBN,NSL)/(XVLX+SMALL)
                CRLN = ABS( VL(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVLX+SMALL)
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
              FGN = AFY(NQY)*VG(1,NQY)
              IF( FGN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVGX = SG(2,NBN)*PORD(2,NBN)
                FCGN = YG(NBN,NSL)/(XVGX+SMALL)
                CRGN = ABS( VG(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
     &        UL,VL,WL,ULX,VLX,WLX,N,MF )
            CALL SHDP( ULX,VLX,WLX,DISPL(IZN),DISPT(IZN),DPLW )
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGX,VGX,WGX,N,MF )
            CALL SHDP( UGX,VGX,WGX,DISPL(IZN),DISPT(IZN),DPGW )
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
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
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
!---        TVD Transport for the boundary surface  ---
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
              IF( FGW.GE.ZERO ) THEN
                UCGX = BCX(1)*FCG*FGW
              ELSEIF( FGW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCGE = YG(NBE,NSL)/(SG(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBE,NSL)*FCGE)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA)*FCGP
     &            + BCX(1)*FGW*THETA*FCG
              ELSEIF( FGW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCGX = C(N,NSL)*FGW*FCGP
              ENDIF
              AW = DLX*FCL + DGX*FCG
              AP = DLX*FCLP + DGX*FCGP
              UC(NPX,NSL) = UC(NPX,NSL) + (UCLX+UCGX)/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCGX
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NQX = NPX+1
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                XVLX = SL(2,NBE)*PORD(2,NBE)
                FCLE = YL(NBE,NSL)/(XVLX+SMALL)
                CRLE = ABS( UL(1,NQX) )*DT/DXGP(NQX)/(XVLX+SMALL)
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
              FGE = AFX(NQX)*UG(1,NQX)
              IF( FGE.GE.ZERO ) THEN
                NBE = N+1
                XVGX = SG(2,NBE)*PORD(2,NBE)
                FCGE = YG(NBE,NSL)/(XVGX+SMALL)
                CRGE = ABS( UG(1,NQX) )*DT/DXGP(NQX)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AW*BCX(1)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLW = MIN( FLW,0.D+0 )
            FGW = MIN( FGW,0.D+0 )
!
!---        TVD Transport for the boundary surface  ---
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
              UCGX = 0.D+0
              IF( FGW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCGE = YG(NBE,NSL)/(SG(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBE,NSL)*FCGE)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA)*FCGP
     &            + BCX(1)*FGW*THETA*FCG
              ELSEIF( FGW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCGX = C(N,NSL)*FGW*FCGP
              ENDIF
              AW = 0.D+0
              AP = 0.D+0
              UC(NPX,NSL) = UC(NPX,NSL) + UCLX/AFX(NPX) + UCGX/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCGX
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NQX = NPX+1
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                XVLX = SL(2,NBE)*PORD(2,NBE)
                FCLE = YL(NBE,NSL)/(XVLX+SMALL)
                CRLE = ABS( UL(1,NQX) )*DT/DXGP(NQX)/(XVLX+SMALL)
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
              FGE = AFX(NQX)*UG(1,NQX)
              IF( FGE.GE.ZERO ) THEN
                NBE = N+1
                XVGX = SG(2,NBE)*PORD(2,NBE)
                FCGE = YG(NBE,NSL)/(XVGX+SMALL)
                CRGE = ABS( UG(1,NQX) )*DT/DXGP(NQX)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AW*BCX(1)
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLW = MAX( FLW,0.D+0 )
            FGW = MAX( FGW,0.D+0 )
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              UCLX = 0.D+0
              IF( FLW.GE.ZERO ) UCLX = BCX(1)*FCL*FLW
              UCGX = 0.D+0
              IF( FGW.GE.ZERO ) UCGX = BCX(1)*FCG*FGW
              AW = 0.D+0
              AP = 0.D+0
              UC(NPX,NSL) = UC(NPX,NSL) + UCLX/AFX(NPX) + UCGX/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCGX
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NQX = NPX+1
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                XVLX = SL(2,NBE)*PORD(2,NBE)
                FCLE = YL(NBE,NSL)/(XVLX+SMALL)
                CRLE = ABS( UL(1,NQX) )*DT/DXGP(NQX)/(XVLX+SMALL)
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
              FGE = AFX(NQX)*UG(1,NQX)
              IF( FGE.GE.ZERO ) THEN
                NBE = N+1
                XVGX = SG(2,NBE)*PORD(2,NBE)
                FCGE = YG(NBE,NSL)/(XVGX+SMALL)
                CRGE = ABS( UG(1,NQX) )*DT/DXGP(NQX)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
     &        UL,VL,WL,ULE,VLE,WLE,N,MF )
            CALL SHDP( ULE,VLE,WLE,DISPL(IZN),DISPT(IZN),DPLE )
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGE,VGE,WGE,N,MF )
            CALL SHDP( UGE,VGE,WGE,DISPL(IZN),DISPT(IZN),DPGE )
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
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
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
!---        TVD Transport for the boundary surface  ---
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
              IF( FGE.LT.ZERO ) THEN
                UCGX = BCX(1)*FCG*FGE
              ELSEIF( FGE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCGW = YG(NBW,NSL)/(SG(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBW,NSL)*FCGW)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                UCGX =  C(N,NSL)*FGE*(1.D+0-THETA)*FCGP
     &          +  BCX(1)*FGE*THETA*FCG
              ELSEIF( FGE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCGX =  C(N,NSL)*FGE*FCGP
              ENDIF
              AE = DLX*FCL + DGX*FCG
              AP = DLX*FCLP + DGX*FCGP
              UC(NQX,NSL) = UC(NQX,NSL) + (UCLX+UCGX)/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCGX
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NPX = NSX(N)
              FLW = AFX(NPX)*UL(1,NPX)
              IF( FLW.LT.ZERO ) THEN
                NBW = N-1
                XVLX = SL(2,NBW)*PORD(2,NBW)
                CRLW = ABS( UL(1,NPX) )*DT/DXGP(NPX)/(XVLX+SMALL)
                FCLW = YL(NBW,NSL)/(XVLX+SMALL)
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
              IF( FGW.LT.ZERO ) THEN
                NBW = N-1
                XVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/DXGP(NPX)/(XVGX+SMALL)
                FCGW = YG(NBW,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AE*BCX(1)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLE = MAX( FLE,0.D+0 )
            FGE = MAX( FGE,0.D+0 )
!
!---        TVD Transport for the boundary surface  ---
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
              UCGX = 0.D+0
              IF( FGE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCGW = YG(NBW,NSL)/(SG(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBW,NSL)*FCGW)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                UCGX =  C(N,NSL)*FGE*(1.D+0-THETA)*FCGP
     &          +  BCX(1)*FGE*THETA*FCG
              ELSEIF( FGE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCGX =  C(N,NSL)*FGE*FCGP
              ENDIF
              AE = 0.D+0
              AP = 0.D+0
              UC(NQX,NSL) = UC(NQX,NSL) + UCLX/AFX(NQX) + UCGX/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCGX
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NPX = NSX(N)
              FLW = AFX(NPX)*UL(1,NPX)
              IF( FLW.LT.ZERO ) THEN
                NBW = N-1
                XVLX = SL(2,NBW)*PORD(2,NBW)
                CRLW = ABS( UL(1,NPX) )*DT/DXGP(NPX)/(XVLX+SMALL)
                FCLW = YL(NBW,NSL)/(XVLX+SMALL)
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
              IF( FGW.LT.ZERO ) THEN
                NBW = N-1
                XVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/DXGP(NPX)/(XVGX+SMALL)
                FCGW = YG(NBW,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AE*BCX(1)
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLE = MIN( FLE,0.D+0 )
            FGE = MIN( FGE,0.D+0 )
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              UCLX = 0.D+0
              IF( FLE.LT.ZERO ) UCLX = BCX(1)*FCL*FLE
              UCGX = 0.D+0
              IF( FGE.LT.ZERO ) UCGX = BCX(1)*FCG*FGE
              AE = 0.D+0
              AP = 0.D+0
              UC(NQX,NSL) = UC(NQX,NSL) + UCLX/AFX(NQX) + UCGX/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCGX
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NPX = NSX(N)
              IF( FGW.LT.ZERO ) THEN
                NBW = N-1
                XVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/DXGP(NPX)/(XVGX+SMALL)
                FCGW = YG(NBW,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
              IF( FGW.LT.ZERO ) THEN
                NBW = N-1
                XVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/DXGP(NPX)/(XVGX+SMALL)
                FCGW = YG(NBW,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
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
!---        TVD Transport for the boundary surface  ---
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
              ELSEIF( FLN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCLY =  C(N,NSL)*FLN*FCLP
              ENDIF
              IF( FGN.LT.ZERO ) THEN
                VCGY = BCX(1)*FCG*FGN
              ELSEIF( FGN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCGS = YG(NBS,NSL)/(SG(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBS,NSL)*FCGS)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                VCGY =  BCX(1)*FGN*THETA*FCG
     &            + C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ELSEIF( FGN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCGY =  C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ENDIF
              AN = DLY*FCL + DGY*FCG
              AP = DLY*FCLP + DGY*FCGP
              VC(NQY,NSL) = VC(NQY,NSL) + (VCLY+VCGY)/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCGY
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NPY = NSY(N)
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/DYGP(NPY)/(XVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(XVLX+SMALL)
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
              FGS = AFY(NPY)*VG(1,NPY)
              IF( FGS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/DYGP(NPY)/(XVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AN*BCX(1)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLN = MAX( FLN,0.D+0 )
            FGN = MAX( FGN,0.D+0 )
!
!---        TVD Transport for the boundary surface  ---
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
              ELSEIF( FLN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCLY =  C(N,NSL)*FLN*FCLP
              ENDIF
              VCGY = 0.D+0
              IF( FGN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCGS = YG(NBS,NSL)/(SG(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBS,NSL)*FCGS)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                VCGY =  BCX(1)*FGN*THETA*FCG
     &            + C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ELSEIF( FGN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCGY =  C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ENDIF
              AN = 0.D+0
              AP = 0.D+0
              VC(NQY,NSL) = VC(NQY,NSL) + VCLY/AFY(NQY) + VCGY/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCGY
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NPY = NSY(N)
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/DYGP(NPY)/(XVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(XVLX+SMALL)
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
              FGS = AFY(NPY)*VG(1,NPY)
              IF( FGS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/DYGP(NPY)/(XVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AN*BCX(1)
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLN = MIN( FLN,0.D+0 )
            FGN = MIN( FGN,0.D+0 )
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              VCLY = 0.D+0
              IF( FLN.LT.ZERO ) VCLY = BCX(1)*FCL*FLN
              VCGY = 0.D+0
              IF( FGN.LT.ZERO ) VCGY = BCX(1)*FCG*FGN
              AN = 0.D+0
              AP = 0.D+0
              VC(NQY,NSL) = VC(NQY,NSL) + VCLY/AFY(NQY) + VCGY/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCGY
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NPY = NSY(N)
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/DYGP(NPY)/(XVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(XVLX+SMALL)
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
              FGS = AFY(NPY)*VG(1,NPY)
              IF( FGS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/DYGP(NPY)/(XVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
     &        UL,VL,WL,ULT,VLT,WLT,N,MF )
            CALL SHDP( WLT,ULT,VLT,DISPL(IZN),DISPT(IZN),DPLT )
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGT,VGT,WGT,N,MF )
            CALL SHDP( WGT,UGT,VGT,DISPL(IZN),DISPT(IZN),DPGT )
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
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
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
!---        TVD Transport for the boundary surface  ---
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
              IF( FGT.LT.ZERO ) THEN
                WCGZ = BCX(1)*FCG*FGT
              ELSEIF( FGT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCGB = YG(NBB,NSL)/(SG(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBB,NSL)*FCGB)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                WCGZ =  C(N,NSL)*FGT*(1.D+0-THETA)*FCGP
     &            + BCX(1)*FGT*THETA*FCG
              ELSEIF( FGT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCGZ =  C(N,NSL)*FGT*FCGP
              ENDIF
              AT = DLZ*FCL + DGZ*FCG
              AP = DLZ*FCLP + DGZ*FCGP
              WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ+WCGZ)/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCGZ
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NPZ = NSZ(N)
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/DZGP(NPZ)/(XVLX+SMALL)
                FCLB = YL(NBB,NSL)/(XVLX+SMALL)
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
              FGB = AFZ(NPZ)*WG(1,NPZ)
              IF( FGB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/DZGP(NPZ)/(XVGX+SMALL)
                FCGB = YG(NBB,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AT*BCX(1)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLT = MAX( FLT,0.D+0 )
            FGT = MAX( FGT,0.D+0 )
!
!---        TVD Transport for the boundary surface  ---
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
              WCGZ = 0.D+0
              IF( FGT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCGB = YG(NBB,NSL)/(SG(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBB,NSL)*FCGB)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                WCGZ =  C(N,NSL)*FGT*(1.D+0-THETA)*FCGP
     &            + BCX(1)*FGT*THETA*FCG
              ELSEIF( FGT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCGZ =  C(N,NSL)*FGT*FCGP
              ENDIF
              AT = 0.D+0
              AP = 0.D+0
              WC(NQZ,NSL) = WC(NQZ,NSL) + WCLZ/AFZ(NQZ) + WCGZ/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCGZ
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NPZ = NSZ(N)
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/DZGP(NPZ)/(XVLX+SMALL)
                FCLB = YL(NBB,NSL)/(XVLX+SMALL)
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
              FGB = AFZ(NPZ)*WG(1,NPZ)
              IF( FGB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/DZGP(NPZ)/(XVGX+SMALL)
                FCGB = YG(NBB,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
            BLU(MP) = BLU(MP) + AT*BCX(1)
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLT = MIN( FLT,0.D+0 )
            FGT = MIN( FGT,0.D+0 )
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              WCLZ = 0.D+0
              IF( FLT.LT.ZERO ) WCLZ = BCX(1)*FCL*FLT
              WCGZ = 0.D+0
              IF( FGT.LT.ZERO ) WCGZ = BCX(1)*FCG*FGT
              AT = 0.D+0
              AP = 0.D+0
              WC(NQZ,NSL) = WC(NQZ,NSL) + WCLZ/AFZ(NQZ) + WCGZ/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCGZ
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              NPZ = NSZ(N)
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/DZGP(NPZ)/(XVLX+SMALL)
                FCLB = YL(NBB,NSL)/(XVLX+SMALL)
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
              FGB = AFZ(NPZ)*WG(1,NPZ)
              IF( FGB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/DZGP(NPZ)/(XVGX+SMALL)
                FCGB = YG(NBB,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
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
!---  End of SBND_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORC_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Compute source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE REACT
      USE PORMED
      USE JACOB
      USE GRID
      USE FDVT
      USE FDVS
      USE FDVP
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
      CHARACTER*132 CHMSGX(2)
      REAL*8 QMGS_X(2)
      REAL*8 XMDLCX(LNGC)
      REAL*8 SRX(8+LNGC+LSOLU),GX(2),XLCZ(LNGC)
      REAL*8 RHOGCX(LNGC),RHOGC(LNGC,LSV,LFD)
      REAL*8 RHOLCX(LNGC),PHILCX(LNGC)
      REAL*8 XGCX(LNGC),XLCX(LNGC),XMGCX(LNGC),XMLCX(LNGC)
      REAL*8 XMDGCX(LNGC),XMGCIX(LNGC)
      REAL*8 VISCGX(LNGC),HGCX(LNGC),UGCX(LNGC),UEGCX(LNGC)
      REAL*8 DHCX(LNGC)
      REAL*8 HCFG(LFZ)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORC_COMP'

      DO IGC = 1,NGC
      
#ifdef ecke
!
!---  Load component sources for linked aqueous species   ---
!
      IF( ISPLK(14+NSPLK+IGC).NE.0 ) THEN
!
!---    Load NGC sources associated with chemical
!       reactions  ---
!
        DO 20 N = 1,NFLD
          SRCGCX = SRCGC(IGC,1,N)*DTI
         DO 10 M = 2,ISVC+2
           SRCGC(IGC,M,N) = SRCGCX
   10    CONTINUE
   20   CONTINUE
!
!---    Zero source terms  ---
!
        IF( NSR.GT.0 ) THEN
          DO 40 N = 1,NFLD
            DO 30 M = 2,ISVC+2
              SRCT(M,N) = 0.D+0
              SRCW(M,N) = 0.D+0
              SRCS(M,N) = 0.D+0
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
              SRCT(M,N) = 0.D+0
              SRCW(M,N) = 0.D+0
              SRCS(M,N) = 0.D+0
              SRCGC(IGC,M,N) = 0.D+0
   50       CONTINUE
   60     CONTINUE
        ENDIF
#ifdef ecke
      ENDIF
#endif
      ENDDO
!
!---  Loop over sources  ---
!
      DO 900 NS = 1,NSR
        IF( TM.LE.SRC(1,1,NS) ) GOTO 900
        SRX(1) = TM
        IF( ISRM(NS).EQ.1 ) THEN
          DO 70 N = 1,8+NGC+NSOLU
            SRX(N) = SRC(N,1,NS)
   70     CONTINUE
        ELSE
          DO 100 M = 2,ISRM(NS)
            IF( TM.LE.SRC(1,M,NS) ) THEN
             DTSR = MIN( SRC(1,M,NS)-TM,DT )
             TFSR = (TM-0.5D+0*DTSR-SRC(1,M-1,NS))/
     &         (SRC(1,M,NS)-SRC(1,M-1,NS))
             DO 80 N = 1,8+NGC+NSOLU
               SRX(N) = SRC(N,M-1,NS) + TFSR*(SRC(N,M,NS)-SRC(N,M-1,NS))
   80        CONTINUE
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 900
        ENDIF
  110   CONTINUE

!
!---  Skip over domain loop for well sources  ---
!
      IF( ISRT(NS).GE.23 .AND. ISRT(NS).LE.25 ) GOTO 510
!
!---    Loop source domain  ---
!
        DO 500 I = ISRDM(1,NS),ISRDM(2,NS)
        DO 500 J = ISRDM(3,NS),ISRDM(4,NS)
        DO 500 K = ISRDM(5,NS),ISRDM(6,NS)
          N = ND(I,J,K)
          N_DB = N
          IZN = IZ(N)
          IF( IXP(N).LE.0 ) GOTO 500
          DO 400 M = 2,ISVC+2
            PGX = PG(M,N) + PATM
            PLX = PL(M,N) + PATM
            TX = T(M,N)
!
!---        Specified power  ---
!
            IF( ISRT(NS).EQ.1 ) THEN
              SRCT(M,N) = SRCT(M,N) + SRX(4)
!
!---        Power density  ---
!
            ELSEIF( ISRT(NS).EQ.2 ) THEN
              SRCT(M,N) = SRCT(M,N) + SRX(4)*VOL(N)
!
!---        Aqueous volumetric rate w/ dissolved-air
!           and dissolved-gas; limit dissolved components
!           to their solubility limit  ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.3 ) THEN
              IF( SRX(4).LT.0.D+0 ) THEN
                DO IGC = 1, NGC
                  SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + SRX(4)*RHOL(M,N)*
     &                  XLC(IGC,M,N)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOL(M,N)*XLW(M,N)
                SRCS(M,N) = SRCS(M,N) + SRX(4)*RHOL(M,N)*XLS(M,N)
                SRCT(M,N) = SRCT(M,N) + SRX(4)*RHOL(M,N)*HL(M,N)
              ELSE
                PORDX = 1.D+0
                SLX = 1.D+0
                TX = SRX(2)
                PX = MAX( PLX,PGX,SRX(3) )
                PCX = 0.D+0
!
!---            Dissolved-salt aqueous concentration  ---
!
                IF( MOD(ISRT(NS),1000)/100.EQ.1 ) THEN
!
!---              Dissolved-gas aqueous concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    RHOLSX = SRX(5)
                    RHOLAX = 0.D+0
                    DO IGC = 1, NGC
                      RHOLCX(IGC) = SRX(5+IGC)
                      RHOLAX = RHOLAX + RHOLCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Conc. and ' //
     &                'Dissolved-Gas Aqu. Conc.  @Node: '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_11( TX,PX,PCX,PGX,RHOLSX,RHOLAX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
!
!---              Dissolved-Gas aqueous relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    RHOLSX = SRX(5)
                    PHILAX = 0.D+0
                    DO IGC = 1, NGC
                      PHILCX(IGC) = SRX(5+NGC)
                      PHILAX = PHILCX(IGC) + PHILAX
                    ENDDO
                    CHMSGX(1) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    CHMSGX(2) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Conc. ' //
     &                'and Dissolved-Gas Aqu. Rel. Sat.  @Node: '
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_12( TX,PX,PCX,PGX,RHOLSX,PHILAX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
!
!---              Dissolved-Gas aqueous mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    RHOLSX = SRX(5)
                    XLAX = 0.D+0
                    DO IGC = 1, NGC
                     XLCX(IGC) = SRX(5+IGC)
                     XLAX = XLAX + XLCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Conc. and ' //
     &                'Dissolved-Gas Aqu. Mass Frac.  @Node: '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_13( TX,PX,PCX,PGX,RHOLSX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
!
!---              No dissolved-Gas  ---
!
                  ELSE
                    RHOLSX = SRX(5)
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Conc. and ' //
     &                'Dissolved-Gas Aqu. Mass Frac.  @Node: '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_13( TX,PX,PCX,PGX,RHOLSX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
                  ENDIF
!
!---            Dissolved-salt aqueous relative saturation  ---
!
                ELSEIF( MOD(ISRT(NS),1000)/100.EQ.2 ) THEN
!
!---              Dissolved-Gas aqueous concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    PHILSX = SRX(5)
                    RHOLAX = 0.D+0
                    DO IGC = 1,NGC
                      RHOLCX(IGC) = SRX(5+IGC)
                      RHOLAX = RHOLAX + RHOLCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Rel. Sat. and ' //
     &                'Dissolved-Gas Aqu. Conc.  @Node: '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_21( TX,PX,PCX,PGX,PHILSX,RHOLAX,
!    &                YLSX,XLSX,XLCX,CHMSGX,XMDGCX )
!
!---              Dissolved-Gas aqueous relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    PHILSX = SRX(5)
                    PHILAX = 0.D+0
                    DO IGC = 1, NGC
                       PHILCX(IGC) = SRX(5+IGC)
                       PHILAX = PHILAX + PHILCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_22( TX,PX,PCX,PGX,PHILSX,PHILAX,
!    &                YLSX,XLSX,XLCX,CHMSGX,XMDGCX )
!
!---              Dissolved-Gas aqueous mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    PHILSX = SRX(5)
                    XLAX = 0.D+0
                    DO IGC = 1, NGC
                       XLCX(IGC)  = SRX(5+IGC)
                       XLAX = XLAX + XLCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_23( TX,PX,PCX,PGX,PHILSX,
!    &                YLSX,XLSX,XLCX,CHMSGX,XMDGCX )
!
!---              No dissolved-Gas  ---
!
                  ELSE
                    PHILSX = SRX(5)
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_23( TX,PX,PCX,PGX,PHILSX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
                  ENDIF
!
!---            Dissolved-salt aqueous mass fraction  ---
!
                ELSEIF( MOD(ISRT(NS),1000)/100.EQ.3 ) THEN
!
!---              Dissolved-Gas aqueous concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    YLSX = SRX(5)
                    RHOLAX = 0.D+0
                    DO IGC = 1, NGC
                      RHOLCX(IGC) = SRX(5+IGC)
                      RHOLAX = RHOLAX + RHOLCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Mass Frac. and ' //
     &                'Dissolved-Gas Aqu. Conc.  @Node: '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_31( TX,PX,PCX,PGX,RHOLAX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
!
!---              Dissolved-Gas aqueous relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    YLSX = SRX(5)
                    PHILAX = 0.D+0
                    DO IGC = 1, NGC
                      PHILCX(IGC) = SRX(5+IGC)
                      PHILAX = PHILAX + PHILCX(IGC)
                    ENDDO

                    CHMSGX(1) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_32( TX,PX,PCX,PGX,PHILAX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
!
!---              Dissolved-Gas aqueous mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    YLSX = SRX(5)
                    XLAX = 0.D+0
                    DO IGC = 1, NGC
                     XLCX(IGC) = SRX(5+IGC)
                     XLAX = XLAX + XLCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_33( TX,PX,PCX,PGX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
!
!---              No dissolved-Gas  ---
!
                  ELSE
                    YLSX = SRX(5)
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_33( TX,PX,PCX,PGX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
                  ENDIF
!
!---            No dissolved-salt  ---
!
                ELSE
!
!---              Dissolved-Gas aqueous concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    YLSX = 0.D+0
                    RHOLAX = 0.D+0
                    DO IGC = 1, NGC
                      RHOLCX(IGC) = SRX(5+IGC)
                      RHOLAX = RHOLAX + RHOLCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Mass Frac. and ' //
     &                'Dissolved-Gas Aqu. Conc.  @Node: '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_31( TX,PX,PCX,PGX,RHOLAX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
!
!---              Dissolved-Gas aqueous relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    YLSX = 0.D+0
                    PHILAX = 0.D+0
                    DO IGC = 1, NGC
                     PHILCX(IGC) = SRX(5+IGC)
                     PHILAX = PHILAX + PHILCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_32( TX,PX,PCX,PGX,PHILAX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
!
!---              Dissolved-Gas aqueous mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    YLSX = 0.D+0
                    XLAX = 0.D+0
                    DO IGC = 1, NGC
                     XLCX(IGC) = SRX(5+IGC)
                     XLAX = XLAX + XLCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_33( TX,PX,PCX,PGX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
!
!---              No dissolved-Gas  ---
!
                  ELSE
                    YLSX = 0.D+0
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Aqu. Vol. Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    DO IGC = 1, NGC
                     XMDGCX(IGC) = XMVGC(IGC,2,N)
                    ENDDO
!                   CALL FLSH_33( TX,PX,PCX,PGX,
!    &                YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
                  ENDIF
                ENDIF
!
!---            Gas enthalpy and internal energy  ---
!
                ISRX = 2
                CALL SP_B( TX,XLSX,PSBX )
                PVBX = PSBX
                DO IGC = 1, NGC
                  XMDGCX(IGC) = XMVGC(IGC,2,N)
                  XMDLCX(IGC) = XMVLC(IGC,2,N)
                ENDDO
                XMDGWX = XMVGW(2,N)
                XMDLWX = XMVLW(2,N)
                CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &            XGCX,XGWX,XLCZ,XLSX,XLWX,XMGCX,XMGWX,XMLCX,
     &            XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
                DO IGC = 1, NGC
                  XMVLC(IGC,2,N) = XMDLCX(IGC) 
                ENDDO
                XMVGW(2,N) = XMDGWX
                XMVLW(2,N) = XMDLWX
                XLAZ = 0.D+0
                DO IGC = 1, NGC
                 XLAZ = XLAZ + XLCZ(IGC)
                ENDDO
                PVAX = MIN( XLAX/XLAZ,1.D+0 )*PGAX
                WTMX = 0.D+0
                DO IGC = 1, NGC
                 WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
                ENDDO
                PVX = PVAX
                CALL ENTH_GC( TX,PVX,XMGCX,XMGWX,WTMX,HGAX,
     &            UGAX,ISRX )
!
!---            Aqueous component fractions and density  ---
!
                XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
                CALL DENS_L( TX,RHOBX,XLCX,RHOLX )
!
!---            Aqueous enthalpy and internal energy  ---
!
                ISRX = 1
                CALL ENTH_W( TX,PX,HLWX,HGWX,ISRX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
!
!---        Change in solubility with temperature ---
!
! values for source or node?
!                DO IGC = 1, NGC
!                  DHCX(IGC) = XMLCX(IGC)
!                END DO
                DTX = 1.D-6
                CALL EQUIL( TX+DTX,PX,PGAX,PGWX,PSBX,PVBX,
     &            XGCX,XGWX,XLCZ,XLSX,XLWX,XMGCX,XMGWX,XMLCX,
     &            XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
                DO IGC = 1, NGC
!                  IF( XMLCX(IGC).GT.EPSL ) THEN
!                    DHCX(IGC) = LOG(DHCX(IGC)/XMLCX(IGC))/DTX
!                  ELSE
                    DHCX(IGC) = 0.D+0
!                  ENDIF
                END DO
                CALL ENTH_L( TX,XLSX,XLCZ,DHCX,HBX,HGAX,HLX )
                DO IGC = 1, NGC
                  SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + SRX(4)*
     &              RHOLX*XLCX(IGC)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOLX*XLWX
                SRCS(M,N) = SRCS(M,N) + SRX(4)*RHOLX*XLSX
                SRCT(M,N) = SRCT(M,N) + SRX(4)*RHOLX*HLX
              ENDIF
!
!---        Gas volumetric rate w/ component relative humidity ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.4 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                TX = SRX(2)
                PX = MAX( PGX,PLX,SRX(3) )
                CALL SP_W( TX,PSWX )
                PVWX = PSWX*SRX(5)
                PVAX = MAX( PX-PVWX,0.D+0 )
                XMGAX = PVAX/PX
                XMGWX = PVWX/PX
                DO IGC = 1,NGC
                  XMGCX(IGC) = XMGAX*SRX(5+IGC)
                ENDDO
                sym = XMGWX*WTMW
                DO IGC = 1,NGC
                  sym = sym + XMGCX(IGC)*GCPP(1,IGC)
                ENDDO
                XGWX = XMGWX*WTMW/sym
                XGAX = 0.D+0
                DO IGC = 1,NGC
                  XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/sym
                  XGAX = XGAX + XGCX(IGC)
                ENDDO
!
!---            Gas density and component fractions  ---
!
                IF( XGWX.GT.EPSL ) THEN
                  ISRX = 2
                  CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
                ELSE
                  RHOGWX = 0.D+0
                ENDIF
                ISRX = 2
                CALL DENS_GC( TX,PVAX,XMGCX,RHOGAX,ISRX,ZGX )
                Z1X = ZGX
                DO IGC = 1,NGC
                  DO JGC = 1,NGC
                    XMGCIX(JGC) = 0.D+0
                    IF (IGC.EQ.JGC) XMGCIX(JGC) = 1.D+0
                  ENDDO
                  CALL DENS_GC( TX,PVAX,XMGCIX,RHOGCX(IGC),ISRX,ZGX )
                ENDDO
!
!---            Gas enthalpy and internal energy  ---
!
                WTMX = 0.D+0
                DO IGC = 1, NGC
                 WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
                ENDDO
                CALL ENTH_GC( TX,PVAX,XMGCX,XMGWX,WTMX,HGAX, 
     &              UGAX,ISRX)
!
!---            Water-vapor enthalpy  ---
!
                IF( XGWX.GT.EPSL ) THEN
                  ISRX = 2
                  CALL ENTH_W( TX,PVWX,HLWX,HGWX,ISRX )
                ELSE
                  HGWX = 0.D+0
                ENDIF
                SRCT(M,N) = SRCT(M,N) + SRX(4)*
     &            (RHOGWX*HGWX+RHOGAX*HGAX)
                DO IGC = 1, NGC
                  SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + 
     &              SRX(4)*RHOGCX(IGC)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOGWX
              ELSE
                HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                SRCT(M,N) = SRCT(M,N) + SRX(4)*RHOG(M,N)*HGX
                DO IGC = 1, NGC
                  SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + SRX(4)*
     &              RHOG(M,N)*XGC(IGC,M,N)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOG(M,N)*XGW(M,N)
              ENDIF
!
!---        Gas volumetric rate w/ component mass fractions ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.5 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                TX = SRX(2)
                PX = MAX( PGX,PLX,SRX(3) )
                XGWX = SRX(5)
                XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                DO IGC = 1,NGC
                  XGCX(IGC) = XGAX*SRX(5+IGC)
                ENDDO

                sym = XGWX/WTMW
                DO IGC = 1,NGC
                  sym = sym + XGCX(IGC)/GCPP(1,IGC)
                ENDDO
                XMGWX = XGWX/WTMW/sym
                DO IGC = 1,NGC
                  XMGCX(IGC) = XGCX(IGC)/GCPP(1,IGC)/sym
                ENDDO

                CALL SP_W( TX,PSWX )
                PVWX = MIN( PX*XMGWX,PSWX )
                PVAX = MAX( PX-PVWX,0.D+0 )
!
!---            Gas density and component fractions  ---
!
                IF( XGWX.GT.EPSL ) THEN
                  ISRX = 2
                  CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
                ELSE
                  RHOGWX = 0.D+0
                ENDIF
                ISRX = 2
                CALL DENS_GC( TX,PVAX,XMGCX,RHOGAX,ISRX,ZGX )
                Z1X = ZGX
                DO IGC = 1,NGC
                  DO JGC = 1,NGC
                    XMGCIX(JGC) = 0.D+0
                    IF (IGC.EQ.JGC) XMGCIX(JGC) = 1.D+0
                  ENDDO
                  CALL DENS_GC( TX,PVAX,XMGCIX,RHOGCX(IGC),ISRX,ZGX )
                ENDDO
!
!---            Gas enthalpy and internal energy  ---
!
                WTMX = 0.D+0
                DO IGC = 1, NGC
                 WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
                ENDDO
                CALL ENTH_GC( TX,PVAX,XMGCX,XMGWX,WTMX,HGAX, 
     &            UGAX,ISRX)
!
!---            Water-vapor enthalpy  ---
!
                IF( XGWX.GT.EPSL ) THEN
                  ISRX = 2
                  CALL ENTH_W( TX,PVWX,HLWX,HGWX,ISRX )
                ELSE
                  HGWX = 0.D+0
                ENDIF
                SRCT(M,N) = SRCT(M,N) + SRX(4)*
     &            (RHOGWX*HGWX+RHOGAX*HGAX)
                DO IGC = 1, NGC
                  SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + 
     &              SRX(4)*RHOGCX(IGC)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOGWX
              ELSE
                HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                SRCT(M,N) = SRCT(M,N) + SRX(4)*RHOG(M,N)*HGX
                DO IGC = 1, NGC
                  SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + SRX(4)*
     &              RHOG(M,N)*XGC(IGC,M,N)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOG(M,N)*XGW(M,N)
              ENDIF
!
!---        Salt source  ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.12 ) THEN
              SRCS(M,N) = SRCS(M,N) + SRX(4)
!
!---        Aqueous mass rate w/ dissolved-air relative
!           saturation, limit dissolved salt
!           sources to solubility limit  ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.7 ) THEN
              IF( SRX(4).LT.0.D+0 ) THEN
                DO IGC = 1, NGC
                  SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + SRX(4)*XLC(IGC,M,N)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XLW(M,N)
                SRCS(M,N) = SRCS(M,N) + SRX(4)*XLS(M,N)
                SRCT(M,N) = SRCT(M,N) + SRX(4)*HL(M,N)
              ELSE
                PORDX = 1.D+0
                SLX = 1.D+0
                TX = SRX(2)
                PX = MAX( PLX,PGX,SRX(3) )
                PCX = 0.D+0
!
!---            Dissolved-salt aqueous concentration  ---
!
                IF( MOD(ISRT(NS),1000)/100.EQ.1 ) THEN
!
!---              Dissolved-CO2 aqueous concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    RHOLSX = SRX(8)
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Conc. and ' //
     &                'Dissolved-CO2 Aqu. Conc.  @Node: '
                    CHMSGX(2) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_11( TX,PX,PCX,PGX,RHOLSX,RHOLAX,
!    &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Dissolved-CO2 aqueous relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    RHOLSX = SRX(8)
!                    PHILAX = SRX(5)
                    PHILAX = 0.D+0
                    DO IGC = 1, NGC
                      PHILCX(IGC) = SRX(8+IGC)
                      PHILAX = PHILAX + PHILCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
                    CHMSGX(2) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Conc. ' //
     &                'and Dissolved-CO2 Aqu. Rel. Sat.  @Node: '
!                   CALL FLSH_12( TX,PX,PCX,PGX,RHOLSX,PHILAX,
!    &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Dissolved-CO2 aqueous mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    RHOLSX = SRX(8)
                    XLAX = 0.D+0
                    DO IGC = 1, NGC
                      XLCX(IGC) = SRX(8+IGC)
                      XLAX = XLAX + XLCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Conc. and ' //
     &                'Dissolved-CO2 Aqu. Mass Frac.  @Node: '
                    CHMSGX(2) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_13( TX,PX,PCX,PGX,RHOLSX,
!    &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              No dissolved-CO2  ---
!
                  ELSE
                    RHOLSX = SRX(8)
!                    DO IGC = 1, NGC
!                      XLCX(IGC) = 0.D+0
!                    ENDDO
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Conc. and ' //
     &                'Dissolved-CO2 Aqu. Mass Frac.  @Node: '
                    CHMSGX(2) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_13( TX,PX,PCX,PGX,RHOLSX,
!    &                YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
!
!---            Dissolved-salt aqueous relative saturation  ---
!
                ELSEIF( MOD(ISRT(NS),1000)/100.EQ.2 ) THEN
!
!---              Dissolved-CO2 aqueous concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    PHILSX = SRX(8)
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Rel. Sat. and ' //
     &                'Dissolved-CO2 Aqu. Conc.  @Node: '
                    CHMSGX(2) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_21( TX,PX,PCX,PGX,PHILSX,RHOLAX,
!    &                YLSX,XLSX,XLCX,CHMSGX )
!
!---              Dissolved-CO2 aqueous relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    PHILSX = SRX(8)
!                    PHILAX = SRX(5)
                    PHILAX = 0.D+0
                    DO IGC = 1, NGC
                      PHILCX(IGC) = SRX(8+IGC)
                      PHILAX = PHILAX + PHILCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_22( TX,PX,PCX,PGX,PHILSX,PHILAX,
!    &                YLSX,XLSX,XLCX,CHMSGX )
!
!---              Dissolved-CO2 aqueous mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    PHILSX = SRX(8)
                    DO IGC = 1, NGC
                      XLCX(IGC) = SRX(8+IGC)
                    ENDDO
!                    XLAX = SRX(5)
                    CHMSGX(1) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_23( TX,PX,PCX,PGX,PHILSX,
!    &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              No dissolved-CO2  ---
!
                  ELSE
                    PHILSX = SRX(8)
!                    DO IGC = 1, NGC
!                      XLCX(IGC) = 0.D+0
!                    ENDDO
                   XLAX = 0.D+0
                    CHMSGX(1) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_23( TX,PX,PCX,PGX,PHILSX,
!    &                YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
!
!---            Dissolved-salt aqueous mass fraction  ---
!
                ELSEIF( MOD(ISRT(NS),1000)/100.EQ.3 ) THEN
!
!---              Dissolved-CO2 aqueous concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    YLSX = SRX(8)
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Mass Frac. and ' //
     &                'Dissolved-CO2 Aqu. Conc.  @Node: '
                    CHMSGX(2) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_31( TX,PX,PCX,PGX,RHOLAX,
!    &                YLSX,XLSX,XLCX,CHMSGX )
!
!---              Dissolved-CO2 aqueous relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    YLSX = SRX(8)
!                    PHILAX = SRX(5)
                    PHILAX = 0.D+0
                    DO IGC = 1, NGC
                      PHILCX(IGC) = SRX(8+IGC)
                      PHILAX = PHILAX + PHILCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_32( TX,PX,PCX,PGX,PHILAX,
!    &                YLSX,XLSX,XLCX,CHMSGX )
!
!---              Dissolved-CO2 aqueous mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    YLSX = SRX(8)
                    XLAX = 0.D+0
                    DO IGC = 1, NGC
                      XLCX(IGC) = SRX(8+IGC)
                      XLAX = XLAX + XLCX(IGC)
                    ENDDO
!                    XLAX = SRX(5)
                    CHMSGX(1) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_33( TX,PX,PCX,PGX,
!    &                YLSX,XLSX,XLCX,CHMSGX )
!
!---              No dissolved-CO2  ---
!
                  ELSE
                    YLSX = SRX(8)
!~                    DO IGC = 1, NGC
!@                      XLCX(IGC) = 0.D+0
!                    ENDDO
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_33( TX,PX,PCX,PGX,
!    &                YLSX,XLSX,XLCX,CHMSGX )
                  ENDIF
!
!---            No dissolved-salt  ---
!
                ELSE
!
!---              Dissolved-CO2 aqueous concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    YLSX = 0.D+0
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Dissolved-Salt Aqu. Mass Frac. and ' //
     &                'Dissolved-CO2 Aqu. Conc.  @Node: '
                    CHMSGX(2) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_31( TX,PX,PCX,PGX,RHOLAX,
!    &                YLSX,XLSX,XLCX,CHMSGX )
!
!---              Dissolved-CO2 aqueous relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    YLSX = 0.D+0
!                    PHILAX = SRX(5)
                    PHILAX = 0.D+0
                    DO IGC = 1, NGC
                     PHILCX = SRX(5)
                     PHILAX = PHILAX + PHILCX(IGC)
                    ENDDO
                    CHMSGX(1) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_32( TX,PX,PCX,PGX,PHILAX,
!    &                YLSX,XLSX,XLCX,CHMSGX )
!
!---              Dissolved-CO2 aqueous mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    YLSX = 0.D+0
                    XLAX = 0.D+0
                    DO IGC = 1, NGC
                      XLCX(IGC) = SRX(5+IGC)
                      XLAX = XLAX + XLCX(IGC)
                    ENDDO
!                    XLAX = SRX(5)
                    CHMSGX(1) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_33( TX,PX,PCX,PGX,
!    &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              No dissolved-CO2  ---
!
                  ELSE
                    YLSX = 0.D+0
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Aqu. Mass Source: Transition ' //
     &                'to Unsaturated Conditions @ Node'
!                   CALL FLSH_33( TX,PX,PCX,PGX,
!    &                YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
                ENDIF
!
!---            CO2 enthalpy and internal energy  ---
!
                ISRX = 2
                CALL SP_B( TX,XLSX,PSBX )
                PVBX = PSBX
                DO IGC = 1, NGC
                  XMDGCX(IGC) = XMVGC(IGC,2,N)
                  XMDLCX(IGC) = XMVLC(IGC,2,N)
                ENDDO
                XMDGWX = XMVGW(2,N)
                XMDLWX = XMVLW(2,N)
                CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &            XGCX,XGWX,XLCZ,XLSX,XLWX,XMGCX,XMGWX,
     &            XMLCX,XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
                DO IGC = 1, NGC
                  XMVLC(IGC,2,N) = XMDLCX(IGC) 
                ENDDO
                XMVGW(2,N) = XMDGWX
                XMVLW(2,N) = XMDLWX
                XLAZ = 0.D+0
                DO IGC = 1, NGC
                 XLAZ = XLAZ + XLCZ(IGC)
                ENDDO
                PVAX = MIN( XLAX/XLAZ,1.D+0 )*PGAX
                WTMX = 0.D+0
                DO IGC = 1, NGC
                 WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
                ENDDO
                PVX = PVAX
                 CALL ENTH_GC( TX,PVX,XMGCX,XMGWX,WTMX,HGAX,
     &             UGAX,ISRX )
!
!---            Aqueous component fraction  ---
!
                XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
!
!---            Aqueous enthalpy and internal energy  ---
!
                ISRX = 1
                CALL ENTH_W( TX,PX,HLWX,HGWX,ISRX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
!
!---        Change in solubility with temperature ---
!
! values for source or node?
!                DO IGC = 1, NGC
!                  DHCX(IGC) = XMLCX(IGC)
!                END DO
                DTX = 1.D-6
                CALL EQUIL( TX+DTX,PX,PGAX,PGWX,PSBX,PVBX,
     &            XGCX,XGWX,XLCZ,XLSX,XLWX,XMGCX,XMGWX,XMLCX,
     &            XMLSX,XMLWX,ZGX,WTMGCX,XMDGCX,XMDGWX,
     &            XMDLCX,XMDLWX)
                DO IGC = 1, NGC
!                  IF( XMLCX(IGC).GT.EPSL ) THEN
!                    DHCX(IGC) = LOG(DHCX(IGC)/XMLCX(IGC))/DTX
!                  ELSE
                    DHCX(IGC) = 0.D+0
!                  ENDIF
                END DO
                CALL ENTH_L( TX,XLSX,XLCZ,DHCX,HBX,HGAX,HLX )
                DO IGC = 1, NGC
                  SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + SRX(4)*XLCX(IGC)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XLWX
                SRCS(M,N) = SRCS(M,N) + SRX(4)*XLSX
                SRCT(M,N) = SRCT(M,N) + SRX(4)*HLX
              ENDIF
!
!---        Gas mass rate w/ component relative humidity ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.8 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                TX = SRX(2)
                PX = MAX( PGX,PLX,SRX(3) )
                CALL SP_W( TX,PSWX )
                PVWX = PSWX*SRX(5)
                PVAX = MAX( PX-PVWX,0.D+0 )
                XMGAX = PVAX/PX
                XMGWX = PVWX/PX
                DO IGC = 1,NGC
                  XMGCX(IGC) = XMGAX*SRX(5+IGC)
                ENDDO
                sym = XMGWX*WTMW
                DO IGC = 1,NGC
                  sym = sym + XMGCX(IGC)*GCPP(1,IGC)
                ENDDO
                XGWX = XMGWX*WTMW/sym
                XGAX = 0.D+0
                DO IGC = 1,NGC
                  XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/sym
                  XGAX = XGAX + XGCX(IGC)
                ENDDO
!
!---            Gas density and component fractions  ---
!
                IF( PVWX.GT.EPSL ) THEN
                  ISRX = 2
                  CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
                ELSE
                  RHOGWX = 0.D+0
                ENDIF
                ISRX = 2
                CALL DENS_GC( TX,PVAX,XMGCX,RHOGAX,ISRX,ZGX )
                Z1X = ZGX
!
!---            Gas enthalpy and internal energy  ---
!
                WTMX = 0.D+0
                DO IGC = 1, NGC
                 WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
                ENDDO
                CALL ENTH_GC( TX,PVAX,XMGCX,XMGWX,WTMX,HGAX,
     &            UGAX,ISRX  )
!
!---            Water-vapor enthalpy  ---
!
                IF( XGWX.GT.EPSL ) THEN
                  ISRX = 2
                  CALL ENTH_W( TX,PVWX,HLWX,HGWX,ISRX )
                ELSE
                  HGWX = 0.D+0
                ENDIF
                SRCT(M,N) = SRCT(M,N) + SRX(4)*(XGWX*HGWX+XGAX*HGAX)
                DO IGC = 1, NGC
                   SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + SRX(4)*XGCX(IGC)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGWX
              ELSE
                HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                SRCT(M,N) = SRCT(M,N) + SRX(4)*HGX
                DO IGC = 1, NGC
                  SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + SRX(4)*XGC(IGC,M,N)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGW(M,N)
              ENDIF
!
!---        Gas mass rate w/ component mass fractions ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.9 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                TX = SRX(2)
                PX = MAX( PGX,PLX,SRX(3) )
                XGWX = SRX(5)
                XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                DO IGC = 1,NGC
                  XGCX(IGC) = XGAX*SRX(5+IGC)
                ENDDO

                sym = XGWX/WTMW
                DO IGC = 1,NGC
                  sym = sym + XGCX(IGC)/GCPP(1,IGC)
                ENDDO
                XMGWX = XGWX/WTMW/sym
                DO IGC = 1,NGC
                  XMGCX(IGC) = XGCX(IGC)/GCPP(1,IGC)/sym
                ENDDO
                
                CALL SP_W( TX,PSWX )
                PVWX = MIN( PX*XMGWX,PSWX )
                PVAX = MAX( PX-PVWX,0.D+0 )
!
!---            Gas density and component fractions  ---
!
                IF( XGWX.GT.EPSL ) THEN
                  ISRX = 2
                  CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
                ELSE
                  RHOGWX = 0.D+0
                ENDIF
                ISRX = 2
                CALL DENS_GC( TX,PVAX,XMGCX,RHOGAX,ISRX,ZGX )
                Z1X = ZGX
!
!---            Gas enthalpy and internal energy  ---
!
                WTMX = 0.D+0
                DO IGC = 1, NGC
                 WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
                ENDDO
                CALL ENTH_GC( TX,PVAX,XMGCX,XMGWX,WTMX,HGAX,
     &               UGAX,ISRX)
!
!---            Water-vapor enthalpy  ---
!
                IF( XGWX.GT.EPSL ) THEN
                  ISRX = 2
                  CALL ENTH_W( TX,PVWX,HLWX,HGWX,ISRX )
                ELSE
                  HGWX = 0.D+0
                ENDIF
                SRCT(M,N) = SRCT(M,N) + SRX(4)*(XGWX*HGWX+XGAX*HGAX)
                DO IGC = 1, NGC
                  SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + SRX(4)*XGCX(IGC)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGWX
              ELSE
                HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                SRCT(M,N) = SRCT(M,N) + SRX(4)*HGX
                DO IGC = 1, NGC
                  SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + SRX(4)*XGC(IGC,M,N)
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGW(M,N)
              ENDIF
!
!---        Salt density source  ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.11 ) THEN
              SRCS(M,N) = SRCS(M,N) + SRX(4)*VOL(N)
!
!---        Z-direction vertical injection well  ---
!
            ELSEIF( MOD(ISRT(NS),100).GE.13 .AND. 
     &        MOD(ISRT(NS),100).LE.15  ) THEN
!
!---          Geometric factors  ---
!
              RDW = SRX(3)
              RDE = SQRT( AFZ(NSZ(N))/GPI/SRX(4) )
              DRD2 = (RDE**2-RDW**2)

!
!---          Well pressure  ---
!
              IF( M.EQ.2 ) THEN
                IF( K.EQ.ISRDM(5,NS) ) THEN
                  PGWX = SRX(2)
                  PX = PGWX+PATM
                  CALL SP_W( SRX(7),PSWX )

                  IF( MOD(ISRT(NS),100).EQ.13 ) THEN
                    PVWX = PSWX*SRX(5)
                    PVAX = MAX( PX-PVWX,0.D+0 )
                    XMGAX = PVAX/PX
                    XMGWX = PVWX/PX
                    DO IGC = 1,NGC
                      XMGCX(IGC) = XMGAX*SRX(7+IGC)
                    ENDDO
                    sym = XMGWX*WTMW
                    DO IGC = 1,NGC
                      sym = sym + XMGCX(IGC)*GCPP(1,IGC)
                    ENDDO
                    XGWX = XMGWX*WTMW/sym
                    XGAX = 0.D+0
                    DO IGC = 1,NGC
                      XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/sym
                      XGAX = XGAX + XGCX(IGC)
                    ENDDO
                  ELSEIF( MOD(ISRT(NS),100).EQ.14 ) THEN
                    XGWX = SRX(5)
                    XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                    DO IGC = 1,NGC
                      XGCX(IGC) = XGAX*SRX(7+IGC)
                    ENDDO
                    sym = XGWX/WTMW
                    DO IGC = 1,NGC
                      sym = sym + XGCX(IGC)/GCPP(1,IGC)
                    ENDDO
                    XMGWX = XGWX/WTMW/sym
                    DO IGC = 1,NGC
                      XMGCX(IGC) = XGCX(IGC)/GCPP(1,IGC)/sym
                    ENDDO
                    PVWX = MIN( PX*XMGWX,PSWX )
                    PVAX = MAX( PX-PVWX,0.D+0 )
                  ENDIF
!
!---              Gas density and component fractions  ---
!
                  IF( XGWX.GT.EPSL ) THEN
                    ISRX = 2
                    CALL DENS_W( SRX(7),PVWX,RHOLWX,RHOGWX,ISRX )
                  ELSE
                    RHOGWX = 0.D+0
                  ENDIF
                  ISRX = 2
                  CALL DENS_GC( SRX(7),PVAX,XMGCX,RHOGAX,ISRX,ZGX )
!
!---              Gas enthalpy and internal energy  ---
!
                  WTMX = 0.D+0
                  DO IGC = 1, NGC
                   WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
                  ENDDO
                  CALL ENTH_GC( SRX(7),PVAX,XMGCX,XMGWX,WTMX,HGAX,
     &                 UGAX,ISRX)
!
!---              Water-vapor enthalpy  ---
!
                  IF( XGWX.GT.EPSL ) THEN
                    ISRX = 2
                    CALL ENTH_W( SRX(7),PVWX,HLWX,HGWX,ISRX )
                  ELSE
                    HGWX = 0.D+0
                  ENDIF
                ELSE
                  PX = PGWX+PATM
                  CALL SP_W( SRX(7),PSWX )
                  IF( MOD(ISRT(NS),100).EQ.13 ) THEN
                    PVWX = PSWX*SRX(5)
                    PVAX = MAX( PX-PVWX,0.D+0 )
                    XMGAX = PVAX/PX
                    XMGWX = PVWX/PX
                    DO IGC = 1,NGC
                      XMGCX(IGC) = XMGAX*SRX(7+IGC)
                    ENDDO
                    sym = XMGWX*WTMW
                    DO IGC = 1,NGC
                      sym = sym + XMGCX(IGC)*GCPP(1,IGC)
                    ENDDO
                    XGWX = XMGWX*WTMW/sym
                    XGAX = 0.D+0
                    DO IGC = 1,NGC
                      XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/sym
                      XGAX = XGAX + XGCX(IGC)
                    ENDDO
                  ELSEIF( MOD(ISRT(NS),100).EQ.14 ) THEN
                    XGWX = SRX(5)
                    XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                    DO IGC = 1,NGC
                      XGCX(IGC) = XGAX*SRX(7+IGC)
                    ENDDO
                    sym = XGWX/WTMW
                    DO IGC = 1,NGC
                      sym = sym + XGCX(IGC)/GCPP(1,IGC)
                    ENDDO
                    XMGWX = XGWX/WTMW/sym
                    DO IGC = 1,NGC
                      XMGCX(IGC) = XGCX(IGC)/GCPP(1,IGC)/sym
                    ENDDO
                    PVWX = MIN( PX*XMGWX,PSWX )
                    PVAX = MAX( PX-PVWX,0.D+0 )
                  ENDIF
!
!---              Gas density and component fractions  ---
!
                  IF( XGWX.GT.EPSL ) THEN
                    ISRX = 2
                    CALL DENS_W( SRX(7),PVWX,RHOLWX,RHOGWX,ISRX )
                  ELSE
                    RHOGWX = 0.D+0
                  ENDIF
                  ISRX = 2
                  CALL DENS_GC( SRX(7),PVAX,XMGCX,RHOGAX,ISRX,ZGX )
!
!---              Gas enthalpy and internal energy  ---
!
                  WTMX = 0.D+0
                  DO IGC = 1, NGC
                   WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
                  ENDDO
                  CALL ENTH_GC( SRX(7),PVAX,XMGCX,XMGWX,WTMX,HGAX,
     &                 UGAX,ISRX)
!
!---              Water-vapor enthalpy  ---
!
                  IF( XGWX.GT.EPSL ) THEN
                    ISRX = 2
                    CALL ENTH_W( SRX(7),PVWX,HLWX,HGWX,ISRX )
                  ELSE
                    HGWX = 0.D+0
                  ENDIF
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
!---          If permeability is too low, skip
!     
              IF ( PERMX.GT.SRX(6) ) THEN
!
!---            Injection  ---
!
                IF( PGWX-PG(M,N).GT.EPSL ) THEN
!
!---              Gas viscosity  ---
!
                  DO IGC = 1, NGC
                    XMGCX(IGC)=XMGC(IGC,M,N)
                  ENDDO
                  CALL VISC_GC( SRX(7),PX,ZGX,XMGCX,RHOGAX,VISCGX,
     &              VISGCX )
                  CALL VISC_W( SRX(7),PVW(M,N),RHOGWX,VISGWX )
                  CALL VISC_G( VISCGX,VISGWX,XMGCX,XMGW(M,N),VISGX )
                  VISGX=VISGX+VISGCX
                  HCGX = 2.D+0*GPI*PERMX*DRD2*DZGF(N)/
     &              (VISGX*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
                  QGX = SRX(4)*(PGWX-PG(M,N))*HCGX
                  SRCT(M,N) = SRCT(M,N) + QGX*(RHOGWX*HGWX+RHOGAX
     &              *HGAX)
                  SRCW(M,N) = SRCW(M,N) + QGX*RHOGWX
                  DO IGC = 1, NGC
                    SRCGC(IGC,M,N) = SRCGC(IGC,M,N) 
     &                + QGX*RHOGAX*XGCX(IGC)
                  ENDDO
!
!---            Withdrawl  ---
!
                ELSE
                  HCGX = 2.D+0*GPI*PERMX*DRD2*DZGF(N)/
     &              (VISG(M,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
                  QGX = SRX(4)*(PGWX-PG(M,N))*RKG(M,N)*HCGX 
                  HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                  SRCT(M,N) = SRCT(M,N) + QGX*RHOG(M,N)*HGX
                  SRCW(M,N) = SRCW(M,N) + QGX*RHOG(M,N)*XGW(M,N)
                  DO IGC = 1, NGC
                    SRCGC(IGC,M,N) = SRCGC(IGC,M,N) 
     &                + QGX*RHOG(M,N)*XGC(IGC,M,N)
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
  400     CONTINUE
  500   CONTINUE
!
!---  Mass injection wells  ---
!
  510   CONTINUE
!
!---    Z-direction Gas Mass Rate injection well  ---
!
        IF( MOD(ISRT(NS),100).EQ.23 .OR. MOD(ISRT(NS),100).EQ.24 ) THEN
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
!---          Gas Mass Rate w/ Component Relative Humidity ---
!
              IF( MOD(ISRT(NS),100).EQ.23 ) THEN
                IF( SRX(4).GE.0.D+0 ) THEN
                  TX = SRX(2)
                  PGWX = SRX(4)*WTX/SRX(7)/HCGX+PG(M,N)
                  PX = MAX( PGX,PLX,PGWX+PATM )
                  CALL SP_W( TX,PSWX )
                  PVWX = PSWX*SRX(5)
                  PVAX = MAX( PX-PVWX,0.D+0 )
                  XMGAX = PVAX/PX
                  XMGWX = PVWX/PX
                  DO IGC = 1,NGC
                    XMGCX(IGC) = XMGAX*SRX(7+IGC)
                  ENDDO
                  sym = XMGWX*WTMW
                  DO IGC = 1,NGC
                    sym = sym + XMGCX(IGC)*GCPP(1,IGC)
                  ENDDO
                  XGWX = XMGWX*WTMW/sym
          XGAX = 0.D+0
                  DO IGC = 1,NGC
                    XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/sym
                    XGAX = XGAX + XGCX(IGC)
          ENDDO
!
!---              Gas density and component fractions  ---
!
                  IF( XGWX.GT.EPSL ) THEN
                    ISRX = 2
                    CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
          ELSE
                    RHOGWX = 0.D+0
          ENDIF
          ISRX = 2
          CALL DENS_GC( TX,PVAX,XMGCX,RHOGAX,ISRX,ZGX )
                Z1X = ZGX
                  DO IGC = 1,NGC
                    DO JGC = 1,NGC
                      XMGCIX(JGC) = 0.D+0
                      IF (IGC.EQ.JGC) XMGCIX(JGC) = 1.D+0
                    ENDDO
                    CALL DENS_GC( TX,PVAX,XMGCIX,RHOGCX(IGC),ISRX,ZGX )
                  ENDDO
!
!---              Gas enthalpy and internal energy  ---
!
              WTMX = 0.D+0
              DO IGC = 1, NGC
                WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
              ENDDO
                  CALL ENTH_GC( TX,PVAX,XMGCX,XMGWX,WTMX,HGAX, 
     &                UGAX,ISRX)
!
!---              Water-vapor enthalpy  ---
!
                  IF( XGWX.GT.EPSL ) THEN
                    ISRX = 2
                    CALL ENTH_W( TX,PVWX,HLWX,HGWX,ISRX )
              ELSE
                    HGWX = 0.D+0
              ENDIF
                  
                SRCT(M,N) = SRCT(M,N) + SRX(4)*(XGWX*HGWX+XGAX*HGAX)
     &                                  *WTX
                DO IGC = 1, NGC
                   SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + SRX(4)*XGCX(IGC)
     &                                               *WTX
                ENDDO
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGWX*WTX
                
                ELSE
!                 DPGX = MIN( SRX(3)-PG(M,N),0.D+0 )
!                 QGX = MAX( DPGX*RKG(M,N)*HCGX,WTX*SRX(4) )
                  QGX = SRX(4)*WTX
                  HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                  SRCT(M,N) = SRCT(M,N) + QGX*HGX
                  DO IGC = 1, NGC
                    SRCGC(IGC,M,N) = SRCGC(IGC,M,N) + QGX*XGC(IGC,M,N)
                  ENDDO
                  SRCW(M,N) = SRCW(M,N) + QGX*XGW(M,N)
            ENDIF
          ENDIF
!
!---  Gas Mass Rate w/ Component Mass Fractions ---
!
              IF( MOD(ISRT(NS),100).EQ.24 .OR. 
     &          MOD(ISRT(NS),100).EQ.25 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    PGWX = SRX(4)*WTX/SRX(7)/HCGX+PG(M,N)
                    PX = MAX( PGX,PLX,PGWX+PATM )
                    TX = SRX(2)
                    XGWX = SRX(5)
                    XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                    DO IGC = 1,NGC
                      XGCX(IGC) = XGAX*SRX(7+IGC)
                    ENDDO

                    sym = XGWX/WTMW
                    DO IGC = 1,NGC
                      sym = sym + XGCX(IGC)/GCPP(1,IGC)
                    ENDDO
                    XMGWX = XGWX/WTMW/sym
                    DO IGC = 1,NGC
                      XMGCX(IGC) = XGCX(IGC)/GCPP(1,IGC)/sym
                    ENDDO

                    CALL SP_W( TX,PSWX )
                    PVWX = MIN( PX*XMGWX,PSWX )
                    PVAX = MAX( PX-PVWX,0.D+0 )
!
!---                Gas density and component fractions  ---
!
                    IF( XGWX.GT.EPSL ) THEN
                      ISRX = 2
                      CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
            ELSE
                      RHOGWX = 0.D+0
            ENDIF
                    ISRX = 2
                    CALL DENS_GC( TX,PVAX,XMGCX,RHOGAX,ISRX,ZGX )
                    Z1X = ZGX
                    DO IGC = 1,NGC
                      DO JGC = 1,NGC
                        XMGCIX(JGC) = 0.D+0
                        IF (IGC.EQ.JGC) XMGCIX(JGC) = 1.D+0
                      ENDDO
                      CALL DENS_GC( TX,PVAX,XMGCIX,RHOGCX(IGC),ISRX,
     &                  ZGX )
                    ENDDO
!
!---                Gas enthalpy and internal energy  ---
!
                    WTMX = 0.D+0
                    DO IGC = 1, NGC
                     WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
                    ENDDO
                    CALL ENTH_GC( TX,PVAX,XMGCX,XMGWX,WTMX,HGAX, 
     &                UGAX,ISRX)
!
!---                Water-vapor enthalpy  ---
!
                    IF( XGWX.GT.EPSL ) THEN
                      ISRX = 2
                      CALL ENTH_W( TX,PVWX,HLWX,HGWX,ISRX )
                    ELSE
                      HGWX = 0.D+0
                    ENDIF
                
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*
     &                (XGWX*HGWX+XGAX*HGAX)*WTX
                    DO IGC = 1, NGC
                      SRCGC(IGC,M,N) = SRCGC(IGC,M,N)
     &                  + SRX(4)*XGCX(IGC)*WTX
                    ENDDO
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*XGWX*WTX

                  ELSE
!                    DPGX = MIN( SRX(3)-PG(M,N),0.D+0 )
!                    QGX = MAX( DPGX*RKG(M,N)*HCGX,WTX*SRX(4) )
                    QGX = SRX(4)*WTX
                    HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                    SRCT(M,N) = SRCT(M,N) + QGX*HGX
                    DO IGC = 1, NGC
                      SRCGC(IGC,M,N) = SRCGC(IGC,M,N) 
     &                + QGX*XGC(IGC,M,N)
                    ENDDO
                    SRCW(M,N) = SRCW(M,N) + QGX*XGW(M,N)
                  ENDIF
                ENDIF
  516       CONTINUE
  518     CONTINUE
        ENDIF
  900 CONTINUE

      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORC_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORIC_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Compute coupled-equation source integrals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE GRID
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
      SUB_LOG(ISUB_LOG) = '/SORIC_COMP'
!
!---  Loop over sources  ---
!
      DO 900 NS = 1,NSR
        IF( TM.LE.SRC(1,1,NS) ) GOTO 900
        SRX(1) = TM
        IF( ISRM(NS).EQ.1 ) THEN
          DO 70 N = 1,8+NSOLU
            SRX(N) = SRC(N,1,NS)
   70     CONTINUE
        ELSE
          DO 100 M = 2,ISRM(NS)
            IF( TM.LE.SRC(1,M,NS) ) THEN
             DTSR = MIN( SRC(1,M,NS)-TM,DT )
             TFSR = (TM-0.5D+0*DTSR-SRC(1,M-1,NS))/
     &         (SRC(1,M,NS)-SRC(1,M-1,NS))
             DO 80 N = 1,8+NSOLU
               SRX(N) = SRC(N,M-1,NS) + TFSR*(SRC(N,M,NS)-SRC(N,M-1,NS))
   80        CONTINUE
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 900
        ENDIF
  110   CONTINUE
!
!---    Integrate total CO2 and water mass injected  ---
!
        IF( MOD(ISRT(NS),100).GE.25 .AND. 
     &    MOD(ISRT(NS),100).LE.27  ) THEN
          DO 600 NC = 1,ISRDM(1,NS)
            N = IWSI(NC,NS)
            IF( IXP(N).LE.0 ) GOTO 600
            SRCP(1,NS) = SRCP(1,NS) + SWSI(1,NC,NS)*SWSI(2,NC,NS)*DT
            SRCP(2,NS) = SRCP(2,NS) + SWSI(1,NC,NS)*SWSI(3,NC,NS)*DT
  600     CONTINUE
        ENDIF
  900 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORIC_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORIT_COMP( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Compute solute transport source integrals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/SORIT_COMP'
!
!---  Loop over sources  ---
!
      DO 600 NS = 1,NSR
        IF( TM.LE.SRC(1,1,NS) ) GOTO 600
        SRX(1) = TM
        IF( ISRM(NS).EQ.1 ) THEN
          DO 70 N = 2,8+NSOLU
            SRX(N) = SRC(N,1,NS)
   70     CONTINUE
        ELSE
          DO 100 M = 2,ISRM(NS)
            IF( TM.LE.SRC(1,M,NS) ) THEN
             DTSR = MIN( SRC(1,M,NS)-TM,DT )
             TFSR = (TM-0.5D+0*DTSR-SRC(1,M-1,NS))/
     &         (SRC(1,M,NS)-SRC(1,M-1,NS))
             DO 80 N = 1,8+NSOLU
               SRX(N) = SRC(N,M-1,NS) + TFSR*(SRC(N,M,NS)-SRC(N,M-1,NS))
   80        CONTINUE
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
              IF( IXP(N).LE.0 ) GOTO 500
!
!---  Aqueous Volumetric Sink  ---
!
              IF( MOD(ISRT(NS),100).EQ.3 .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YL(N,NSL)*DT/(PORD(2,N)*SL(2,N))
!
!---  Gas Volumetric Sink  ---
!
              ELSEIF( (MOD(ISRT(NS),100).EQ.4 .OR. 
     &          MOD(ISRT(NS),100).EQ.5) .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YG(N,NSL)*DT/(PORD(2,N)*SG(2,N))
!
!---  Aqueous Mass Sink  ---
!
              ELSEIF( MOD(ISRT(NS),100).EQ.7 .AND. 
     &          SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YL(N,NSL)*DT/(RHOL(2,N)*PORD(2,N)*SL(2,N))
!
!---  Gas Mass Sink  ---
!
              ELSEIF( (MOD(ISRT(NS),100).EQ.8 .OR. 
     &          MOD(ISRT(NS),100).EQ.9) .AND. SRX(4).LT.0.D+0 ) THEN
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
!---  End of SORIT_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT_COMP( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Compute solute transport source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/SORT_COMP'
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
          IF( IXP(N).LE.0 ) GOTO 500
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
          IF( MOD(ISRT(NS),100).EQ.3 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YL(N,NSL)/(PORD(2,N)*SL(2,N))
!
!---      Gas Volumetric Sink  ---
!
          ELSEIF( MOD(ISRT(NS),100).EQ.4 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YG(N,NSL)/(PORD(2,N)*SG(2,N))
!
!---      Gas Volumetric Sink  ---
!
          ELSEIF( MOD(ISRT(NS),100).EQ.5 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YG(N,NSL)/(PORD(2,N)*SG(2,N))
!
!---      Aqueous Mass Sink  ---
!
          ELSEIF( MOD(ISRT(NS),100).EQ.7 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YL(N,NSL)/(RHOL(2,N)*PORD(2,N)*SL(2,N))
!
!---      Gas Mass Sink  ---
!
          ELSEIF( MOD(ISRT(NS),100).EQ.8 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YG(N,NSL)/(RHOG(2,N)*PORD(2,N)*SG(2,N))
!
!---      Gas Mass Sink  ---
!
          ELSEIF( MOD(ISRT(NS),100).EQ.9 .AND. SRX(4).LT.0.D+0 ) THEN
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
!---  End of SORT_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPRP_COMP( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Calculates the aqueous- and gas-phase solute
!     mole fractions from user-specified partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/SPRP_COMP'
!
!---  Loop over all nodes  ---
!
      DO 900 N = 1,NFLD
        N_DB = N
        IF( IXP(N).LE.0 ) GOTO 900
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
!---  End of SPRP_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TPORT_COMP( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Solute/Reactive Species Transport Shell.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!
#ifdef lis
!
!----------------------PETSc Modules-----------------------------------!
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
      USE SOLTN
      USE JACOB
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TPORT_COMP'
!
!---  Zero Jacobian matrix  ---
!
#ifdef petsc
      IF( ILES.EQ.5 ) CALL MatZeroEntries(T_MAT,IERR)
#endif
      INDX = 1
      CALL JCBZ( ISVT,MUT,MLT,MKT,INDX )
!
!---  Compute solute sources ---
!
      CALL SORT_COMP( NSL )
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
      CALL SBND_COMP( NSL )

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
!---  Update solute concentrations ---
!
      CALL UPDTC( NSL )
!
!---  Compute solute aqueous-phase fluxes (interior nodes)  ---
!
      CALL SFXL( NSL )
!
!---  Compute solute gas-phase fluxes (interior nodes)  ---
!
      CALL SFXG( NSL )
!
!---  Compute solute aqueous and gas fluxes (boundary surfaces)  ---
!
      CALL SFXB32( NSL )
!
!---  Integrate solute sources  ---
!
      CALL SORIT_COMP( NSL )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TPORT_COMP group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_COMP
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     H2O-NaCl-CO2-Energy Mode
!
!     Update the primary variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE OUTPU
      USE JACOB
      USE HYST
      USE GRID
      USE FILES
      USE FDVS
      USE FDVP
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
      CHARACTER*64 PH_CND(3)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Saturated w/ Dissolved CO2',
     &             'Unsaturated', 'Saturated w/ Trapped Gas'/
      REAL*8 DPC(LNGC)
      INTEGER MPC(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_COMP'

      IF( ICNV.EQ.1 ) GOTO 300
      IERR = 0
!
!---  Update primary variables
!
      DO 200 N = 1,NFLD
        N_DB = N
        IF( IXP(N).LE.0 ) GOTO 200
        NMD = IXP(N)
        MPL = IM(IEQW,NMD)
!        MPG = IM(IEQA,NMD)
        DPL = BLU(MPL)
        DO IGC = 1, NGC
          MPC(IGC) = IM(IEQGC(IGC),NMD)
          DPC(IGC) = BLU(MPC(IGC))
        ENDDO
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
!
!---      Limit changes in temperature to 0.5 C, 
!         isothermal option   ---
!
          MPT = IM(IEQT,NMD)
          DPT = BLU(MPT)
          DPT = SIGN( MIN( 5.D-1,ABS(DPT) ),DPT )
          T(2,N) = MAX( (T(2,N)+DPT),1.D-1 )
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          MPS = IM(IEQS,NMD)
          DPS = BLU(MPS)
        ENDIF
!
!---    Saturated system w/o entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - dissolved-CO2 aqueous mass fraction
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
!
!---      Limit changes in pressure  ---
!
          DPX = 2.5D-1*MAX(PL(2,N)+PATM,1.D+6)
          DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
          PL(2,N) = PL(2,N) + DPL
!
!---      Zero negative corrections for zero dissolved CO2  ---
!
          DO IGC = 1, NGC
           IF( XLC(IGC,2,N)/EPSL.LT.EPSL .AND. BLU(MPC(IGC))/EPSL.LT.
     &             EPSL ) THEN
            BLU(MPC(IGC)) = 0.D+0
            DPC(IGC) = 0.D+0
           ENDIF
          ENDDO
          XLA(2,N) = 0.D+0
          DO IGC = 1, NGC
            XLC(IGC,2,N) = MAX( (XLC(IGC,2,N)+DPC(IGC)),0.D+0 )
            XLA(2,N) = XLA(2,N) + XLC(IGC,2,N)
          ENDDO
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Limit salt mass fraction changes to 0.25 of the
!           maximum value if salt mass fraction is less than
!           the maximum   ---
!
            CALL SOL_LS( T(2,N),XLSMX )
            IF( YLS(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
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
!       Water mass - aqueous pressure
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
!
!---      Limit changes in pressure  ---
!
          DPX = 2.5D-1*MAX(PL(2,N)+PATM,PG(2,N)+PATM,1.D+6)
          DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
          DO IGC = 1, NGC
            DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
          ENDDO
          PL(2,N) = PL(2,N) + DPL
          PG(2,N) = PG(2,N) + DPC(1)
          DO IGC = 2,NGC
            PVC(IGC,2,N) = MAX( PVC(IGC,2,N)+DPC(IGC),EPSL )
          ENDDO
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Limit salt mass fraction changes to 0.25 of the
!           maximum value if salt mass fraction is less than
!           the maximum   ---
!
            CALL SOL_LS( T(2,N),XLSMX )
            IF( YLS(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
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
!---      Maintain the gas pressure above or at the water vapor
!         pressure  ---
!
          CALL SP_B( T(2,N),XLS(2,N),PSBX )
          CALL DENS_B( T(2,N),PSBX,XLS(2,N),RHOBX )
          PCX = MAX( PSBX-PL(2,N),0.D+0 )
          IF( ISLC(44).EQ.1 ) THEN
            PVBX = PSBX
          ELSE
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(2,N) )
          ENDIF
          PG(2,N) = MAX( PG(2,N),(PVBX-PATM) )
!
!---    Saturated system w/ entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - trapped gas saturation
!       NaCl mass - total NaCl aqueous mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
!
!---      Limit changes in pressure  ---
!
          DPX = 2.5D-1*MAX(PL(2,N)+PATM,1.D+6)
          DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
          PL(2,N) = PL(2,N) + DPL
!
!---      Relax disappearance of trapped gas  ---
!
          IF( (SG(2,N)+DPC(1)).LT.0.D+0 ) THEN
            SG(2,N) = MAX( (SG(2,N)+6.D-1*DPC(1)),0.D+0 )
          ELSE
            SG(2,N) = MAX( (SG(2,N)+DPC(1)),0.D+0 )
          ENDIF
          IF( SG(2,N).LT.EPSL ) SG(2,N) = 0.D+0
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Limit salt mass fraction changes to 0.25 of the
!           maximum value if salt mass fraction is less than
!           the maximum   ---
!
            CALL SOL_LS( T(2,N),XLSMX )
            IF( YLS(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS(2,N)/EPSL.LT.EPSL .AND. BLU(MPS)/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS(2,N) = YLS(2,N)+DPS
            IF( YLS(2,N).LT.EPSL ) YLS(2,N) = 0.D+0
            XLS(2,N) = MIN( YLS(2,N),XLSMX )
          ENDIF
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
          IF( NPHAZ(2,N).EQ.1 ) THEN
            DO IGC = 1,NGC
              NCHX = INDEX( GCNM(IPNVGC)(1:),'  ' ) - 1
              WRITE(ISC,'(4X,3A,1PE12.5)')'Dissolved ',
     &          GCNM(IPNVGC)(1:NCHX),' Mass Fraction = ',XLC(IGC,2,N)
              WRITE(IWR,'(4X,3A,1PE12.5)')'Dissolved ',
     &          GCNM(IPNVGC)(1:NCHX),' Mass Fraction = ',XLC(IGC,2,N)
            END DO
          ELSE
            WRITE(ISC,'(4X,A,1PE12.5,A,I6)')
     &        'Gas Pressure = ',PG(2,N)+PATM
            WRITE(IWR,'(4X,A,1PE12.5,A,I6)')
     &        'Gas Pressure = ',PG(2,N)+PATM
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
            XLA(2,N) = XLA(1,N)
            SG(2,N) = SG(1,N)
            YLS(2,N) = YLS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
            DO IGC = 1, NGC
             XLC(IGC,2,N) = XLC(IGC,1,N)
             PVC(IGC,2,N) = PVC(IGC,1,N)
            ENDDO
  400     CONTINUE
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
          DO 410 N = 1,NFLD
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            XLA(2,N) = XLA(1,N)
            SG(2,N) = SG(1,N)
            YLS(2,N) = YLS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
            DO IGC = 1, NGC
             XLC(IGC,2,N) = XLC(IGC,1,N)
             PVC(IGC,2,N) = PVC(IGC,1,N)
            ENDDO
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
!---  End of UPDT_COMP group
!
        RETURN
       END

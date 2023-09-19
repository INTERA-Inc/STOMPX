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
!     STOMP-CO2
!
!     This engineering program numerically simulates the transport
!     of H2O, NaCl and CO2 through multifluid subsurface environments
!     under isothermal conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      USE UCODE
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE OUTPU
      USE JACOB
      USE GEO_MECH
      USE FDVP
      USE FDVS
      USE FILES
      USE COUP_WELL
#ifdef valuerss
      USE IFPORT
#endif
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
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
      INTEGER, EXTERNAL :: VALUERSS
!
!----------------------Executable Lines--------------------------------!
!
#if defined(intel) || defined(gnu)
!
!---  Initialize system clock ---
!
      CALL SYSTEM_CLOCK( ICOUNT1,ICRATE,ICMAX )
#endif
#ifdef petsc
!
!---  Initialize PETSc ---
!
      CALL PetscInitialize( PETSC_NULL_CHARACTER,IERR )
#endif
#ifdef lis
!
!---  Initialize Lis ---
!
      CALL lis_initialize( IERR )
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
      IF( ISLC(67).EQ.1 .AND. LFD.GT.12000000 ) THEN
        INDX = 3
        CHMSG = 'Number of Grid Cells Yields Binary Files > 2.4 GB'
        CALL WRMSGS( INDX )
      ENDIF
      ISUB_LOG = 1
      SUB_LOG(1) = 'STOMP-CO2'
      ICODE = 32
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
      CALL RDINPT_CO2
!
!---  Read process ID and check for status file  ---
!
      IF( ISLC(79).EQ.1 ) THEN
#ifdef valuerss
        IPID = GETPID()
#else
        WRITE(6,'(/,A,/)')'Enter the Process ID: '
        READ(5,*) IPID
#endif     
        WRITE(PID_CHAR,'(I8)') IPID
        PIDFNM='/proc/'//TRIM(ADJUSTL(PID_CHAR))//'/status'
        INQUIRE (FILE=PIDFNM,EXIST=IFXST)
        IF (.NOT.IFXST) THEN
          PRINT *,'System File Does Not Exist: Memory Checking Off: ',
     &      PIDFNM
          ISLC(79) = 0
        ENDIF
      ENDIF
!
!---  Create a node connection map  ---
!
      CALL CONNMAP
!
!---  Define leaky well nodes, check leaky well trajectory,
!     set geometric parameters for leaky well nodes  ---
!
      IF( L_LW.EQ.1 ) THEN
        CALL CHK_LEAK_WELL
      ENDIF
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
      CALL CHK_CO2
!
!---  Write bin input files for parallel processing  ---
!
      IF( ISLC(67).EQ.1 ) THEN
        CALL INCRM_CO2
        CALL PROP_CO2
        IF( L_CW.EQ.1 ) THEN
          CALL CHK_COUP_WELL        
          CALL WR_WELL
        ENDIF
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
        CALL WRITE_BIN_CO2
        WRITE(IWR,'(A)') 'NOTE: Preprocessing for STOMPX-CO2'
        WRITE(ISC,'(A)') 'NOTE: Preprocessing for STOMPX-CO2'
!
!---    End parallel pre-processing unless output of the flow or
!       transport linear system is requested  ---
!
        IF( ISLC(34).EQ.0 ) THEN
          WRITE(IWR,'(/,A,/)') '---  End of STOMPX-CO2 Preprocessing' //
     &      '  ---'
          WRITE(ISC,'(/,A,/)') '---  End of STOMPX-CO2 Preprocessing' //
     &      '  ---'
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
!---  Compute primary variable increments  ---
!
      CALL INCRM_CO2
      IF( ICNV.EQ.4 ) GOTO 900
!
!---  Initial hydrologic and thermodynamic properties for active 
!     field nodes ---
!
      CALL PROP_CO2
!
!---  Initial hydrologic and thermodynamic properties for
!     boundary surfaces  ---
!
      CALL BCP_CO2
!
!---  Define coupled-well nodes, check coupled-well trajectory,
!     initialize coupled-well pressure, and increment coupled-well
!     primary variables, skip for parallel preprocessing executions  ---
!
      IF( L_CW.EQ.1 ) THEN
        IF( ISLC(67).EQ.0 ) CALL CHK_COUP_WELL
        CALL INCRM_COUP_WELL
      ENDIF
!
!---  PEST output for CO2 Mineralization at Wallula, store initial
!     pressure of coupled-well #3 and time offset  ---
!
      IF( ISLC(80).EQ.2 ) THEN
        R_OBDT(1,1) = TM
        R_OBDT(4,1) = P_CW(2,3)
        R_OBDT(5,1) = P_CW(2,3)-R_OBDT(4,1)
        R_OBDT(6,1) = P_CW(2,3)-R_OBDT(4,1)
!
!---    Open observational data file for PEST  ---
!
        OPEN(UNIT=46,FILE='wallula_sim.dat',STATUS='UNKNOWN',
     &    FORM='FORMATTED')
        CLOSE(UNIT=46,STATUS='DELETE')
        OPEN(UNIT=46, FILE='wallula_sim.dat', STATUS='NEW',
     &    FORM='FORMATTED')
        WRITE(46,'(A)') 'Pressure Differential for Well #3, MPa'
      ENDIF
!
!---  Write coupled-well and leaky well data to well.dat,
!     skip for parallel preprocessing executions  ---
!
      IF( L_CW.EQ.1 .OR. L_LW.EQ.1 ) THEN
        IF( ISLC(67).EQ.0 ) CALL WR_WELL
      ENDIF
#ifdef ecke
!
!---  Sequence reaction equations  ---
!
      IF( ISLC(40).EQ.1 .AND. ISLC(67).EQ.0 ) CALL SEQEQ
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
!---  Compute initial solute concentrations  ---
!
      CALL CISC_CO2
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
!---  Create PETSc matrix, solver, and solution and problem vectors  ---
!
      IF( ILES.EQ.5 ) THEN
!
!---    PETSc solver for coupled flow and heat transport  ---
!
        INDX = 0
        CALL STOMP_PETSC_CREATE( F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC,INDX )
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
          CALL STOMP_PETSC_CREATE( T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC,
     &      INDX )
        ENDIF
!
!---    Solver for geomechanics  ---
!
        IF( ISLC(50).NE.0 ) THEN
          INDX = 2
          CALL STOMP_PETSC_CREATE( G_KSP,G_MAT,G_RHS_VEC,G_SOL_VEC,
     &      INDX )
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
      CALL DFFGW
      CALL DFFLA
!
!---  Isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS
      CALL BCF_CO2
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
#ifdef ecke
!
!---  Normalize mineral species concentrations after initial
!     output for normal simulations  ---
!
      IF( (NSTEP-NRST).EQ.0 ) CALL NMNSP
#endif
!
!---  Update porosity and permeability in response to geomechanical
!     stress  ---
!
      IF( ISLC(50).NE.0 .AND. NSTEP.EQ.0 ) THEN
        CALL PORSTY_GM
        CALL PERMRF_GM
      ENDIF
!
!---  Load old time step arrays  ---
!
      CALL LDO_CO2
!
!---  Load old time step arrays for the coupled-well model  ---
!
      IF( L_CW.EQ.1 ) THEN
        CALL LDO_COUP_WELL
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
        CALL BCP_CO2
        GOTO 600
      ENDIF
!
!---  Reset the time step reduction counter  ---
!
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
      CALL BCP_CO2
!
!---  Compute coupled-well fluxes  ---
!
      IF( L_CW.EQ.1 ) THEN
!
!---    Enhanced coupled-well model  ---
!
        IF( ISLC(70).EQ.1 ) THEN
          DO 310 NCW = 1,N_CW
!
!---        Injection well  ---
!
            IF( IT_CW(NCW).GT.0 ) THEN
              CALL INJ_COUP_WELL( NCW )
            ENDIF
  310     CONTINUE
!
!---    Standard coupled-well model  ---
!
        ELSE
          CALL FLUX_COUP_WELL
        ENDIF
      ENDIF
!      PRINT *,'P_CW(1,1) = ',P_CW(1,1),' P_CW(2,1) = ',P_CW(2,1)
!
!---  Compute leaky well fluxes  ---
!
      IF( L_LW.EQ.1 ) THEN
        CALL FLUX_LEAK_WELL
      ENDIF
!
!---  Compute source contributions  ---
!
      CALL SORC_CO2
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
      CALL DFFGW
!
!---  Compute aqueous-CO2 diffusion flux through the aqueous phase
!     (interior surfaces)  ---
!
      CALL DFFLA
!
!---  Compute aqueous-NaCl diffusion flux through the aqueous phase
!     (interior surfaces), isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS
!
!---  Compute aqueous-phase volumetric flux, gas-phase volumetric flux,
!     water vapor mass flux, aqueous-phase salt flux
!     (boundary surfaces)  ---
!
      CALL BCF_CO2
!
!---  Zero Jacobian matrix  ---
!
      INDX = 0
      CALL JCBZ( ISVC,MUC,MLC,MKC,INDX )
!
!---  Load Jacobian matrix for the water equation for field nodes
!     (zero flux boundary)  ---
!
      CALL JCBW33
!
!---  Load Jacobian matrix for the CO2 equation for field nodes
!     (zero flux boundary)  ---
!
      CALL JCBA33
!
!---  Load Jacobian matrix for the salt equation for field nodes
!     (zero flux boundary)  ---
!
      IF( ISLC(32).EQ.0 ) CALL JCBS33
!
!---  Modify the Jacobian matrix for boundary conditions  ---
!
      CALL BCJ_CO2
!
!---  Modify Jacobian matrix for coupled-well equations  ---
!
      IF( L_CW.EQ.1 ) THEN
        CALL JCB_COUP_WELL
      ENDIF
!
!---  Modify Jacobian matrix for leaky well equations  ---
!
      IF( L_LW.EQ.1 ) THEN
        CALL JCB_LEAK_WELL
      ENDIF
      IF( ISLC(79).EQ.1 ) THEN
        INDX = VALUERSS()
        PRINT *,'Pre-Solve: ',INDX
      ENDIF
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
        CALL STOMP_LIS_SOLVE( ISVC,F_KSP,F_MAT,F_RHS_VEC,
     &    F_SOL_VEC,INDX )
#endif
#ifdef petsc
      ELSEIF( ILES.EQ.5 ) THEN
        INDX = 0
        CALL STOMP_PETSC_SOLVE( F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC,INDX )
#endif
      ENDIF
      IF( ISLC(79).EQ.1 ) THEN
        INDX = VALUERSS()
        PRINT *,'Post-Solve: ',INDX
      ENDIF
!
!---  Update primary variables for coupled wells  ---
!
      IF( L_CW.EQ.1 ) THEN
        CALL UPDT_COUP_WELL
      ENDIF
!
!---  Update primary variables for nodes  ---
!
      CALL UPDT_CO2
!
!---  Convergence check for couped wells  ---
!
      IF( L_CW.EQ.1 ) THEN
        CALL RSDL_COUP_WELL
      ENDIF
!
!---  Convergence check for nodes  ---
!
      CALL RSDL_CO2
!
!---  Compute primary variable increments, saturation,
!     relative permeability, porosity, tortuosity,
!     thermodynamic properties for interior nodes,
!     except immediately after a new time step  ---
!
      CALL INCRM_CO2
      IF( ICNV.EQ.4 ) GOTO 900
      CALL PROP_CO2
!
!---  For geomechanics simulations alter permeability with
!     porosity  --
!
      IF( ISLC(50).NE.0 ) CALL PERMRF_GM
!
!---  Increment coupled-well primary variables  ---
!
      IF( L_CW.EQ.1 ) THEN
        CALL INCRM_COUP_WELL
      ENDIF
!      PRINT *,'NSTEP = ',NSTEP,' NITER = ',NITER
!      DO M = 1,ISVC
!        PRINT *,'RSD(',M,') = ',RSD(M),
!     &    'NSD(',M,') = ',NSD(M),' ICNV = ',ICNV
!      ENDDO
!      PRINT *,'RSD_CW = ',RSD_CW,' ID_CW(8,1) = ',ID_CW(8,1),
!     &  'P_CW(2,1) = ',P_CW(2,1),'P_CW(3,1) = ',P_CW(3,1)
      GOTO( 200,300,600,900 ) ICNV
  600 CONTINUE
!
!---  PEST output for CO2 Mineralization at Wallula, write pressure
!     differential at well #3 at observation times  ---
!
      IF( ISLC(80).EQ.2 ) CALL WROBDA_WALLUA
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
!        PRINT *,'RSD_GM = ',RSD_GM,'RSDM_GM(IEPD) = ',RSDM_GM(IEPD)
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
!
!---    Update porosity and permeability for geomechical stress  ---
!
        CALL PORSTY_GM
        CALL PERMRF_GM
      ENDIF
!
!---  Integrate coupled-equation sources ---
!
      CALL SORIC_CO2
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
      CALL DFFGW
!
!---  Compute aqueous CO2 diffusion flux through the aqueous phase
!     (interior surfaces)  ---
!
      CALL DFFLA
!
!---  Compute aqueous-NaCl diffusion flux through the aqueous phase
!     (interior surfaces), isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS
!
!---  Compute aqueous-phase volumetric flux, gas-phase volumetric flux,
!     water vapor diffusive, aqueous-CO2 diffusive, and
!     aqueous-NaCl diffusive fluxes (boundary surfaces)  ---
!
      CALL BCF_CO2
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
          CALL SPRP_CO2( NSL )
!
!---      Solute transport ---
!
          CALL TPORT_CO2( NSL )
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
!---        Skip transport for linked aqueous CO2   ---
!
            IF( ISPLK(6).EQ.NSL ) GOTO 720
!
!---        Mobile conservation component fractions   ---
!
            CALL MOBCF( NEQ )
!
!---        Solute transport ---
!
            CALL TPORT_CO2( NSL )
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
            IF( ISPLK(6).EQ.NSL ) GOTO 740
! 
!---        Mobile kinetic component fractions   ---
!
            CALL MOBKF( NEQ )
!
!---        Solute transport ---
!
            CALL TPORT_CO2( NSL )
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
!
!---  PEST output for CO2 Mineralization at Wallula, close
!     observation data file  ---
!
      IF( ISLC(80).EQ.2 ) CLOSE(UNIT=46)
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
!---  Destroy solution vector, problem vector, matrix, solver, and
!     solution-vector scatter scheme for the PETSc solver for coupled
!     flow and heat transport  ---
!
      CALL STOMP_PETSC_DESTROY( F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC )
!
!---  Destroy solution vector, problem vector, matrix, solver, and
!     solution-vector scatter scheme for the PETSc solver for 
!     solute transport  ---
!
      IF( IEQC.NE.0 .OR. ISLC(40).NE.0 ) THEN
        CALL STOMP_PETSC_DESTROY( T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC )
      ENDIF
!
!---  Destroy solution vector, problem vector, matrix, solver, and
!     solution-vector scatter scheme for the PETSc solver for 
!     geomechanics  ---
!
      IF( ISLC(50).NE.0 ) THEN
        CALL STOMP_PETSC_DESTROY( G_KSP,G_MAT,G_RHS_VEC,G_SOL_VEC )
      ENDIF
!
!---  Finalize PETSc execution  ---
!
      CALL PetscFinalize( IERR )
#endif
      STOP
!
!---  End of STOMP program  ---
!
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCF_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Compute boundary surface fluxes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
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
      SUB_LOG(ISUB_LOG) = '/BCF_CO2'
!
!---  Zero boundary fluxes  ---
!
      DO 70 NB = 1,NBC
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
          DO 10 M = 1,ISVF
            WL(M,NPZ) = 0.D+0
            WG(M,NPZ) = 0.D+0
            WDGW(M,NPZ) = 0.D+0
            WDLA(M,NPZ) = 0.D+0
            WS(M,NPZ) = 0.D+0
            WDS(M,NPZ) = 0.D+0
   10     CONTINUE
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          DO 20 M = 1,ISVF
            VL(M,NPY) = 0.D+0
            VG(M,NPY) = 0.D+0
            VDGW(M,NPY) = 0.D+0
            VDLA(M,NPY) = 0.D+0
            VS(M,NPY) = 0.D+0
            VDS(M,NPY) = 0.D+0
   20     CONTINUE
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          DO 30 M = 1,ISVF
            UL(M,NPX) = 0.D+0
            UG(M,NPX) = 0.D+0
            UDGW(M,NPX) = 0.D+0
            UDLA(M,NPX) = 0.D+0
            US(M,NPX) = 0.D+0
            UDS(M,NPX) = 0.D+0
   30     CONTINUE
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          DO 40 M = 1,ISVF
            UL(M,NQX) = 0.D+0
            UG(M,NQX) = 0.D+0
            UDGW(M,NQX) = 0.D+0
            UDLA(M,NQX) = 0.D+0
            US(M,NQX) = 0.D+0
            UDS(M,NQX) = 0.D+0
   40     CONTINUE
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          DO 50 M = 1,ISVF
            VL(M,NQY) = 0.D+0
            VG(M,NQY) = 0.D+0
            VDGW(M,NQY) = 0.D+0
            VDLA(M,NQY) = 0.D+0
            VS(M,NQY) = 0.D+0
            VDS(M,NQY) = 0.D+0
   50     CONTINUE
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          DO 60 M = 1,ISVF
            WL(M,NQZ) = 0.D+0
            WG(M,NQZ) = 0.D+0
            WDGW(M,NQZ) = 0.D+0
            WDLA(M,NQZ) = 0.D+0
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
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 110 M = 1,ISVF
              WL(M,NPZ) = BCX(3)
  110       CONTINUE
            CALL DFFLAB( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLB( N,NB )
            CALL DFFLAB( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO 111 M = 1,ISVF
              WG(M,NPZ) = BCX(4)
  111       CONTINUE
            CALL DFFGWB( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGB( N,NB )
            CALL DFFGWB( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSB( N,NB )
            ENDIF
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 120 M = 1,ISVF
              VL(M,NPY) = BCX(3)
  120       CONTINUE
            CALL DFFLAS( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLS( N,NB )
            CALL DFFLAS( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO 121 M = 1,ISVF
              VG(M,NPY) = BCX(4)
  121       CONTINUE
            CALL DFFGWS( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGS( N,NB )
            CALL DFFGWS( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSS( N,NB )
            ENDIF
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 130 M = 1,ISVF
              UL(M,NPX) = BCX(3)
  130       CONTINUE
            CALL DFFLAW( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLW( N,NB )
            CALL DFFLAW( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO 131 M = 1,ISVF
              UG(M,NPX) = BCX(4)
  131       CONTINUE
            CALL DFFGWW( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGW( N,NB )
            CALL DFFGWW( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSW( N,NB )
            ENDIF
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 140 M = 1,ISVF
              UL(M,NQX) = BCX(3)
  140       CONTINUE
            CALL DFFLAE( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLE( N,NB )
            CALL DFFLAE( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO 141 M = 1,ISVF
              UG(M,NQX) = BCX(4)
  141       CONTINUE
            CALL DFFGWE( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGE( N,NB )
            CALL DFFGWE( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSE( N,NB )
            ENDIF
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 150 M = 1,ISVF
              VL(M,NQY) = BCX(3)
  150       CONTINUE
            CALL DFFLAN( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLN( N,NB )
            CALL DFFLAN( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO 151 M = 1,ISVF
              VG(M,NQY) = BCX(4)
  151       CONTINUE
            CALL DFFGWN( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGN( N,NB )
            CALL DFFGWN( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSN( N,NB )
            ENDIF
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO 160 M = 1,ISVF
              WL(M,NQZ) = BCX(3)
  160       CONTINUE
            CALL DFFLAT( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLT( N,NB )
            CALL DFFLAT( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO 161 M = 1,ISVF
              WG(M,NQZ) = BCX(4)
  161       CONTINUE
            CALL DFFGWT( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGT( N,NB )
            CALL DFFGWT( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt Boundary Conditions  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLST( N,NB )
            ENDIF
          ENDIF
        ENDIF
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCF_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCJ_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Modify the Jacobian matrix for boundary conditions
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      SUB_LOG(ISUB_LOG) = '/BCJ_CO2'
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
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWB( N,NB,NPZ )
            CALL JCBLAB( N,NB,NPZ )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWB( N,NB,NPZ )
            CALL JCBGAB( N,NB,NPZ )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSB( N,NB,NPZ )
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWS( N,NB,NPY )
            CALL JCBLAS( N,NB,NPY )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWS( N,NB,NPY )
            CALL JCBGAS( N,NB,NPY )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSS( N,NB,NPY )
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWW( N,NB,NPX )
            CALL JCBLAW( N,NB,NPX )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWW( N,NB,NPX )
            CALL JCBGAW( N,NB,NPX )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSW( N,NB,NPX )
          ENDIF
!
!---    East boundary
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWE( N,NB,NQX )
            CALL JCBLAE( N,NB,NQX )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWE( N,NB,NQX )
            CALL JCBGAE( N,NB,NQX )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSE( N,NB,NQX )
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWN( N,NB,NQY )
            CALL JCBLAN( N,NB,NQY )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWN( N,NB,NQY )
            CALL JCBGAN( N,NB,NQY )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSN( N,NB,NQY )
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWT( N,NB,NQZ )
            CALL JCBLAT( N,NB,NQZ )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWT( N,NB,NQZ )
            CALL JCBGAT( N,NB,NQZ )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt  ---
!
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBST( N,NB,NQZ )
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCJ_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCP_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Compute saturation, relative permeability and thermodynamic
!     properties for boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      USE FDVG
      USE CONST
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
      REAL*8 BCX(LBCV),YLSX(3)
      REAL*8 GX(3),RX(2),RPX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCP_CO2'
!
!---  Assign values for initial condition type boundary conditions  ---
!
      IF( NSTEP-NRST.LE.1 .AND. NITER.LE.1 ) THEN
        DO 50 NB = 1,NBC
          N = IBCN(NB)
          IF( IBCT(IEQW,NB).EQ.12 .OR. IBCT(IEQA,NB).EQ.12 ) THEN
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
          IF( IBCT(IEQW,NB).EQ.12 ) THEN
            PLB(1,NB) = PL(2,N) + RHOL(2,N)*GB
            PLB(2,NB) = PLB(1,NB)
          ENDIF
          IF( IBCT(IEQA,NB).EQ.12 ) THEN
            IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.3 .OR.
     &        NPHAZ(2,N).EQ.5 .OR. NPHAZ(2,N).EQ.7 ) THEN
              PGB(1,NB) = PG(2,N) + RHOL(2,N)*GB
            ELSE
              PGB(1,NB) = PG(2,N) + RHOG(2,N)*GB
            ENDIF
            PGB(2,NB) = PGB(1,NB)
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 )  THEN
            IF( IBCT(IEQS,NB).EQ.12 ) YLSB(1,NB) = YLS(2,N)
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
             IF( IBCT(IEQW,NB).EQ.2 ) THEN
               BCX(3) = BCX(3)-5.D-1*DTBC*(BC(3,M,MB)-BC(3,M-1,MB))/TDBC
             ELSEIF( IBCT(IEQW,NB).EQ.12 ) THEN
               BCX(3) = PLB(1,NB)
             ENDIF
             IF( IBCT(IEQA,NB).EQ.2 ) THEN
               BCX(4) = BCX(4)-5.D-1*DTBC*(BC(4,M,MB)-BC(4,M-1,MB))/TDBC
             ELSEIF( IBCT(IEQA,NB).EQ.12 ) THEN
               BCX(4) = PGB(1,NB)
             ENDIF
!
!---         Isobrine option  ---
!
             IF( ISLC(32).EQ.0 )  THEN
              IF( IBCT(IEQS,NB).EQ.2 ) THEN
               BCX(6) = BCX(6)-5.D-1*DTBC*(BC(6,M,MB)-BC(6,M-1,MB))/TDBC
              ELSEIF( IBCT(IEQS,NB).EQ.12 ) THEN
                BCX(6) = YLSB(1,NB)
              ENDIF
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
        IF( IBCT(IEQA,NB).EQ.12 ) BCX(4) = PGB(1,NB)
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 )  THEN
          IF( IBCT(IEQS,NB).EQ.12 ) BCX(6) = YLSB(1,NB)
        ENDIF
        N = IBCN(NB)
        N_DB = -NB
        IBD = ABS(IBCD(NB))
        IZN = IZ(N)
!#ifdef ecke
!        POR0(1,N) = POR0(1,N)
!        POR0(2,N) = POR0(2,N)
!#endif
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
           DB = 0.5D+0*DZGF(N)
           DB = DZGP(NPZ)
           GB = GRVZ(NPZ)*DB
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
           DB = 0.5D+0*DYGF(N)*RP(ID(N))
           DB = DYGP(NPY)*RP(ID(N))
           GB = GRVY(NPY)*DB
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
           DB = 0.5D+0*DXGF(N)
           DB = DXGP(NPX)
           GB = GRVX(NPX)*DB
        ELSEIF( IBCD(NB).EQ.1 ) THEN
           DB = -0.5D+0*DXGF(N)
           DB = -DXGP(NQX)
           GB = GRVX(NQX)*DB
        ELSEIF( IBCD(NB).EQ.2 ) THEN
           DB = -0.5D+0*DYGF(N)*RP(ID(N))
           DB = -DYGP(NQY)*RP(ID(N))
           GB = GRVY(NQY)*DB
        ELSEIF( IBCD(NB).EQ.3 ) THEN
           DB = -0.5D+0*DZGF(N)
           DB = -DZGP(NQZ)
           GB = GRVZ(NQZ)*DB
        ENDIF
!
!---    Loop over secondary variable indices  ---
!
        DO 300 M = 2,ISVC+2
          TX = T(2,N)
          PLX = PL(M,N)
          PGX = PG(M,N)
!
!---      Gas Dirichlet  ---
!
          IF( IBCT(IEQA,NB).EQ.1 .OR. IBCT(IEQA,NB).EQ.12 ) THEN
            PGX = BCX(4)
!
!---      Gas Neumann  ---
!
          ELSEIF( IBCT(IEQA,NB).EQ.2 ) THEN
            PGX = PGX + BCX(4)*DB*VISG(M,N)/PERM(IBD,IZN)
     &        + RHOG(M,N)*GB
!
!---      Gas Zero Flux  ---
!
          ELSEIF( IBCT(IEQA,NB).EQ.3 ) THEN
            IF( ABS(BCX(4)+PATM).GT.EPSL ) THEN
              PGX = BCX(4)
            ELSE
              PGX = PGX + RHOG(M,N)*GB
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
!---      Aqueous Hydrostatic  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.11 ) THEN
            IF( IBCM(NB).EQ.1 .AND. (NSTEP-NRST).GT.1 ) GOTO 400
            IF( M.EQ.2 ) THEN
              N_DB = -NB
              CALL HYDST_BC_CO2( BCX,PGX,PLX,TX,XLAB(M,NB),XLSB(M,NB),
     &          YLSB(M,NB),ZP(N) )
            ELSE
              TX = TB(2,NB)
              PLX = PLB(2,NB)
              PGX = PGB(2,NB)
              XLAB(M,NB) = XLAB(2,NB)
              XLSB(M,NB) = XLSB(2,NB)
              YLSB(M,NB) = YLSB(2,NB)
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
!---      Hydrostatic boundary conditions  ---
!
          IF( IBCT(IEQW,NB).EQ.11 ) THEN
            CALL SP_B( TX,YLSB(M,NB),PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( TX,PX,YLSB(M,NB),RHOBX )
            CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),IZN )
            PVWB(M,NB) = PVBX
            PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
            XLWB(M,NB) = MAX( 1.D+0-XLAB(M,NB)-XLSB(M,NB),0.D+0 )
            WTMLX = 1.D+0/(XLAB(M,NB)/WTMA + XLSB(M,NB)/WTMS + 
     &        XLWB(M,NB)/WTMW)
            XMLAB(M,NB) = XLAB(M,NB)*WTMLX/WTMA
            XMLSB(M,NB) = XLSB(M,NB)*WTMLX/WTMS
            XMLWB(M,NB) = XLWB(M,NB)*WTMLX/WTMW
            XGAB(M,NB) = 0.D+0
            XGWB(M,NB) = 1.D+0
            XMGAB(M,NB) = 0.D+0
            XMGWB(M,NB) = 1.D+0
          ELSE
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Aqueous-salt concentration boundary condition  ---
!
            IF( IBCT(IEQS,NB).EQ.40 .OR. IBCT(IEQS,NB).EQ.41 ) THEN
              RHOLSX = BCX(6)
              ISRX = 1
              CALL DENS_W( TX,PX,RHOLWX,RHOGWX,ISRX )
!
!---          Guess brine salt mass fraction  ---
!
              YLSX(2) = RHOLSX/(RHOLWX+RHOLSX)
              CALL SOL_LS( TX,XLSMX )
              IF( YLSX(2).GE.XLSMX ) THEN
                YLSB(M,NB) = XLSMX
                GOTO 204
              ENDIF
              DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSX(2)))
              NC = 0
!
!---          Single-variable Newton-Raphson loop (YLS)  ---
!
  200         CONTINUE
              NC = NC + 1
              IF( NC.GT.32 ) THEN
                INDX = 17
                RLMSG = BCX(6)
                N_DB = -NB
                CHMSG = 'Unconverged Boundary Condition: '
     &            // 'Salt Aqu. Concentration = '
                CALL WRMSGS( INDX )
              ENDIF
              DO 202 L = 2,3
                YLSX(L) = YLSX(2)
                IF( L.EQ.3 ) YLSX(L) = YLSX(2) + DYLSX
                CALL SP_B( TX,YLSX(L),PSBX )
                PX = MAX( PGX,PSBX )
                CALL DENS_B( TX,PX,YLSX(L),RHOBX )
                CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSX(L),IZN )
                PVWB(M,NB) = BCX(5)*PVBX
                CALL SOL_LS( TX,XLSMX )
                XLSSX = MIN( YLSX(L),XLSMX )
                CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAX,XGWX,
     &            XLASX,XLSSX,XLWSX,XMGAX,XMGWX,XMLAX,XMLSX,XMLWX )
                IF( BCX(7).LT.0 ) THEN
                  XLAX = MIN(-BCX(7),XLASX) ! mass fraction
                ELSE
                  XLAX = BCX(7)*XLASX       ! relative saturation
                ENDIF
                XLSX = YLSX(L) + (XLSSX-YLSX(L))*(XLAX/XLASX)
                CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
                GX(L) = XLSX - RHOLSX/RHOLX
  202         CONTINUE
              RX(1) = (GX(3)-GX(2))/DYLSX
              RPX(1) = -GX(2)
              CYLSX = RPX(1)/RX(1)
              YLSX(2) = MAX( YLSX(2)+CYLSX,0.D+0 )
              IF( YLSX(2).GE.XLSMX ) THEN
                YLSB(M,NB) = XLSMX
                GOTO 204
              ENDIF
              IF( ABS(CYLSX).GT.(1.D-9*XLSMX) ) GOTO 200
              YLSB(M,NB) = YLSX(2)
  204         CONTINUE
              CALL SP_B( TX,YLSB(M,NB),PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( TX,PX,YLSB(M,NB),RHOBX )
              CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),IZN )
              PVWB(M,NB) = BCX(5)*PVBX
              PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
              CALL SOL_LS( TX,XLSMX )
              XLSSX = MIN( YLSB(M,NB),XLSMX )
              CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAB(M,NB),
     &          XGWB(M,NB),XLASX,XLSSX,XLWSX,XMGAB(M,NB),XMGWB(M,NB),
     &          XMLAX,XMLSX,XMLWX )
              IF( BCX(7).LT.0 ) THEN
                XLAB(M,NB) = MIN(-BCX(7),XLASX) ! mass fraction
              ELSE
                XLAB(M,NB) = BCX(7)*XLASX       ! relative saturation
              ENDIF
              XLSB(M,MB) = YLSB(M,NB) + 
     &          (XLSSX-YLSB(M,NB))*(XLAB(M,NB)/XLASX)
              XLWB(M,NB) = MAX( 1.D+0-XLAB(M,NB)-XLSB(M,NB),0.D+0 )
              WTMLX = 1.D+0/(XLAB(M,NB)/WTMA + XLSB(M,NB)/WTMS + 
     &          XLWB(M,NB)/WTMW)
              XMLAB(M,NB) = XLAB(M,NB)*WTMLX/WTMA
              XMLSB(M,NB) = XLSB(M,NB)*WTMLX/WTMS
              XMLWB(M,NB) = XLWB(M,NB)*WTMLX/WTMW
!
!---        Aqueous-salt relative saturation boundary condition  ---
!
            ELSEIF( IBCT(IEQS,NB).EQ.34 .OR. IBCT(IEQS,NB).EQ.35 ) THEN
              PHILSX = BCX(6)
              CALL SOL_LS( TX,XLSMX )
              YLSB(M,NB) = PHILSX*XLSMX
              CALL SP_B( TX,YLSB(M,NB),PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( TX,PX,YLSB(M,NB),RHOBX )
              CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),IZN )
              PVWB(M,NB) = BCX(5)*PVBX
              PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
              CALL SOL_LS( TX,XLSMX )
              XLSSX = MIN( YLSB(M,NB),XLSMX )
              CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAB(M,NB),
     &          XGWB(M,NB),XLASX,XLSSX,XLWSX,XMGAB(M,NB),XMGWB(M,NB),
     &          XMLAX,XMLSX,XMLWX )
              IF( BCX(7).LT.0 ) THEN
                XLAB(M,NB) = MIN(-BCX(7),XLASX) ! mass fraction
              ELSE
                XLAB(M,NB) = BCX(7)*XLASX       ! relative saturation
              ENDIF
              XLSB(M,NB) = BCX(6)
              XLWB(M,NB) = MAX( 1.D+0-XLAB(M,NB)-XLSB(M,NB),0.D+0 )
              WTMLX = 1.D+0/(XLAB(M,NB)/WTMA + XLSB(M,NB)/WTMS + 
     &          XLWB(M,NB)/WTMW)
              XMLAB(M,NB) = XLAB(M,NB)*WTMLX/WTMA
              XMLSB(M,NB) = XLSB(M,NB)*WTMLX/WTMS
              XMLWB(M,NB) = XLWB(M,NB)*WTMLX/WTMW
!
!---        Aqueous-salt mass fraction boundary condition  ---
!
            ELSEIF( IBCT(IEQS,NB).EQ.36 .OR. IBCT(IEQS,NB).EQ.37 .OR.
     &        IBCT(IEQS,NB).EQ.12 ) THEN
!
!---          Guess brine salt mass fraction  ---
!
              YLSX(2) = BCX(6)
              CALL SOL_LS( TX,XLSMX )
              IF( YLSX(2).GE.XLSMX ) THEN
                YLSB(M,NB) = XLSMX
                GOTO 214
              ENDIF
              DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSX(2)))
              NC = 0
!
!---          Single-variable Newton-Raphson loop (YLS)  ---
!
  210         CONTINUE
              NC = NC + 1
              IF( NC.GT.32 ) THEN
                INDX = 17
                RLMSG = BCX(6)
                N_DB = -NB
                CHMSG = 'Unconverged Boundary Condition: '
     &            // 'Salt Aqu. Mass Fraction = '
                CALL WRMSGS( INDX )
              ENDIF
              DO 212 L = 2,3
                YLSX(L) = YLSX(2)
                IF( L.EQ.3 ) YLSX(L) = YLSX(2) + DYLSX
                CALL SP_B( TX,YLSX(L),PSBX )
                PX = MAX( PGX,PSBX )
                CALL DENS_B( TX,PX,YLSX(L),RHOBX )
                CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSX(L),IZN )
                PVWB(M,NB) = BCX(5)*PVBX
                CALL SOL_LS( TX,XLSMX )
                XLSSX = MIN( YLSX(L),XLSMX )
                CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAX,XGWX,
     &            XLASX,XLSSX,XLWSX,XMGAX,XMGWX,XMLAX,XMLSX,XMLWX )
                IF( BCX(7).LT.0 ) THEN
                  XLAX = MIN(-BCX(7),XLASX) ! mass fraction
                ELSE
                  XLAX = BCX(7)*XLASX       ! relative saturation
                ENDIF
                XLSX = YLSX(L) + (XLSSX-YLSX(L))*(XLAX/XLASX)
                GX(L) = BCX(6) - XLSX
  212         CONTINUE
              RX(1) = (GX(3)-GX(2))/DYLSX
              RPX(1) = -GX(2)
              CYLSX = RPX(1)/RX(1)
              YLSX(2) = MAX( YLSX(2)+CYLSX,0.D+0 )
              IF( YLSX(2).GE.XLSMX ) THEN
                YLSB(M,NB) = XLSMX
                GOTO 214
              ENDIF
              IF( ABS(CYLSX).GT.(1.D-9*XLSMX) ) GOTO 210
              YLSB(M,NB) = YLSX(2)
  214         CONTINUE
              CALL SP_B( TX,YLSB(M,NB),PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( TX,PX,YLSB(M,NB),RHOBX )
              CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),IZN )
              PVWB(M,NB) = BCX(5)*PVBX
              PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
              CALL SOL_LS( TX,XLSMX )
              XLSSX = MIN( YLSB(M,NB),XLSMX )
              CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAB(M,NB),
     &          XGWB(M,NB),XLASX,XLSSX,XLWSX,XMGAB(M,NB),XMGWB(M,NB),
     &          XMLAX,XMLSX,XMLWX )
              IF( BCX(7).LT.0 ) THEN
                XLAB(M,NB) = MIN(-BCX(7),XLASX) ! mass fraction
              ELSE
                XLAB(M,NB) = BCX(7)*XLASX       ! relative saturation
              ENDIF
              XLSB(M,NB) = BCX(6)
              XLWB(M,NB) = MAX( 1.D+0-XLAB(M,NB)-XLSB(M,NB),0.D+0 )
              WTMLX = 1.D+0/(XLAB(M,NB)/WTMA + XLSB(M,NB)/WTMS + 
     &          XLWB(M,NB)/WTMW)
              XMLAB(M,NB) = XLAB(M,NB)*WTMLX/WTMA
              XMLSB(M,NB) = XLSB(M,NB)*WTMLX/WTMS
              XMLWB(M,NB) = XLWB(M,NB)*WTMLX/WTMW
!
!---        Aqueous-salt molality boundary condition  ---
!
            ELSEIF( IBCT(IEQS,NB).EQ.38 .OR. IBCT(IEQS,NB).EQ.39 ) THEN
!
!---          Guess brine salt mass fraction  ---
!
              YLSX(2) = BCX(6)*WTMS/(1.D+0 + BCX(6)*WTMS)
              CALL SOL_LS( TX,XLSMX )
              IF( YLSX(2).GE.XLSMX ) THEN
                YLSB(M,NB) = XLSMX
                GOTO 224
              ENDIF
              DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSX(2)))
              NC = 0
!
!---          Single-variable Newton-Raphson loop (YLS)  ---
!
  220         CONTINUE
              NC = NC + 1
              IF( NC.GT.32 ) THEN
                INDX = 17
                RLMSG = BCX(6)
                N_DB = -NB
                CHMSG = 'Unconverged Boundary Condition: '
     &            // 'Salt Aqu. Molality = '
                CALL WRMSGS( INDX )
              ENDIF
              DO 222 L = 2,3
                YLSX(L) = YLSX(2)
                IF( L.EQ.3 ) YLSX(L) = YLSX(2) + DYLSX
                CALL SP_B( TX,YLSX(L),PSBX )
                PX = MAX( PGX,PSBX )
                CALL DENS_B( TX,PX,YLSX(L),RHOBX )
                CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSX(L),IZN )
                PVWB(M,NB) = BCX(5)*PVBX
                CALL SOL_LS( TX,XLSMX )
                XLSSX = MIN( YLSX(L),XLSMX )
                CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAX,XGWX,
     &            XLASX,XLSSX,XLWSX,XMGAX,XMGWX,XMLAX,XMLSX,XMLWX )
                IF( BCX(7).LT.0 ) THEN
                  XLAX = MIN(-BCX(7),XLASX) ! mass fraction
                ELSE
                  XLAX = BCX(7)*XLASX       ! relative saturation
                ENDIF
                XLSX = YLSX(L) + (XLSSX-YLSX(L))*(XLAX/XLASX)
                GX(L) = (BCX(6)*WTMS/(1.D+0 + BCX(6)*WTMS)) - XLSX
  222         CONTINUE
              RX(1) = (GX(3)-GX(2))/DYLSX
              RPX(1) = -GX(2)
              CYLSX = RPX(1)/RX(1)
              YLSX(2) = MAX( YLSX(2)+CYLSX,0.D+0 )
              IF( YLSX(2).GE.XLSMX ) THEN
                YLSB(M,NB) = XLSMX
                GOTO 224
              ENDIF
              IF( ABS(CYLSX).GT.(1.D-9*XLSMX) ) GOTO 220
              YLSB(M,NB) = YLSX(2)
  224         CONTINUE
              CALL SP_B( TX,YLSB(M,NB),PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( TX,PX,YLSB(M,NB),RHOBX )
              CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),IZN )
              PVWB(M,NB) = BCX(5)*PVBX
              PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
              CALL SOL_LS( TX,XLSMX )
              XLSSX = MIN( YLSB(M,NB),XLSMX )
              CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAB(M,NB),
     &          XGWB(M,NB),XLASX,XLSSX,XLWSX,XMGAB(M,NB),XMGWB(M,NB),
     &          XMLAX,XMLSX,XMLWX )
              IF( BCX(7).LT.0 ) THEN
                XLAB(M,NB) = MIN(-BCX(7),XLASX) ! mass fraction
              ELSE
                XLAB(M,NB) = BCX(7)*XLASX       ! relative saturation
              ENDIF
              XLSB(M,NB) = BCX(6)*WTMS/(1.D+0 + BCX(6)*WTMS)
              XLWB(M,NB) = MAX( 1.D+0-XLAB(M,NB)-XLSB(M,NB),0.D+0 )
              WTMLX = 1.D+0/(XLAB(M,NB)/WTMA + XLSB(M,NB)/WTMS + 
     &          XLWB(M,NB)/WTMW)
              XMLAB(M,NB) = XLAB(M,NB)*WTMLX/WTMA
              XMLSB(M,NB) = XLSB(M,NB)*WTMLX/WTMS
              XMLWB(M,NB) = XLWB(M,NB)*WTMLX/WTMW
!
!---        Zero-flux or outflow boundary condition  ---
!
            ELSEIF( IBCT(IEQS,NB).EQ.3 .OR. IBCT(IEQS,NB).EQ.7 ) THEN
              YLSB(M,NB) = YLS(M,N)
              CALL SP_B( TX,YLSB(M,NB),PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( TX,PX,YLSB(M,NB),RHOBX )
              CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),IZN )
              PVWB(M,NB) = BCX(5)*PVBX
              PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
              CALL SOL_LS( TX,XLSMX )
              XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
              CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAB(M,NB),
     &          XGWB(M,NB),XLASX,XLSB(M,NB),XLWSX,XMGAB(M,NB),
     &          XMGWB(M,NB),XMLAX,XMLSX,XMLWX )
              IF( BCX(7).LT.0 ) THEN
                XLAB(M,NB) = MIN(-BCX(7),XLASX) ! mass fraction
              ELSE
                XLAB(M,NB) = BCX(7)*XLASX       ! relative saturation
              ENDIF
              XLWB(M,NB) = MAX( 1.D+0-XLAB(M,NB)-XLSB(M,NB),0.D+0 )
              WTMLX = 1.D+0/(XLAB(M,NB)/WTMA + XLSB(M,NB)/WTMS + 
     &          XLWB(M,NB)/WTMW)
              XMLAB(M,NB) = XLAB(M,NB)*WTMLX/WTMA
              XMLSB(M,NB) = XLSB(M,NB)*WTMLX/WTMS
              XMLWB(M,NB) = XLWB(M,NB)*WTMLX/WTMW
            ENDIF
!
!---      No salt  ---
!
          ELSE
            YLSB(M,NB) = 0.D+0
            CALL SP_B( TX,YLSB(M,NB),PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( TX,PX,YLSB(M,NB),RHOBX )
            CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),IZN )
            PVWB(M,NB) = BCX(5)*PVBX
            PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
            CALL SOL_LS( TX,XLSMX )
            XLSSX = MIN( YLSB(M,NB),XLSMX )
            CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAB(M,NB),
     &        XGWB(M,NB),XLASX,XLSSX,XLWSX,XMGAB(M,NB),XMGWB(M,NB),
     &        XMLAX,XMLSX,XMLWX )
              IF( BCX(7).LT.0 ) THEN
                XLAB(M,NB) = MIN(-BCX(7),XLASX) ! mass fraction
              ELSE
                XLAB(M,NB) = BCX(7)*XLASX       ! relative saturation
              ENDIF
            XLSB(M,NB) = BCX(6)*WTMS/(1.D+0 + BCX(6)*WTMS)
            XLWB(M,NB) = MAX( 1.D+0-XLAB(M,NB)-XLSB(M,NB),0.D+0 )
            WTMLX = 1.D+0/(XLAB(M,NB)/WTMA + XLSB(M,NB)/WTMS + 
     &        XLWB(M,NB)/WTMW)
            XMLAB(M,NB) = XLAB(M,NB)*WTMLX/WTMA
            XMLSB(M,NB) = XLSB(M,NB)*WTMLX/WTMS
            XMLWB(M,NB) = XLWB(M,NB)*WTMLX/WTMW
          ENDIF
          ENDIF
!!
!!---      Gas composition  ---
!!
!          CALL SP_B( TX,XLSX,PSBX )
!          PX = MAX( PGX,PSBX )
!          CALL DENS_B( TX,PX,XLSX,RHOBX )
!          CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
!          PVBX = BCX(5)*PVBX
!          PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
!          IF( PVAB(M,NB).GT.EPSL ) THEN
!            CALL SOL_LS( TX,XLSMX )
!            XLSX = MIN( YLSB(M,NB),XLSMX )
!            XLSB(M,NB) = XLSX
!            CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,XGAX,XGWX,XLAX,
!     &        XLSB(M,NB),XLWX,XMGAX,XMGWX,XMLAX,XMLSX,XMLWX )
!            XGAB(M,NB) = XGAX
!            XGWB(M,NB) = XGWX
!          ELSE
!            XGAB(M,NB) = 0.D+0
!            XGWB(M,NB) = 1.D+0
!          ENDIF
!
!---      Porous-media porosity  ---
!
          CALL PORSTY_CO2E( N,PX,PCMP(N),PORDB(M,NB),PORTB(M,NB) )
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
          CALL KSP_CO2( N,PGX,PLX,SGB(M,NB),SGTX,SLB(M,NB),SLFX,SLMX,
     &      RKLB(1,M,NB),RKGB(M,NB),ASLX,ASLMINX,ESGTX,ESGTMX,
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
            CALL DENS_A( TX,PVAB(M,NB),RHOGAX,I_VX )
            IF( ISLC(9).EQ.1 ) THEN
              RHOGB(M,NB) = RHOGI
            ELSE
              RHOGB(M,NB) = XGAB(M,NB)*RHOGAX + XGWB(M,NB)*RHOGWX
            ENDIF
            WTMGX = 1.D+0/(XGAB(M,NB)/WTMA + XGWB(M,NB)/WTMW)
            RHOMGB(M,NB) = RHOGB(M,NB)/WTMGX
          ELSE
            IF( ISLC(9).EQ.1 ) THEN
              RHOGB(M,NB) = RHOGI
            ELSE
              RHOGB(M,NB) = RHOGWX
            ENDIF
            WTMGX = WTMW
            RHOMGB(M,NB) = RHOGB(M,NB)/WTMGX
          ENDIF
!
!---      Gas viscosity  ---
!
          CALL VISC_A( TX,RHOGAX,VISGAX )
          CALL VISC_W( TX,PVBX,RHOGWX,VISGWX )
          CALL VISC_G( VISGAX,VISGWX,XMGAB(M,NB),XMGWB(M,NB),
     &      VISGB(M,NB) )
!
!---      Water-vapor diffusion coefficient  ---
!
          CALL DIFC_GW( TX,PX,DFGWB(M,NB) )
!
!---      Aqueous component fractions and density  ---
!
          CALL DENS_L( TX,RHOBX,XLAB(M,NB),RHOLB(M,NB) )
          WTMLX = 1.D+0/(XLAB(M,NB)/WTMA + XLSB(M,NB)/WTMS + 
     &      XLWB(M,NB)/WTMW)
          RHOMLB(M,NB) = RHOLB(M,NB)/WTMLX
!
!---      Aqueous viscosity  ---
!
          ISRX = 1
          CALL DENS_W( TX,PX,RHOLWX,RHOX,ISRX )
          CALL VISC_W( TX,PX,RHOLWX,VISLWX )
          CALL VISC_B( TX,XLSB(M,NB),VISLWX,VISBX )
          CALL VISC_L( XMLAB(M,NB),VISBX,VISGAX,VISLB(M,NB) )
!
!---      Aqueous-CO2 and -NaCl diffusion coefficient  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLAB(M,NB) = DFLAC
            DFLSB(M,NB) = DFLSC
          ELSEIF( ISLC(4).EQ.2 ) THEN
!            CALL DIFC_LA( VISLB(M,NB),VISGAX,DFLAB(M,NB) )
            CALL DIFC_LA( TX,XLSB(M,NB),DFLAB(M,NB) )
            CALL DIFC_LS( TX,XLSB(M,NB),VISLB(M,NB),DFLSB(M,NB) )
          ENDIF
!
!---      Assign boundary primary variables  ---
!
          TB(M,NB) = TX
          PLB(M,NB) = PLX - PATM
          PGB(M,NB) = PGX - PATM
          PVWB(M,NB) = PVBX
!
!---      Gas pressure for Darcy velocity calculation  ---
!
          PSOB(M,NB) = PGB(M,NB)
  300   CONTINUE
  400 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCP_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CAP_CO2( ASLMINX,BTGLX,CPGL,SLX,SGTX,N )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Compute the gas/aqueous capillary pressure from the aqueous
!     saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 GX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CAP_CO2'
      IZN = IZ(N)
      IF( SLX.LT.EPSL ) THEN
        CPGL = SCHR(12,IZN)*RHORL*GRAV/BTGLX
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
     &        SCHRV(1,N)
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
     &        SCHRV(1,N)
          ELSE
            HDGL = SCHR(12,IZN)
            SLX = SLRX + 1.D-6
          ENDIF
        ENDIF
        CPGL = HDGL*RHORL*GRAV/BTGLX
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
            HDGL = SCHRV(1,N)*(1.D+0/ESLX)**(1.D+0/CL)
          ENDIF
!
!---    Brooks and Corey saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ---
!
        ELSE
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          IF( (1.D+0-SLX)/EPSL.LT.EPSL ) THEN
            HDGL = SCHRV(1,N)
            GOTO 222
          ENDIF
          CL = MAX( SCHR(3,IZN),SMALL )
          IF( SLX.GT.SCHR(4,IZN) ) THEN
            HDGL = SCHRV(1,N)*(1.D+0/ESLX)**(1.D+0/CL)
          ELSE
            HDGL = SCHR(12,IZN)
            SLX = SLRX + 1.D-9
          ENDIF
  222     CONTINUE
        ENDIF
        CPGL = HDGL*RHORL*GRAV/BTGLX
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
        PORD_MX = (1.D+0-POR(4,IZN))*POR0(2,N)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR0(2,N) + SMALL )
        PORD_FX = POR(4,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR0(2,N) + SMALL )
!
!---    Use matrix properties to generate a guess for 
!       capillary head  ---
!
        IF( SLX.GT.SCHR(4,IZN) ) THEN
          ASLX = (SLX-SCHR(4,IZN))/(1.D+0-SCHR(4,IZN))
          HDGL = (((1.D+0/ASLX)**(1.D+0/CMM)-1.D+0)**(1.D+0/CNM))/
     &      SCHRV(1,N)
        ELSE
          HDGL = SCHR(12,IZN)
          SLX = SCHR(4,IZN) + 1.D-6
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  300   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,IZN)),ZERO )*REALX
!
!---    Matrix saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        DSLRX = -SCHR(4,IZN)/(HDGL*LOG(SCHR(12,IZN)))*REALX
        ASLX = 1.D+0/((1.D+0 + (SCHRV(1,N)*HDGL)**CNM)**CMM)
        DASLX = -CMM*SCHRV(1,N)*CNM*((SCHRV(1,N)*HDGL)**(CNM-1.D+0))
     &  /((1.D+0 + (SCHRV(1,N)*HDGL)**CNM)**(CMM+1.D+0))
        SLMZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLMZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
!
!---    Fracture saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(7,IZN),ZERO )
        DSLRX = -SCHR(7,IZN)/(HDGL*LOG(SCHR(12,IZN)))*REALX
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
          CHMSG = 'Dual Porosity van Genuchten: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &         'Aqueous Saturation = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 300
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Dual porosity Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.4 ) THEN
        IF( (1.D+0-SLX)/EPSL.LT.EPSL ) THEN
          HDGL = SCHRV(1,N)
          GOTO 410
        ENDIF
        CLM = MAX( SCHR(3,IZN),SMALL )
        CLF = MAX( SCHR(6,IZN),SMALL )
        PORD_MX = (1.D+0-POR(4,IZN))*POR0(2,N)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR0(2,N) + SMALL )
        PORD_FX = POR(4,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR0(2,N) + SMALL )
!
!---    Use matrix properties to generate a guess for 
!       capillary head  ---
!
        IF( SLX.GT.SCHR(4,IZN) ) THEN
          ASLX = (SLX-SCHR(4,IZN))/(1.D+0-SCHR(4,IZN))
          HDGL = SCHRV(1,N)*(1.D+0/ASLX)**(1.D+0/CLM)
        ELSE
          HDGL = SCHR(12,IZN)
          SLX = SCHR(4,IZN) + 1.D-9
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  400   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,IZN)),ZERO )*REALX
!
!---    Matrix saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        DSLRX = -SCHR(4,IZN)*REALX/(HDGL*LOG(SCHR(12,IZN)))
        HDGLX = MAX( SCHRV(1,N),HDGL )
        ASLX = (SCHRV(1,N)/HDGLX)**CLM
        DASLX = -CLM*(SCHRV(1,N)/(HDGLX**2))
     &    *(SCHRV(1,N)/HDGLX)**(CLM-1.D+0)
        SLMZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLMZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
!
!---    Fracture saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(7,IZN),ZERO )
        DSLRX = -SCHR(7,IZN)*REALX/(HDGL*LOG(SCHR(12,IZN)))
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
          CHMSG = 'Dual Porosity Brooks and Corey: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &         'Aqueous Saturation = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 400
  410   CONTINUE
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Haverkamp saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.5 ) THEN
        IF( SLX.GT.SCHR(4,IZN) ) THEN
          ASLX = (SLX-SCHR(4,IZN))/(1.D+0-SCHR(4,IZN))
          HDGL = SCHRV(1,N) + SCHR(5,IZN)*
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
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,IZN)),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        DSLRX = -SCHR(4,IZN)*REALX/(HDGL*LOG(SCHR(12,IZN)))
        ASLX = SCHR(2,IZN)/(SCHR(2,IZN)+((HDGL-SCHRV(1,N))/
     &    SCHR(5,IZN))**SCHR(3,IZN))
        DASLX = -(SCHR(2,IZN)*SCHR(3,IZN)*
     &    (((HDGL-SCHRV(1,N))/SCHR(5,IZN))**(SCHR(3,IZN)-1.D+0))
     &    /SCHR(5,IZN))/((SCHR(2,IZN)+((HDGL-SCHRV(1,N))/SCHR(5,IZN))
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
          CHMSG = 'Haverkamp: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &         'Aqueous Saturation = '
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
!---  Log-linear interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.11 ) THEN
        ITBX = 0
        HDGL = FNTBLX( SLX,ISLTBL(1,IZN),ISLTBL(2,IZN),ITBX )
        HDGL = EXP(HDGL)
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Hysteretic linear interpolation function, using the drainage
!     curve as an initial guess  ---
!
      ELSEIF( ISCHR(IZN).EQ.12 ) THEN
        SLDX = MAX( SLX-SGTX,0.D+0 )
        ITBX = 0
        HDGL = FNTBLX( SLDX,ISLTBL(1,IZN),ISLTBL(2,IZN),ITBX )
        CPGL = HDGL*RHORL*GRAV/BTGLX
        
!
!---  Hysteretic log-linear interpolation function, using the drainage
!     curve  ---
!
      ELSEIF( ISCHR(IZN).EQ.13 ) THEN
        SLDX = MAX( SLX-SGTX,0.D+0 )
        ITBX = 0
        HDGL = FNTBLX( SLDX,ISLTBL(1,IZN),ISLTBL(2,IZN),ITBX )
        HDGL = EXP(HDGL)
        CPGL = HDGL*RHORL*GRAV/BTGLX
!!
!!---  Cubic spline interpolation function  ---
!!
!      ELSEIF( ISCHR(IZN).EQ.11 ) THEN
!        HDGL = FSPLNX( SLX,ISLTBL(1,IZN),ISLTBL(2,IZN) )
!        CPGL = HDGL*RHORL*GRAV/BTGLX
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
     &        SCHRV(1,N)
          ENDIF
!
!---    van Genuchten saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          HDGL = (((1.D+0/ASLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))/
     &      SCHRV(1,N)
        ENDIF
        CPGL = HDGL*RHORL*GRAV/BTGLX
!        SLRX = SCHR(4,IZN)
!!
!!---    van Genuchten saturation function w/ gas entrapment, 
!!       w/ Webb extension ASL = SL + SGT  ---
!!
!        IF( ISM(IZN).EQ.2 ) THEN
!!
!!---      Check/correct the trapped gas saturation  ---
!!
!          ASLM = 0.D+0
!          SLRX = SCHR(4,IZN)
!          ESGTMX = SCHR(15,IZN)
!          ASLX = SLX + SGTX
!          IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
!            R = 1.D+0/ESGTMX - 1.D+0
!            ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
!     &        (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
!          ELSE
!            ESGTX = 0.D+0
!          ENDIF
!          SLX = MAX( MIN( SLX,(1.D+0-SGTX) ),0.D+0 )
!          SGTMX = ESGTX
!          SGTX = MIN( SGTX,SGTMX )
!          SLGTX = SLX + SGTX
!          CN = MAX( SCHR(3,IZN),SMALL )
!          IF( SCHR(14,IZN).LE.ZERO ) THEN
!            IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
!              CM = 1.D+0 - 2.D+0/CN
!            ELSE
!              CM = 1.D+0 - 1.D+0/CN
!            ENDIF
!          ELSE
!            CM = SCHR(14,IZN)
!          ENDIF
!          SMPX = SCHR(8,IZN)
!!
!!---      Aqueous saturation below the matching point,
!!         use Webb extension  ---
!!
!          IF( SLGTX.LT.SMPX ) THEN
!            HMPX = SCHR(9,IZN)
!            DMPX = -(LOG10(SCHR(12,IZN))-LOG10(HMPX))/SMPX
!            HDGL = 1.D+1**(DMPX*(SLGTX-SMPX) + LOG10(HMPX))
!!
!!---      Aqueous saturation at or above the matching point,
!!         use van Genuchten function
!!
!          ELSE
!            ASLX = (SLGTX-SLRX)/(1.D+0-SLRX)
!            HDGL = (((1.D+0/ASLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))/
!     &        SCHRV(1,N)
!          ENDIF
!          CPGL = HDGL*RHORL*GRAV/BTGLX
!!
!!---    van Genuchten saturation function w/ gas entrapment, 
!!       w/o Webb extension ASL = ESL + ESGT  ---
!!
!        ELSE
!!
!!---      Check/correct the trapped gas saturation  ---
!!
!          ASLM = 0.D+0
!          SLRX = SCHR(4,IZN)
!          ESGTMX = SCHR(15,IZN)/(1.D+0-SLRX)
!          ASLX = (SLX+SGTX-SLRX)/(1.D+0-SLRX)
!          IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
!            R = 1.D+0/ESGTMX - 1.D+0
!            ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
!     &        (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
!          ELSE
!            ESGTX = 0.D+0
!          ENDIF
!          SLX = MAX( MIN( SLX,(1.D+0-SGTX) ),0.D+0 )
!          SGTMX = ESGTX*(1.D+0-SLRX)
!          SGTX = MIN( SGTX,SGTMX )
!          SLGTX = SLX + SGTX
!          CN = MAX( SCHR(3,IZN),SMALL )
!          IF( SCHR(14,IZN).LE.ZERO ) THEN
!            IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
!              CM = 1.D+0 - 2.D+0/CN
!            ELSE
!              CM = 1.D+0 - 1.D+0/CN
!            ENDIF
!          ELSE
!            CM = SCHR(14,IZN)
!          ENDIF
!          IF( SLGTX.GT.SCHR(4,IZN) ) THEN
!            ASLX = (SLGTX-SLRX)/(1.D+0-SLRX)
!            HDGL = (((1.D+0/ASLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))/
!     &        SCHRV(1,N)
!          ELSE
!            HDGL = SCHR(12,IZN)
!            SLGTX = SLRX + SGTX + 1.D-6
!          ENDIF
!        ENDIF
!        CPGL = HDGL*RHORL*GRAV/BTGLX
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
            HDGL = SCHRV(1,N)*(1.D+0/ESLX)**(1.D+0/CL)
          ENDIF
!
!---    Brooks and Corey saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          HDGL = SCHRV(1,N)*(1.D+0/ASLX)**(1.D+0/CL)
        ENDIF
        CPGL = HDGL*RHORL*GRAV/BTGLX
!!
!!---    Brooks and Corey saturation function w/ gas entrapment, 
!!       w/ Webb extension ASL = SL + SGT  ---
!!
!        IF( ISM(IZN).EQ.2 ) THEN
!!
!!---      Check/correct the trapped gas saturation  ---
!!
!          ASLM = 0.D+0
!          SLRX = SCHR(4,IZN)
!          ESGTMX = SCHR(15,IZN)
!          ASLX = SLX + SGTX
!          IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
!            R = 1.D+0/ESGTMX - 1.D+0
!            ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
!     &        (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
!          ELSE
!            ESGTX = 0.D+0
!          ENDIF
!          SLX = MAX( MIN( SLX,(1.D+0-SGTX) ),0.D+0 )
!          SGTMX = ESGTX
!          SGTX = MIN( SGTX,SGTMX )
!          SLGTX = SLX + SGTX
!          CL = MAX( SCHR(3,IZN),SMALL )
!          SMPX = SCHR(8,IZN)
!!
!!---      Aqueous saturation below the matching point,
!!         use Webb extension  ---
!!
!          IF( SLGTX.LT.SMPX ) THEN
!            HMPX = SCHR(9,IZN)
!            DMPX = -(LOG10(SCHR(12,IZN))-LOG10(HMPX))/SMPX
!            HDGL = 1.D+1**(DMPX*(SLGTX-SMPX) + LOG10(HMPX))
!!
!!---      Aqueous saturation at or above the matching point,
!!         use Brooks and Corey function
!!
!          ELSE
!            ASLX = (SLGTX-SLRX)/(1.D+0-SLRX)
!            HDGL = SCHRV(1,N)*(1.D+0/ASLX)**(1.D+0/CL)
!          ENDIF
!          CPGL = HDGL*RHORL*GRAV/BTGLX
!!
!!---    Brooks and Corey saturation function w/ gas entrapment, 
!!       w/o Webb extension ASL = ESL + ESGT  ---
!!
!        ELSE
!!
!!---      Check/correct the trapped gas saturation  ---
!!
!          ASLM = 0.D+0
!          SLRX = SCHR(4,IZN)
!          ESGTMX = SCHR(15,IZN)/(1.D+0-SLRX)
!          ASLX = (SLX+SGTX-SLRX)/(1.D+0-SLRX)
!          IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
!            R = 1.D+0/ESGTMX - 1.D+0
!            ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
!     &        (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
!          ELSE
!            ESGTX = 0.D+0
!          ENDIF
!          SLX = MAX( MIN( SLX,(1.D+0-SGTX) ),0.D+0 )
!          SGTMX = ESGTX*(1.D+0-SLRX)
!          SGTX = MIN( SGTX,SGTMX )
!          SLGTX = SLX + SGTX
!          IF( (1.D+0-SLGTX)/EPSL.LT.EPSL ) THEN
!            HDGL = SCHRV(1,N)
!            GOTO 722
!          ENDIF
!          CL = MAX( SCHR(3,IZN),SMALL )
!          IF( SLGTX.GT.SCHR(4,IZN) ) THEN
!            ASLX = (SLGTX-SCHR(4,IZN))/(1.D+0-SCHR(4,IZN))
!            HDGL = SCHRV(1,N)*(1.D+0/ASLX)**(1.D+0/CL)
!          ELSE
!            HDGL = SCHR(12,IZN)
!            SLGTX = SCHR(4,IZN) + SGTX + 1.D-9
!          ENDIF
!  722     CONTINUE
!          CPGL = HDGL*RHORL*GRAV/BTGLX
!        ENDIF
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
     &      SCHRV(1,N)
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
          RLMSG = SLX
          CHMSG = 'van Genuchten Drainage-Imbibition: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &         'Aqueous Saturation = '
          CALL WRMSGS( INDX )
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
              ASLDZ = (1.D+0/(1.D+0 + (SCHRV(1,N)*HDGLZ)**CND))**CMD
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
            ASLDZ = (1.D+0/(1.D+0 + (SCHRV(1,N)*HDGLZ)**CND))**CMD
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
        CPGL = HDGL*RHORL*GRAV/BTGLX
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
            HDGLD = SCHRV(1,N)*(1.D+0/ASLZ)**(1.D+0/CLD)
!
!---      w/o Webb extension  ---
!
          ELSE
            HDGLD = SCHRV(1,N)*(1.D+0/ASLX)**(1.D+0/CLD)
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
          RLMSG = SLX
          CHMSG = 'Brooks and Corey Drainage-Imbibition: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &         'Aqueous Saturation = '
          CALL WRMSGS( INDX )
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
              IF( HDGLZ.LE.SCHRV(1,N) ) THEN
                ASLDZ = 1.D+0
              ELSE
                ASLDZ = (SCHRV(1,N)/HDGLZ)**CLD
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
            IF( HDGLZ.LE.SCHRV(1,N) ) THEN
              ASLDZ = 1.D+0
            ELSE
              ASLDZ = (SCHRV(1,N)/HDGLZ)**CLD
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
        CPGL = HDGL*RHORL*GRAV/BTGLX
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Check the thermodynamic and hydrologic states declared through
!     user inputs.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE NAPL
      USE LEAK_WELL
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
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHK_CO2'
      EPSLX = 1.D-4
!
!---  Initialize variable properties  ---
!
      DO 10 N = 1,NFLD+NWN_LW
        IZN = IZ(N)
        POR0(1,N) = POR(1,IZN)
        POR0(2,N) = POR(2,IZN)
        PERMV(1,N) = PERM(1,IZN)
        PERMV(2,N) = PERM(2,IZN)
        PERMV(3,N) = PERM(3,IZN)
        SCHRV(1,N) = SCHR(1,IZN)
        PG(1,N) = 0.D+0
        PL(1,N) = 0.D+0
   10 CONTINUE
!
!---  Set CO2 molecular weight, critical temperature, critical
!     pressure, triple-point pressure, and triple-point 
!     temperature  ---
!
      WTMA = 44.010D+0
      TCRA = 304.1282D+0
      PCRA = 73.773D+5
      PTPA = 517.95D+3
      TTPA = 216.592D+0
      IF( LPTA.LE.0 ) THEN
        INDX = 22
        CHMSG = 'Parameter LPTA = 0 for CO2 Operational Mode'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Read CO2 property file  ---
!
      CALL RDPF_A
!
!---  Set oven-dried head to Webb value (1.e+7 cm)  ---
!
      HDOD = 1.D+5
!
!---  Check initial temperature, aqueous pressure, gas pressure,
!     and aqueous saturation  ---
!
      INDX = 0
      DO 100 N = 1,NFLD+NWN_LW
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
!---    Active field node  ---
!
        IZN = IZ(N)
        N_DB = N
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) GOTO 100
        IF( T(2,N).GT.374.14D+0 .OR. T(2,N).LT.0.01D+0 ) THEN
          INDX = 17
          N_DB = N
          RLMSG = T(2,N)
          CHMSG = 'Out of Range Initial Temperature(C) = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( PL(2,N).GT.8.D+8-PATM ) THEN
          INDX = 17
          N_DB = N
          RLMSG = PL(2,N)+PATM
          CHMSG = 'Out of Range Initial Aqueous Pressure(Pa) = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( PG(2,N).GT.8.D+8-PATM .OR. PG(2,N).LT.6.1125D+2-PATM ) THEN
          INDX = 17
          N_DB = N
          RLMSG = PG(2,N)+PATM
          CHMSG = 'Out of Range Initial Gas Pressure(Pa) = '
          CALL WRMSGS( INDX )
        ENDIF
        IF( SL(2,N).GT.1.D+0 .OR. SL(2,N).LT.0.D+0 ) THEN
          INDX = 17
          N_DB = N
          RLMSG = SL(2,N)
          CHMSG = 'Out of Range Initial Aqueous Saturation = '
          CALL WRMSGS( INDX )
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
              N_DB = N
              RLMSG = RPGC(3,IZN)
              CHMSG = 'Modified-Corey Gas Relative Permeability: ' // 
     &          'Non-zero Residual Gas Saturation with '  //
     &          'Gas Entrapment in Saturation Function'
              CALL WRMSGS( INDX )
            ENDIF
!
!---      Free-Corey gas relative permeability function  ---
!
          ELSEIF( IRPG(IZN).EQ.7 ) THEN
            IF( RPGC(4,IZN).GT.EPSL ) THEN
              INDX = 17
              N_DB = N
              RLMSG = RPGC(4,IZN)
              CHMSG = 'Free-Corey Gas Relative Permeability: ' // 
     &          'Non-zero Residual Gas Saturation with '  //
     &          'Gas Entrapment in Saturation Function'
              CALL WRMSGS( INDX )
            ENDIF         
!
!---       van Genuchten gas relative permeability function  ---
!
          ELSEIF( IRPG(IZN).EQ.9 ) THEN
            IF( RPGC(3,IZN).GT.EPSL ) THEN
              INDX = 17
              N_DB = N
              RLMSG = RPGC(4,IZN)
              CHMSG = 'van Genuchten Gas Relative Permeability: ' // 
     &          'Non-zero Residual Gas Saturation with '  //
     &          'Gas Entrapment in Saturation Function'
              CALL WRMSGS( INDX )
            ENDIF
!
!---      Classical-Corey gas relative permeability function  ---
!
          ELSEIF( IRPG(IZN).EQ.17 ) THEN
            IF( RPGC(4,IZN).GT.EPSL ) THEN
              INDX = 17
              N_DB = N
              RLMSG = RPGC(3,IZN)
              CHMSG = 'Classical-Corey Gas Relative Permeability: ' // 
     &          'Non-zero Residual Gas Saturation with '  //
     &          'Gas Entrapment in Saturation Function'
              CALL WRMSGS( INDX )
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
            ALPHAX = SCHRV(1,N)
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
!---      Brooks and Corey moisture retension function  ---
!
          ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR. 
     &      ISCHR(IZN).EQ.202 ) THEN
            SRX = SCHR(4,IZN)
            PSIX = SCHRV(1,N)
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
          ENDIF
        ENDIF
!#ifdef ecke
!!
!!---    Load reactive transport total and diffusive porosity  ---
!!
!        POR0(1,N) = POR0(1,N)
!        POR0(2,N) = POR0(2,N)
!#endif
  100 CONTINUE
      IF( INDX.GT.0 ) STOP
!
!---  Return for restart simulations  ---
!
      IF( (IEO.EQ.2)
     &  .AND. (ISLC(21).EQ.0 .OR. ISLC(21).EQ.IOM) ) THEN
        ISUB_LOG = ISUB_LOG-1
!
!---    Gas saturation scaling  ---
!
        DO N = 1,NFLD+NWN_LW
          IF( POSM(2,N).GE.0.D+0 ) THEN
            SG(2,N) = POSM(2,N)*SG(2,N)
            SGT(2,N) = POSM(2,N)*SGT(2,N)
            SL(2,N) = 1.D+0 - SG(2,N)
!
!---        Minimum aqueous saturation computed by CAP_CO2  ---
!
            ASLMINX = -1.D+0
            CALL CAP_CO2( ASLMINX,BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
            ASLMIN(2,N) = ASLMINX
            PL(2,N) = PG(2,N)-PCX
          ENDIF
        ENDDO
        RETURN
      ENDIF
!
!---    Establish reference pressure for soil compressibility  ---
!
      DO 200 N = 1,NFLD+NWN_LW
!
!---    Skip inactive nodes  ---
!
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
!---  Hydrostatic initial conditions, establish aqueous pressure,
!     temperature, and salt mass fraction at node 1  ---
!
      IF( ISIC.EQ.4 ) THEN
!
!---    Assign properties at reference pressure elevation  ---
!
        PLX = VSLC(1) + PATM
        ZPX = VSLC(2)
        PGX = PLX
        TX = VSLC(3) + (ZPX-VSLC(4))*VSLC(5)
        YLSX = VSLC(6) + (ZPX-VSLC(7))*VSLC(8)
        XLAX = 0.D+0
!
!---    Check for out-of-range salt mass fraction  ---
!
        CALL SOL_LS( TX,XLSMX )
        IF( YLSX.LT.0.D+0 ) THEN
          INDX = 17
          N_DB = 1
          RLMSG = YLSX
          CHMSG = 'Hydrostatic Initial Condition: Salt ' //
     &      'Mass Fraction < 0.0 : XLS = '
           CALL WRMSGS( INDX )
        ELSEIF( YLSX.GT.XLSMX ) THEN
          INDX = 17
          N_DB = 1
          RLMSG = YLSX
          CHMSG = 'Hydrostatic Initial Condition: Salt ' //
     &      'Mass Fraction > Solubility Limit : XLS = '
           CALL WRMSGS( INDX )
        ENDIF
        PX = MAX( PGX,PLX )
        CHMSGX(1) = 'Unconverged Hydrostatic Initial Conditions: ' //
     &    'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
        CHMSGX(2) = 'Hydrostatic Initial Condition Transition: ' //
     &    'Vapor Pressure > Gas Pressure: PVB + PVA = '
        CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
        CALL DENS_B( TX,PLX,XLSX,RHOBX )
        CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Hydrostatic initial conditions, with salt mass fraction
!     read from an external file  ---
!
      ELSEIF( ISIC.EQ.5 ) THEN
!
!---    Locate the node nearest to the reference pressure location  ---
!
        DISTX = 1.D+20
        NX = 0
        DO 220 N = 1,NFLD+NWN_LW
          DX = SQRT( (XP(N)-VSLC(6))**2 + (YP(N)-VSLC(7))**2 + 
     &      (ZP(N)-VSLC(8))**2 )
          IF( DX.LT.DISTX ) THEN
            DISTX = DX
            NX = N
          ENDIF
  220   CONTINUE
!
!---    Check for reference pressure location within nearest node  ---
!
        CALL WITHIN( VSLC(6),VSLC(7),VSLC(8),ICWX,NX )
        IF( ICWX.EQ.0 ) THEN
         INDX = 24
         IMSG = NX
         CHMSG = 'Hydrostatic Initial Condition: Reference Pressure' //
     &     'Coordinates Outside of Nearest Node: '
         CALL WRMSGS( INDX )
        ENDIF
!
!---    Translate reference pressure to nearest node  ---
!
        PLX = VSLC(1) + PATM
        PGX = PLX
        XLAX = 0.D+0
        TX = VSLC(3) + (ZP(NX)-VSLC(4))*VSLC(5)
        YLSX = YLS(2,NX)
!
!---    Check for out-of-range salt mass fraction  ---
!
        CALL SOL_LS( TX,XLSMX )
        IF( YLSX.LT.0.D+0 ) THEN
          INDX = 17
          N_DB = NX
          RLMSG = YLSX
          CHMSG = 'Hydrostatic Initial Condition: Salt ' //
     &      'Mass Fraction < 0.0 : XLS = '
          CALL WRMSGS( INDX )
        ELSEIF( YLSX.GT.XLSMX ) THEN
          INDX = 17
          N_DB = NX
          RLMSG = YLSX
          CHMSG = 'Hydrostatic Initial Condition: Salt ' //
     &      'Mass Fraction > Solubility Limit : XLS = '
           CALL WRMSGS( INDX )
        ENDIF
        PX = MAX( PGX,PLX )
        CHMSGX(1) = 'Unconverged Hydrostatic Initial Conditions: ' //
     &    'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
        CHMSGX(2) = 'Hydrostatic Initial Condition Transition: ' //
     &    'Vapor Pressure > Gas Pressure: PVB + PVA = '
        CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
        CALL DENS_B( TX,PLX,XLSX,RHOBX )
        CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
        RHOLSX = RHOLX
        PLX = PLX + RHOLX*((VSLC(8)-ZP(NX))*GRAVZ + 
     &    (VSLC(7)-YP(NX))*GRAVY + (VSLC(6)-XP(NX))*GRAVX)
        PGX = PLX
        PL(2,NX) = PLX - PATM
!
!---    Loop from nearest node to domain bottom  ---
!
        NPX = NX
        RHOLPX = RHOLSX
        DO 240 N = NX,1,-1
!
!---      Skip inactive nodes  ---
!
          IF( IXP(N).EQ.0 ) GOTO 240
          IZN = IZ(N)
!
!---      Assign gas-entry pressure for non Brooks-Corey;
!         Brooks-Corey; Brooks-Corey, Dual Porosity; and
!         Brooks-Corey, Entrapment  ---
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
          TX = VSLC(3) + (ZP(N)-VSLC(4))*VSLC(5)
          YLSX = YLS(2,N)
          XLAX = 0.D+0
!
!---      Check for out-of-range salt mass fraction  ---
!
          CALL SOL_LS( TX,XLSMX )
          IF( YLSX.LT.0.D+0 ) THEN
            INDX = 17
            N_DB = N
            RLMSG = YLSX
            CHMSG = 'Hydrostatic Initial Condition: Salt ' //
     &        'Mass Fraction < 0.0 : XLS = '
            CALL WRMSGS( INDX )
          ELSEIF( YLSX.GT.XLSMX ) THEN
            INDX = 17
            N_DB = N
            RLMSG = YLSX
            CHMSG = 'Hydrostatic Initial Condition: Salt ' //
     &        'Mass Fraction > Solubility Limit : XLS = '
             CALL WRMSGS( INDX )
          ENDIF
          PX = PLX
          CHMSGX(1) = 'Unconverged Hydrostatic Initial ' //
     &      'Conditions: CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
          CHMSGX(2) = 'Hydrostatic Initial Condition Transition: ' //
     &      'Vapor Pressure > Gas Pressure: PVB + PVA = '
          CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
          CALL DENS_B( TX,PLX,XLSX,RHOBX )
          CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
          PLX = PLX - (ZP(N)-ZP(NPX))*GRAVZ*
     &          (RHOLPX*DXGF(N)+RHOLX*DXGF(NPX))/(DXGF(N)+DXGF(NPX))
          PL(2,N) = PLX - PATM
          PG(2,N) = PL(2,N)+ENPR-EPSLX
          T(2,N) = TX
          YLS(2,N) = YLSX
          XLS(2,N) = XLSX
          XLA(2,N) = XLAX
!
!---      Assign initial phase condition  ---
!
          NPHAZ(2,N) = 1
!
!---      Small concentration limits  ---
!
          IF( YLS(2,N).LT.1.D-15 ) YLS(2,N) = 0.D+0
          IF( XLS(2,N).LT.1.D-15 ) XLS(2,N) = 0.D+0
          IF( XLA(2,N).LT.1.D-15 ) XLA(2,N) = 0.D+0
!
!---      Initialize reference pressure for compressibility  ---
!
          IF( CMP(3,IZN).GT.PATM ) THEN
            PCMP(N) = CMP(3,IZN)
          ELSEIF( ISLC(61).EQ.0 ) THEN
            PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
          ENDIF
!#ifdef ecke
!          POR0(1,N) = POR0(1,N)
!          POR0(2,N) = POR0(2,N)
!#endif
          NPX = N
          RHOLPX = RHOLX
  240   CONTINUE
!
!---    Loop from nearest node to domain top  ---
!
        NPX = NX
        RHOLPX = RHOLSX
        DO 250 N = NX+1,NFLD
!
!---      Skip inactive nodes  ---
!
          IF( IXP(N).EQ.0 ) GOTO 250
          IZN = IZ(N)
!
!---      Assign gas-entry pressure for non Brooks-Corey;
!         Brooks-Corey; Brooks-Corey, Dual Porosity; and
!         Brooks-Corey, Entrapment  ---
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
          TX = VSLC(3) + (ZP(N)-VSLC(4))*VSLC(5)
          YLSX = YLS(2,N)
          XLAX = 0.D+0
!
!---      Check for out-of-range salt mass fraction  ---
!
          CALL SOL_LS( TX,XLSMX )
          IF( YLSX.LT.0.D+0 ) THEN
            INDX = 17
            N_DB = N
            RLMSG = YLSX
            CHMSG = 'Hydrostatic Initial Condition: Salt ' //
     &        'Mass Fraction < 0.0 : XLS = '
            CALL WRMSGS( INDX )
          ELSEIF( YLSX.GT.XLSMX ) THEN
            INDX = 17
            N_DB = N
            RLMSG = YLSX
            CHMSG = 'Hydrostatic Initial Condition: Salt ' //
     &        'Mass Fraction > Solubility Limit : XLS = '
             CALL WRMSGS( INDX )
          ENDIF
          PX = PLX
          CHMSGX(1) = 'Unconverged Hydrostatic Initial ' //
     &      'Conditions: CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
          CHMSGX(2) = 'Hydrostatic Initial Condition Transition: ' //
     &      'Vapor Pressure > Gas Pressure: PVB + PVA = '
          CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
          CALL DENS_B( TX,PLX,XLSX,RHOBX )
          CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
          PLX = PLX - (ZP(N)-ZP(NPX))*GRAVZ*
     &          (RHOLPX*DXGF(N)+RHOLX*DXGF(NPX))/(DXGF(N)+DXGF(NPX))
          PL(2,N) = PLX - PATM
          PG(2,N) = PL(2,N)+ENPR-EPSLX
          T(2,N) = TX
          YLS(2,N) = YLSX
          XLS(2,N) = XLSX
          XLA(2,N) = XLAX
!
!---      Assign initial phase condition  ---
!
          NPHAZ(2,N) = 1
!
!---      Small concentration limits  ---
!
          IF( YLS(2,N).LT.1.D-15 ) YLS(2,N) = 0.D+0
          IF( XLS(2,N).LT.1.D-15 ) XLS(2,N) = 0.D+0
          IF( XLA(2,N).LT.1.D-15 ) XLA(2,N) = 0.D+0
!
!---      Initialize reference pressure for compressibility  ---
!
          IF( CMP(3,IZN).GT.PATM ) THEN
            PCMP(N) = CMP(3,IZN)
          ELSEIF( ISLC(61).EQ.0 ) THEN
            PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
          ENDIF
!#ifdef ecke
!          POR0(1,N) = POR0(1,N)
!          POR0(2,N) = POR0(2,N)
!#endif
          NPX = N
          RHOLPX = RHOLX
  250   CONTINUE
        GOTO 3010
      ENDIF
!
!---  Convert initial aqueous-salt and aqueous-CO2 inputs into
!     aqueous mass fractions, considering initial saturation
!     and surface tension affects  ---
!
      DO 3000 N = 1,NFLD+NWN_LW
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 .AND. ISIC.NE.4 ) GOTO 3000
        IZN = IZ(N)
        N_DB = N
!
!---    Assign gas-entry pressure for non Brooks-Corey;
!       Brooks-Corey; Brooks-Corey, Dual Porosity; and
!       Brooks-Corey, Entrapment  ---
!
        IF( IXP(N).EQ.0 ) THEN
          ENPR = 0.D+0
          ESGTMX = 0.D+0
        ELSE
          IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &      ISCHR(IZN).EQ.202 ) THEN
            ENPR = SCHRV(1,N)*RHORL*GRAV
          ELSEIF( ISCHR(IZN).EQ.4 ) THEN
            ENPR = MIN( SCHRV(1,N),SCHR(5,IZN) )*RHORL*GRAV
          ELSE
            ENPR = 0.D+0
          ENDIF
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
        ENDIF
!
!---    Hydrostatic initial conditions  ---
!
        IF( ISIC.EQ.4 ) THEN
          IF( ABS(ZP(N)-ZPX).LT.EPSL ) THEN
            PL(2,N) = PLX - PATM
            PG(2,N) = PL(2,N)+ENPR-EPSLX
            T(2,N) = TX
            XLA(2,N) = 0.D+0
            ICBRN(N) = 3
            ICAIR(N) = 3
!
!---        Assign initial phase condition  ---
!
            NPHAZ(2,N) = 1
!
!---        Aqueous-CO2 mass fraction  ---
!
            YLS(2,N) = YLSX
            XLS(2,N) = XLSX
            RHOL(2,N) = RHOLX
          ELSE 
            PLX = VSLC(1) + PATM
            KC = MAX( 1,INT(ABS(ZP(N)-VSLC(2))) )
            DISTZ = (ZP(N)-VSLC(2))/REAL(KC)
            ZX = VSLC(2)
            XLAX = 0.D+0
            DO M = 1,KC
              ZX = ZX + 5.D-1*DISTZ
              TX = VSLC(3) + (ZX-VSLC(4))*VSLC(5)
              YLSX = VSLC(6) + (ZX-VSLC(7))*VSLC(8)
              ZX = ZX + 5.D-1*DISTZ
!
!---          Check for out-of-range salt mass fraction  ---
!
              CALL SOL_LS( TX,XLSMX )
              IF( YLSX.LT.0.D+0 ) THEN
                INDX = 17
                N_DB = N
                RLMSG = YLSX
                CHMSG = 'Hydrostatic Initial Condition: Salt ' //
     &            'Mass Fraction < 0.0 : XLS = '
                 CALL WRMSGS( INDX )
              ELSEIF( YLSX.GT.XLSMX ) THEN
                INDX = 17
                N_DB = N
                RLMSG = YLSX
                CHMSG = 'Hydrostatic Initial Condition: Salt ' //
     &            'Mass Fraction > Solubility Limit : XLS = '
                 CALL WRMSGS( INDX )
              ENDIF
              PX = PLX
              CHMSGX(1) = 'Unconverged Hydrostatic Initial Cond' //
     &          'itions: CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
              CHMSGX(2) = 'Hydrostatic Initial Condition Transi' //
     &          'tion: Vapor Pressure > Gas Pressure: PVB + PVA = '
              CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
              CALL DENS_B( TX,PLX,XLSX,RHOBX )
              CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
              PLX = PLX - RHOLX*GRAVZ*DISTZ
            ENDDO
            PL(2,N) = PLX - PATM
            PG(2,N) = PL(2,N)+ENPR-EPSLX
            T(2,N) = VSLC(3) + (ZP(N)-VSLC(4))*VSLC(5)
            YLS(2,N) = VSLC(6) + (ZP(N)-VSLC(7))*VSLC(8)
            XLA(2,N) = 0.D+0
            ICBRN(N) = 3
            ICAIR(N) = 3
!
!---        Assign initial phase condition  ---
!
            NPHAZ(2,N) = 1
!
!---        Aqueous-CO2 mass fraction  ---
!
            PLX = PL(2,N)+PATM
            PGX = PG(2,N)+PATM
            PX = MAX( PGX,PLX )
            CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &        'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &        'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_33( T(2,N),PX,PGX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
            CALL DENS_B( T(2,N),PLX,XLS(2,N),RHOBX )
            CALL DENS_L( T(2,N),RHOBX,XLA(2,N),RHOL(2,N) )
!
!---        Assign properties for current elevation  ---
!
            ZPX = ZP(N)
            PLX = PL(2,N) + PATM
            TX = T(2,N)
            YLSX = YLS(2,N)
            XLSX =  XLS(2,N)
            RHOLX = RHOL(2,N)
          ENDIF
!
!---    Unsaturated conditions  ---
!
        ELSEIF( (ISIC.EQ.1 .AND. (1.D+0-SL(2,N)).GT.1.D-6) .OR.
     &      (ISIC.EQ.2 .AND. (1.D+0-SL(2,N)).GT.1.D-6) .OR.
     &      (ISIC.EQ.3 .AND. (PG(2,N)-PL(2,N)-ENPR).GT.1.D-6) ) THEN
!
!---    Gas pressure and aqueous saturation  ---
!
        IF( ISIC.EQ.1 ) THEN
!
!---      Aqueous-salt concentration  ---
!
          IF( ICBRN(N).EQ.1 ) THEN
            RHOLSX = TMS(2,N)
            PGX = PG(2,N)+PATM
            ISRX = 1
            CALL DENS_W( T(2,N),PGX,RHOLWX,RHOGWX,ISRX )
!
!---        Guess aqueous-salt mass fractions  ---
!
            YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
            CALL SOL_LS( T(2,N),XLSMX )
            DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLS(2,N)))
            NC = 0
!
!---        Single-variable Newton-Raphson loop (YLS)  ---
!
  310       CONTINUE
            NC = NC + 1
            IF( NC.GT.32 ) THEN
              INDX = 17
              RLMSG = RHOLSX
              N_DB = N
              CHMSG = 'Unconverged Initial Conditions: Gas Press., ' //
     &          'Aqu. Sat., Salt Aqu. Conc. = '
              CALL WRMSGS( INDX )
            ENDIF
            DO 320 M = 2,3
              YLS(M,N) = YLS(2,N)
              IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
              XLSX = MIN( YLS(M,N),XLSMX )
              XLS(M,N) = YLS(M,N)
              CALL SFT_L( T(2,N),XLSX,SFTLX )
              BTGL(M,N) = 1.D+0
              IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &          BTGL(M,N) = SCHR(16,IZN)/SFTLX
!
!---          Minimum aqueous saturation computed by CAP_CO2  ---
!
              ASLMINX = -1.D+0
              CALL CAP_CO2( ASLMINX,BTGL(M,N),PCX,SL(2,N),SGT(2,N),N )
              CALL SP_B( T(2,N),XLSX,PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
              CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
              IF( PVBX.GT.PGX .AND. M.EQ.2 ) THEN
                INDX = 17
                N_DB = N
                RLMSG = PVBX
                CHMSG = 'Initial Condition Transition: Vapor ' //
     &            'Pressure > Gas Pressure: PVB = '
                CALL WRMSGS( INDX )
              ENDIF
              PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
              XLSSX = XLS(M,N)
              CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &          XMLSX,XMLWX )
              XLA(M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
              XLS(M,N) = XLS(M,N) + (XLSSX-XLS(M,N))*(XLA(M,N)/XLAX)
              CALL DENS_L( T(2,N),RHOBX,XLA(M,N),RHOLX )
              GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
  320       CONTINUE
            RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
            RPX(1) = -GX(2,1)
            CYLSX = RPX(1)/RX(1,1)
            YLS(2,N) = MAX( YLS(2,N)+CYLSX,0.D+0 )
            IF( ABS(CYLSX).GT.(1.D-9*XLSMX) ) GOTO 310
!
!---        Limit salt mass fraction  ---
!
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &          BTGL(2,N) = SCHR(16,IZN)/SFTLX
!
!---        Minimum aqueous saturation computed by CAP_CO2  ---
!
            ASLMINX = -1.D+0
            CALL CAP_CO2( ASLMINX,BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
            PL(2,N) = PG(2,N)-PCX
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVB = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
!
!---      Aqueous-salt relative saturation  ---
!
          ELSEIF( ICBRN(N).EQ.2 ) THEN
            PHILSX = TMS(2,N)
            PGX = PG(2,N)+PATM
!
!---        Guess aqueous-salt mass fractions  ---
!
            CALL SOL_LS( T(2,N),XLSMX )
            YLS(2,N) = PHILSX*XLSMX
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
!
!---        Minimum aqueous saturation computed by CAP_CO2  ---
!
            ASLMINX = -1.D+0
            CALL CAP_CO2( ASLMINX,BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
            PL(2,N) = PG(2,N)-PCX
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX .AND. M.EQ.2 ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVBX = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
!
!---      Aqueous-salt mass fraction  ---
!
          ELSEIF( ICBRN(N).EQ.3 .OR. ICBRN(N).EQ.0 ) THEN
            PGX = PG(2,N)+PATM
!
!---        Guess aqueous-salt mass fractions  ---
!
            CALL SOL_LS( T(2,N),XLSMX )
            YLS(2,N) = TMS(2,N)
            DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLS(2,N)))
            NC = 0
!
!---        Single-variable Newton-Raphson loop (YLS)  ---
!
  350       CONTINUE
            NC = NC + 1
            IF( NC.GT.32 ) THEN
              INDX = 17
              RLMSG = TMS(2,N)
              N_DB = N
              CHMSG = 'Unconverged Initial Conditions: Gas Press., ' //
     &          'Aqu. Sat., Salt Aqu. Mass Frac. = '
              CALL WRMSGS( INDX )
            ENDIF
            DO 360 M = 2,3
              YLS(M,N) = YLS(2,N)
              IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
              XLSX = MIN( YLS(M,N),XLSMX )
              XLS(M,N) = YLS(M,N)
              CALL SFT_L( T(2,N),XLSX,SFTLX )
              BTGL(M,N) = 1.D+0
              IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &          BTGL(M,N) = SCHR(16,IZN)/SFTLX
!
!---          Minimum aqueous saturation computed by CAP_CO2  ---
!
              ASLMINX = -1.D+0
              CALL CAP_CO2( ASLMINX,BTGL(M,N),PCX,SL(2,N),SGT(2,N),N )
              CALL SP_B( T(2,N),XLSX,PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
              CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
              IF( PVBX.GT.PGX .AND. M.EQ.2 ) THEN
                INDX = 17
                N_DB = N
                RLMSG = PVBX
                CHMSG = 'Initial Condition Transition: Vapor ' //
     &            'Pressure > Gas Pressure: PVBX = '
                CALL WRMSGS( INDX )
              ENDIF
              PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
              XLSSX = XLS(M,N)
              CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &          XMLSX,XMLWX )
              XLA(M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
              XLS(M,N) = XLS(M,N) + (XLSSX-XLS(M,N))*(XLA(M,N)/XLAX)
              CALL DENS_L( T(2,N),RHOBX,XLA(M,N),RHOLX )
              GX(M,1) = XLS(M,N) - TMS(2,N)
  360       CONTINUE
            RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
            RPX(1) = -GX(2,1)
            CYLSX = RPX(1)/RX(1,1)
            YLS(2,N) = MAX( YLS(2,N)+CYLSX,0.D+0 )
            IF( ABS(CYLSX).GT.(1.D-9*XLSMX) ) GOTO 350
!
!---        Limit salt mass fraction  ---
!
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
!
!---        Minimum aqueous saturation computed by CAP_CO2  ---
!
            ASLMINX = -1.D+0
            CALL CAP_CO2( ASLMINX,BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
            PL(2,N) = PG(2,N)-PCX
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVBX = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
!
!---      Aqueous-salt molality (mol salt/kg water)  ---
!
          ELSEIF( ICBRN(N).EQ.4 ) THEN
            PGX = PG(2,N)+PATM
            YLS(2,N) = TMS(2,N)*WTMS/(1.D+0 + TMS(2,N)*WTMS)
            CALL SOL_LS( T(2,N),XLSMX )
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
!
!---        Minimum aqueous saturation computed by CAP_CO2  ---
!
            ASLMINX = -1.D+0
            CALL CAP_CO2( ASLMINX,BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVBX = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
            PL(2,N) = PG(2,N)-PCX
          ENDIF
!
!---      Assign initial phase condition  ---
!
          CALL DENS_A( T(2,N),PGAX,RHOAX,I_VX )
          NPHAZ(2,N) = (I_VX+1)*2
!
!---    Aqueous pressure and aqueous saturation  ---
!
        ELSEIF( ISIC.EQ.2 ) THEN
!
!---      Aqueous-salt concentration  ---
!
          IF( ICBRN(N).EQ.1 ) THEN
            RHOLSX = TMS(2,N)
            CALL SP_W( T(2,N),PSWX )
            ISRX = 1
            CALL DENS_W( T(2,N),PSWX,RHOLWX,RHOGWX,ISRX )
!
!---        Guess aqueous-salt mass fractions  ---
!
            YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
            CALL SOL_LS( T(2,N),XLSMX )
            DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLS(2,N)))
            NC = 0
!
!---        Single-variable Newton-Raphson loop (YLS)  ---
!
  410       CONTINUE
            NC = NC + 1
            IF( NC.GT.32 ) THEN
              INDX = 17
              RLMSG = RHOLSX
              N_DB = N
              CHMSG = 'Unconverged Initial Conditions: Aqu. Press., ' //
     &          'Aqu. Sat., Salt Aqu. Conc. = '
              CALL WRMSGS( INDX )
            ENDIF
            DO 420 M = 2,3
              YLS(M,N) = YLS(2,N)
              IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
              XLSX = MIN( YLS(M,N),XLSMX )
              XLS(M,N) = YLS(M,N)
              CALL SFT_L( T(2,N),XLSX,SFTLX )
              BTGL(M,N) = 1.D+0
              IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &          BTGL(M,N) = SCHR(16,IZN)/SFTLX
!
!---          Minimum aqueous saturation computed by CAP_CO2  ---
!
              ASLMINX = -1.D+0
              CALL CAP_CO2( ASLMINX,BTGL(M,N),PCX,SL(2,N),SGT(2,N),N )
              PGX = PL(2,N)+PCX+PATM
              CALL SP_B( T(2,N),XLSX,PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
              CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
              IF( PVBX.GT.PGX .AND. M.EQ.2 ) THEN
                INDX = 17
                N_DB = N
                RLMSG = PVBX
                CHMSG = 'Initial Condition Transition: Vapor ' //
     &            'Pressure > Gas Pressure: PVBX = '
                CALL WRMSGS( INDX )
              ENDIF
              PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
              XLSSX = XLS(M,N)
              CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &          XMLSX,XMLWX )
              XLA(M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
              XLS(M,N) = XLS(M,N) + (XLSSX-XLS(M,N))*(XLA(M,N)/XLAX)
              CALL DENS_L( T(2,N),RHOBX,XLA(M,N),RHOLX )
              GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
  420       CONTINUE
            RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
            RPX(1) = -GX(2,1)
            CYLSX = RPX(1)/RX(1,1)
            YLS(2,N) = MAX( YLS(2,N)+CYLSX,0.D+0 )
            IF( ABS(CYLSX).GT.(1.D-9*XLSMX) ) GOTO 410
!
!---        Limit salt mass fraction  ---
!
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
!
!---        Minimum aqueous saturation computed by CAP_CO2  ---
!
            ASLMINX = -1.D+0
            CALL CAP_CO2( ASLMINX,BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
            PG(2,N) = PL(2,N)+PCX
            PGX = PL(2,N)+PCX+PATM
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVBX = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
!
!---      Aqueous-salt relative saturation  ---
!
          ELSEIF( ICBRN(N).EQ.2 ) THEN
            PHILSX = TMS(2,N)
!
!---        Guess aqueous-salt mass fractions  ---
!
            CALL SOL_LS( T(2,N),XLSMX )
            YLS(2,N) = PHILSX*XLSMX
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
!
!---        Minimum aqueous saturation computed by CAP_CO2  ---
!
            ASLMINX = -1.D+0
            CALL CAP_CO2( ASLMINX,BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
            PG(2,N) = PL(2,N)+PCX
            PGX = PL(2,N)+PCX+PATM
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX .AND. M.EQ.2 ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVBX = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
!
!---      Aqueous-salt mass fraction  ---
!
          ELSEIF( ICBRN(N).EQ.3 .OR. ICBRN(N).EQ.0 ) THEN
!
!---        Guess aqueous-salt mass fractions  ---
!
            CALL SOL_LS( T(2,N),XLSMX )
            YLS(2,N) = TMS(2,N)
            DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLS(2,N)))
            NC = 0
!
!---        Single-variable Newton-Raphson loop (YLS)  ---
!
  450       CONTINUE
            NC = NC + 1
            IF( NC.GT.32 ) THEN
              INDX = 17
              RLMSG = TMS(2,N)
              N_DB = N
              CHMSG = 'Unconverged Initial Conditions: Aqu. Press., ' //
     &          'Aqu. Sat., Salt Aqu. Mass Frac. = '
              CALL WRMSGS( INDX )
            ENDIF
            DO 460 M = 2,3
              YLS(M,N) = YLS(2,N)
              IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
              XLSX = MIN( YLS(M,N),XLSMX )
              XLS(M,N) = YLS(M,N)
              CALL SFT_L( T(2,N),XLSX,SFTLX )
              BTGL(M,N) = 1.D+0
              IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &          BTGL(M,N) = SCHR(16,IZN)/SFTLX
!
!---          Minimum aqueous saturation computed by CAP_CO2  ---
!
              ASLMINX = -1.D+0
              CALL CAP_CO2( ASLMINX,BTGL(M,N),PCX,SL(2,N),SGT(2,N),N )
              PGX = PL(2,N)+PCX+PATM
              CALL SP_B( T(2,N),XLSX,PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
              CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
              IF( PVBX.GT.PGX .AND. M.EQ.2 ) THEN
                INDX = 17
                N_DB = N
                RLMSG = PVBX
                CHMSG = 'Initial Condition Transition: Vapor ' //
     &            'Pressure > Gas Pressure: PVBX = '
                CALL WRMSGS( INDX )
              ENDIF
              PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
              XLSSX = XLS(M,N)
              CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &          XMLSX,XMLWX )
              XLA(M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
              XLS(M,N) = XLS(M,N) + (XLSSX-XLS(M,N))*(XLA(M,N)/XLAX)
              CALL DENS_L( T(2,N),RHOBX,XLA(M,N),RHOLX )
              GX(M,1) = XLS(M,N) - TMS(2,N)
  460       CONTINUE
            RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
            RPX(1) = -GX(2,1)
            CYLSX = RPX(1)/RX(1,1)
            YLS(2,N) = MAX( YLS(2,N)+CYLSX,0.D+0 )
            IF( ABS(CYLSX).GT.(1.D-9*XLSMX) ) GOTO 450
!
!---        Limit salt mass fraction  ---
!
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
!
!---        Minimum aqueous saturation computed by CAP_CO2  ---
!
            ASLMINX = -1.D+0
            CALL CAP_CO2( ASLMINX,BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
            PG(2,N) = PL(2,N)+PCX
            PGX = PL(2,N)+PCX+PATM
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVBX = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
!
!---      Aqueous-salt molality (mol salt/kg water)  ---
!
          ELSEIF( ICBRN(N).EQ.4 ) THEN
            YLS(2,N) = TMS(2,N)*WTMS/(1.D+0 + TMS(2,N)*WTMS)
            CALL SOL_LS( T(2,N),XLSMX )
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
!
!---        Minimum aqueous saturation computed by CAP_CO2  ---
!
            ASLMINX = -1.D+0
            CALL CAP_CO2( ASLMINX,BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
            PGX = PL(2,N)+PCX+PATM
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX .AND. M.EQ.2 ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVBX = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
            PG(2,N) = PL(2,N)+PCX
          ENDIF
!
!---      Assign initial phase condition  ---
!
          CALL DENS_A( T(2,N),PGAX,RHOAX,I_VX )
          NPHAZ(2,N) = (I_VX+1)*2
!
!---    Gas pressure and aqueous pressure  ---
!
        ELSEIF( ISIC.EQ.3 ) THEN
!
!---      Aqueous-salt concentration  ---
!
          IF( ICBRN(N).EQ.1 ) THEN
            RHOLSX = TMS(2,N)
            PGX = PG(2,N)+PATM
            PCX = PG(2,N)-PL(2,N)
            ISRX = 1
            CALL DENS_W( T(2,N),PGX,RHOLWX,RHOGWX,ISRX )
!
!---        Guess aqueous-salt mass fractions  ---
!
            YLS(2,N) = RHOLSX/(RHOLWX+RHOLSX)
            CALL SOL_LS( T(2,N),XLSMX )
            DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLS(2,N)))
            NC = 0
!
!---        Single-variable Newton-Raphson loop (YLS)  ---
!
  510       CONTINUE
            NC = NC + 1
            IF( NC.GT.32 ) THEN
              INDX = 17
              RLMSG = RHOLSX
              N_DB = N
              CHMSG = 'Unconverged Initial Conditions: Gas. Press., ' //
     &          'Aqu. Press., Salt Aqu. Conc. = '
              CALL WRMSGS( INDX )
            ENDIF
            DO 520 M = 2,3
              YLS(M,N) = YLS(2,N)
              IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
              XLSX = MIN( YLS(M,N),XLSMX )
              XLS(M,N) = YLS(M,N)
              CALL SP_B( T(2,N),XLSX,PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
              CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
              IF( PVBX.GT.PGX .AND. M.EQ.2 ) THEN
                INDX = 17
                N_DB = N
                RLMSG = PVBX
                CHMSG = 'Initial Condition Transition: Vapor ' //
     &            'Pressure > Gas Pressure: PVBX = '
                CALL WRMSGS( INDX )
              ENDIF
              PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
              XLSSX = XLS(M,N)
              CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &          XMLSX,XMLWX )
              XLA(M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
              XLS(M,N) = XLS(M,N) + (XLSSX-XLS(M,N))*(XLA(M,N)/XLAX)
              CALL DENS_L( T(2,N),RHOBX,XLA(M,N),RHOLX )
              GX(M,1) = XLS(M,N) - RHOLSX/RHOLX
  520       CONTINUE
            RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
            RPX(1) = -GX(2,1)
            CYLSX = RPX(1)/RX(1,1)
            YLS(2,N) = MAX( YLS(2,N)+CYLSX,0.D+0 )
            IF( ABS(CYLSX).GT.(1.D-9*XLSMX) ) GOTO 510
!
!---        Limit salt mass fraction  ---
!
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
            INDX = 0
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVBX = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
            IF( SGT(2,N).GT.EPSL ) THEN
              SG(2,N) = SGT(2,N)
              INDX = 1
            ENDIF
            CALL KSP_CO2( N,PG(2,N),PL(2,N),SG(2,N),SGT(2,N),SL(2,N),
     &        SDPF(N),SDPM(N),RKL(1,2,N),RKG(2,N),ASL(N),ASLMIN(2,N),
     &        ASGT(N),ESGTMX,SLRX,BTGL(2,N),INDX )
!
!---      Aqueous-salt relative saturation  ---
!
          ELSEIF( ICBRN(N).EQ.2 ) THEN
            PHILSX = TMS(2,N)
            PGX = PG(2,N)+PATM
            PCX = PG(2,N)-PL(2,N)
            ISRX = 1
            CALL DENS_W( T(2,N),PGX,RHOLWX,RHOGWX,ISRX )
!
!---        Guess aqueous-salt mass fractions  ---
!
            CALL SOL_LS( T(2,N),XLSMX )
            YLS(2,N) = PHILSX*XLSMX
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX .AND. M.EQ.2 ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVBX = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
            INDX = 0
            IF( SGT(2,N).GT.EPSL ) THEN
              SG(2,N) = SGT(2,N)
              INDX = 1
            ENDIF
            CALL KSP_CO2( N,PG(2,N),PL(2,N),SG(2,N),SGT(2,N),SL(2,N),
     &        SDPF(N),SDPM(N),RKL(1,2,N),RKG(2,N),ASL(N),ASLMIN(2,N),
     &        ASGT(N),ESGTMX,SLRX,BTGL(2,N),INDX )

!---      Aqueous-salt mass fraction  ---
!
          ELSEIF( ICBRN(N).EQ.3 .OR. ICBRN(N).EQ.0 ) THEN
            PGX = PG(2,N)+PATM
            PCX = PG(2,N)-PL(2,N)
            ISRX = 1
            CALL DENS_W( T(2,N),PGX,RHOLWX,RHOGWX,ISRX )
!
!---        Guess aqueous-salt mass fractions  ---
!
            CALL SOL_LS( T(2,N),XLSMX )
            YLS(2,N) = TMS(2,N)
            DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLS(2,N)))
            NC = 0
!
!---        Single-variable Newton-Raphson loop (YLS)  ---
!
  550       CONTINUE
            NC = NC + 1
            IF( NC.GT.32 ) THEN
              INDX = 17
              RLMSG = TMS(2,N)
              N_DB = N
              CHMSG = 'Unconverged Initial Conditions: Gas. Press., ' //
     &          'Aqu. Press., Salt Aqu. Mass Frac. = '
              CALL WRMSGS( INDX )
            ENDIF
            DO 560 M = 2,3
              YLS(M,N) = YLS(2,N)
              IF( M.EQ.3 ) YLS(M,N) = YLS(2,N) + DYLSX
              XLSX = MIN( YLS(M,N),XLSMX )
              XLS(M,N) = YLS(M,N)
              CALL SP_B( T(2,N),XLSX,PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
              CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
              IF( PVBX.GT.PGX .AND. M.EQ.2 ) THEN
                INDX = 17
                N_DB = N
                RLMSG = PVBX
                CHMSG = 'Initial Condition Transition: Vapor ' //
     &            'Pressure > Gas Pressure: PVBX = '
                CALL WRMSGS( INDX )
              ENDIF
              PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
              XLSSX = XLS(M,N)
              CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &          XMLSX,XMLWX )
              XLA(M,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
              XLS(M,N) = XLS(M,N) + (XLSSX-XLS(M,N))*(XLA(M,N)/XLAX)
              CALL DENS_L( T(2,N),RHOBX,XLA(M,N),RHOLX )
              GX(M,1) = XLS(M,N) - TMS(2,N)
  560       CONTINUE
            RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
            RPX(1) = -GX(2,1)
            CYLSX = RPX(1)/RX(1,1)
            YLS(2,N) = MAX( YLS(2,N)+CYLSX,0.D+0 )
            IF( ABS(CYLSX).GT.(1.D-9*XLSMX) ) GOTO 550
!
!---        Limit salt mass fraction  ---
!
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
            INDX = 0
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVBX = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
            IF( SGT(2,N).GT.EPSL ) THEN
              SG(2,N) = SGT(2,N)
              INDX = 1
            ENDIF
            CALL KSP_CO2( N,PG(2,N),PL(2,N),SG(2,N),SGT(2,N),SL(2,N),
     &        SDPF(N),SDPM(N),RKL(1,2,N),RKG(2,N),ASL(N),ASLMIN(2,N),
     &        ASGT(N),ESGTMX,SLRX,BTGL(2,N),INDX )
!
!---      Aqueous-salt molality (mol salt/kg water)  ---
!
          ELSEIF( ICBRN(N).EQ.4 ) THEN
            PGX = PG(2,N)+PATM
            PCX = PG(2,N)-PL(2,N)
            ISRX = 1
            CALL DENS_W( T(2,N),PGX,RHOLWX,RHOGWX,ISRX )
            YLS(2,N) = TMS(2,N)*WTMS/(1.D+0 + TMS(2,N)*WTMS)
            CALL SOL_LS( T(2,N),XLSMX )
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            IF( PVBX.GT.PGX ) THEN
              INDX = 17
              N_DB = N
              RLMSG = PVBX
              CHMSG = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVBX = '
              CALL WRMSGS( INDX )
            ENDIF
            PVA(2,N) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = XLS(2,N)
            CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &        XMLSX,XMLWX )
            XLA(2,N) = MIN( PVA(2,N)/PGAX,1.D+0 )*XLAX
            XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
            INDX = 0
            IF( SGT(2,N).GT.EPSL ) THEN
              SG(2,N) = SGT(2,N)
              INDX = 1
            ENDIF
            CALL KSP_CO2( N,PG(2,N),PL(2,N),SG(2,N),SGT(2,N),SL(2,N),
     &        SDPF(N),SDPM(N),RKL(1,2,N),RKG(2,N),ASL(N),ASLMIN(2,N),
     &        ASGT(N),ESGTMX,SLRX,BTGL(2,N),INDX )
          ENDIF
!
!---      Assign initial phase condition  ---
!
          CALL DENS_A( T(2,N),PGAX,RHOAX,I_VX )
          NPHAZ(2,N) = (I_VX+1)*2
        ENDIF
!
!---    Saturated conditions  ---
!
        ELSE
!
!---    Assign initial phase condition  ---
!
        NPHAZ(2,N) = 1
!
!---    Aqueous-salt concentration  ---
!
        IF( ICBRN(N).EQ.1 ) THEN
!
!---      Aqueous-CO2 concentration  ---
!
          IF( ICAIR(N).EQ.1 ) THEN
            RHOLSX = TMS(2,N)
            RHOLAX = PVA(2,N)
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
            CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &        'CO2 Aqu. Conc., Salt Aqu. Conc. = '
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &        'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_11( T(2,N),PX,PGX,RHOLSX,RHOLAX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
!
!---      Aqueous-CO2 relative saturation  ---
!
          ELSEIF( ICAIR(N).EQ.2 ) THEN
            RHOLSX = TMS(2,N)
            PHILAX = PVA(2,N)
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
            CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &        'CO2 Aqu. Rel. Sat., Salt Aqu. Conc. = '
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &        'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_12( T(2,N),PX,PGX,RHOLSX,PHILAX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
!
!---      Aqueous-CO2 mass fraction  ---
!
          ELSEIF( ICAIR(N).EQ.3 ) THEN
            RHOLSX = TMS(2,N)
            XLA(2,N) = PVA(2,N)
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
            CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &        'CO2 Aqu. Mass Frac., Salt Aqu. Conc. = '
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &        'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_13( T(2,N),PX,PGX,RHOLSX,YLS(2,N),XLS(2,N),
     &        XLA(2,N),CHMSGX )
          ENDIF
!
!---    Aqueous-salt relative saturation  ---
!
        ELSEIF( ICBRN(N).EQ.2 ) THEN
!
!---      Aqueous-CO2 concentration  ---
!
          IF( ICAIR(N).EQ.1 ) THEN
            PHILSX = TMS(2,N)
            RHOLAX = PVA(2,N)
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
            CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &        'Salt Aqu. Rel. Conc., CO2 Aqu. Conc.  = '
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &        'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_21( T(2,N),PX,PGX,PHILSX,RHOLAX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
!
!---      Aqueous-CO2 relative saturation  ---
!
          ELSEIF( ICAIR(N).EQ.2 ) THEN
            PHILSX = TMS(2,N)
            PHILAX = PVA(2,N)
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
            CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &        'Salt Aqu. Rel. Sat., CO2 Aqu. Rel. Sat. = '
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &        'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_22( T(2,N),PX,PGX,PHILSX,PHILAX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
!
!---      Aqueous-CO2 mass fraction  ---
!
          ELSEIF( ICAIR(N).EQ.3 ) THEN
            PHILSX = TMS(2,N)
            XLA(2,N) = PVA(2,N)
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
            CHMSGX(1) = 'null'
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_23( T(2,N),PX,PGX,PHILSX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
          ENDIF
!
!---    Aqueous-salt mass fraction  ---
!
        ELSEIF( ICBRN(N).EQ.3 .OR. ICBRN(N).EQ.0 ) THEN
!
!---      Aqueous-CO2 concentration  ---
!
          IF( ICAIR(N).EQ.1 ) THEN
            YLS(2,N) = TMS(2,N)
            RHOLAX = PVA(2,N)
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
            CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &        'CO2 Aqu. Conc., Salt Aqu. Mass Frac. = '
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &        'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_31( T(2,N),PX,PGX,RHOLAX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
!
!---      Aqueous-CO2 relative saturation  ---
!
          ELSEIF( ICAIR(N).EQ.2 ) THEN
            YLS(2,N) = TMS(2,N)
            PHILAX = PVA(2,N)
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
            CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &        'CO2 Aqu. Rel. Sat., Salt Aqu. Mass Frac. = '
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &        'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_32( T(2,N),PX,PGX,PHILAX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
!
!---      Aqueous-CO2 mass fraction  ---
!
          ELSEIF( ICAIR(N).EQ.3 ) THEN
            YLS(2,N) = TMS(2,N)
            XLA(2,N) = PVA(2,N)
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
            CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &          'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &            'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_33( T(2,N),PX,PGX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
          ENDIF
!
!---    Aqueous-salt molality (mol salt/kg water)  ---
!
        ELSEIF( ICBRN(N).EQ.4 ) THEN
!
!---      Aqueous-CO2 concentration  ---
!
          IF( ICAIR(N).EQ.1 ) THEN
            YLS(2,N) = TMS(2,N)
            RHOLAX = PVA(2,N)
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
            CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &          'CO2 Aqu. Conc., Salt Aqu. Molality = '
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &            'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_41( T(2,N),PX,PGX,RHOLAX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
!
!---      Aqueous-CO2 relative saturation  ---
!
          ELSEIF( ICAIR(N).EQ.2 ) THEN
            YLS(2,N) = TMS(2,N)
            PHILAX = PVA(2,N)
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
            CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &          'CO2 Aqu. Rel. Sat., Salt Aqu. Molality = '
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_42( T(2,N),PX,PGX,PHILAX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
!
!---      Aqueous-CO2 mass fraction  ---
!
          ELSEIF( ICAIR(N).EQ.3 ) THEN
            XLA(2,N) = PVA(2,N)
            YLS(2,N) = TMS(2,N)
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
            CHMSGX(1) = 'null'
            CHMSGX(2) = 'Initial Condition Transition: Vapor ' //
     &          'Pressure > Gas Pressure: PVB + PVA = '
            CALL FLSH_43( T(2,N),PX,PGX,YLS(2,N),
     &        XLS(2,N),XLA(2,N),CHMSGX )
          ENDIF
        ENDIF
        ENDIF
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
!#ifdef ecke
!        POR0(1,N) = POR0(1,N)
!        POR0(2,N) = POR0(2,N)
!#endif
 3000 CONTINUE
 3010 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CISC_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Compute initial solute concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE LEAK_WELL
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
      IF( IEQC.EQ.0 .AND. ISLC(40).EQ.0 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CISC_CO2'
      DO 140 NSL = 1,NSOLU
        DO 110 N = 1,NFLD+NWN_LW
          N_DB = N
!
!---      Skip inactive nodes  ---
!
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
!---  End of CISC_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HYDST_BC_CO2( BCX,PGX,PLX,TX,XLAX,XLSX,YLSX,ZPX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2e
!
!     Establish hydrostatic boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 March 2016
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
      CHARACTER*132 CHMSGX(2)
      REAL*8 BCX(*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HYDST_BC_CO2'
!
!---  Hydrostatic initial conditions, establish aqueous pressure,
!     temperature, and salt mass fraction at boundary surface  ---
!
      ZX = BCX(3)
      KC = MAX( 1,INT(ABS(ZX-ZPX)) )
      DISTZ = (ZX-ZPX)/REAL(KC)
      PLX = BCX(2) + PATM
      PGX = PLX
      DO 100 K = 1,KC
        ZX = ZX - 5.D-1*DISTZ
        TX = BCX(4) + (ZX-BCX(5))*BCX(6)
        YLSX = BCX(7) + (ZX-BCX(8))*BCX(9)
        XLAX = BCX(10) + (ZX-BCX(11))*BCX(12)
        ZX = ZX - 5.D-1*DISTZ
!
!---    Check for out-of-range salt mass fraction  ---
!
        CALL SOL_LS( TX,XLSMX )
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
        PX = MAX( PGX,PLX )
        CHMSGX(1) = 'Unconverged Hydrostatic Boundary Conditions: ' //
     &    'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
        CHMSGX(2) = 'Hydrostatic Boundary Condition Transition: ' //
     &    'Vapor Pressure > Gas Pressure: PVB + PVA = '
        CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
        CALL DENS_B( TX,PLX,XLSX,RHOBX )
        CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
        PLX = PLX + RHOLX*GRAVZ*DISTZ
        PGX = PLX
  100 CONTINUE
      PLX = PLX - PATM
      PGX = PLX
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HYDST_BC_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Compute primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE LEAK_WELL
      USE JACOB
      USE HYST
      USE GRID
      USE FDVS
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
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM_CO2'
      PETA = 1.D-1
      EPSLX = 1.D-4
!
!---  Phase options, compute phase condition   ---
!
      DO 300 N = 1,NFLD+NWN_LW
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) GOTO 300
        IZN = IZ(N)
        N_DB = N
!
!---    Minimum apparent aqueous saturation  ---
!
        ASLMINX = ASLMIN(2,N)
!
!---    Surface tension and saturation  ---
!
        CALL SOL_LS( T(2,N),XLSMX )
        XLSX = MIN( YLS(2,N),XLSMX )
        CALL SFT_L( T(2,N),XLSX,SFTLX )
        BTGLX = 1.D+0
        IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &    BTGLX = SCHR(16,IZN)/SFTLX
!
!---    Assign gas entry pressure and minimum gas saturation
!       for transition to unsaturated conditions---
!
        SGMNX = 1.D-3
        ENPR = 0.D+0
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 .OR.
     &    ISCHR(IZN).EQ.201 ) THEN
          SGMNX = 1.D+1**(-3.D+0+LOG10(1.D+0/SCHRV(1,N)))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 ) THEN
          ENPR = SCHRV(1,N)*RHORL*GRAV
          SGMNX = 1.D+1**(-3.D+0+LOG10(SCHRV(1,N)))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ELSEIF( ISCHR(IZN).EQ.3 ) THEN
          SCHRVX = MAX( SCHRV(1,N),SCHR(5,IZN) )
          SGMNX = 1.D+1**(-3.D+0+LOG10(1.D+0/SCHRVX))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          ENPR = MIN( SCHRV(1,N),SCHR(5,IZN) )*RHORL*GRAV
          SCHRVX = MIN( SCHRV(1,N),SCHR(5,IZN) )
          SGMNX = 1.D+1**(-3.D+0+LOG10(SCHRVX))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ENDIF
!
!---    Maximum effective trapped gas saturation  ---
!
        IF( ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.201 .OR. ISCHR(IZN).EQ.202 ) THEN
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
!---    Saturated system w/o trapped gas
!       Water mass - aqueous pressure
!       CO2 mass - aqueous-CO2 mass fraction
!       NaCl mass - total NaCl brine mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          XLS(2,N) = XLSX
          CALL SP_B( T(2,N),XLSX,PSBX )
          PVBX = PSBX
          PLX = PL(2,N)+PATM
          PX = MAX( PLX,PVBX )
          XLSSX = XLS(2,N)
          CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
          PGX = PGAX + PGWX
!
!---      Aqueous CO2 mass fraction exceeds solubility
!         limit transition to two-phase conditions  ---
!
          IF( XLA(2,N).GT.XLAX ) THEN
!
!---        Estimate gas saturation and transition to unsaturated
!           conditions if the gas saturation exceeds 0.001  ---
!
            CALL DENS_A( T(2,N),PGAX,RHOGAX,I_VX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL DENS_L( T(2,N),RHOBX,XLAX,RHOLX )
            SGX = (XLA(2,N)-XLAX)*RHOLX/RHOGAX
            IF( SGX.LT.SGMNX ) THEN
              PG(2,N) = PL(2,N)+(ENPR/BTGLX)-EPSLX
              NPHAZ(2,N) = 1
            ELSE
              SLX = 1.D+0-MIN( SGX,1.D-1 )
              SGTX = 0.D+0
              CALL SFT_L( T(2,N),XLSX,SFTLX )
              BTGL(2,N) = 1.D+0
              IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &          BTGL(2,N) = SCHR(16,IZN)/SFTLX
              CALL CAP_CO2( ASLMIN(2,N),BTGL(2,N),PCX,SLX,SGTX,N )
              PCX = MIN( PCX,(ENPR/BTGLX)+1.D+5 )
              PG(2,N) = PL(2,N)+PCX
              NPHAZ(2,N) = (I_VX+1)*2
            ENDIF
!
!---      Water-vapor pressure plus CO2 bubbling pressure less than
!         gas-entry pressure, remain as saturated condition  ---
!
          ELSE
            PG(2,N) = PL(2,N)+(ENPR/BTGLX)-EPSLX
            NPHAZ(2,N) = 1
          ENDIF
!
!---    Unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 .OR. NPHAZ(2,N).EQ.4 
     &    .OR. NPHAZ(2,N).EQ.6 ) THEN
          PGX = PG(2,N) + PATM
          PCX = MAX( PG(2,N)-PL(2,N),0.D+0 )
          XLS(2,N) = XLSX
          CALL SP_B( T(2,N),XLSX,PSBX )
          PX = MAX( PGX,PSBX )
          CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
          CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVW(2,N),XLSX,IZN )
          PVAX = MAX( PX-PVW(2,N),0.D+0 )
          CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVW(2,N),
     &      XGAX,XGWX,XLAX,XLS(2,N),XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          INDX = 0
          CALL KSP_CO2( N,PG(2,N),PL(2,N),SGX,SGTX,SLX,SLFX,SLMX,
     &      RKLX,RKGX,ASLX,ASLMINX,ESGTX,ESGTMX,
     &      SLRX,BTGL(2,N),INDX )
!         
!---      Apparent aqueous saturation (i.e., aqueous saturation +
!         trapped gas saturation = 1.0), transition to saturation
!         conditions  ---
!
          IF( ABS(1.D+0-ASLX).LT.EPSL ) THEN
!
!---        Trapped gas exists transition to saturated condition
!           w/ entrapped gas  ---
!
            IF( ESGTMX.GT.EPSL .AND. SGTX.GT.EPSL ) THEN
              CALL CAP_CO2( ASLMIN(2,N),BTGL(2,N),PCX,SLX,SGTX,N )
              PG(2,N) = PL(2,N) + PCX
              SG(2,N) = SGX
              SGT(2,N) = SGX
              CALL DENS_A( T(2,N),PGAX,RHOAX,I_VX )
              IF( I_VX.EQ.0 ) THEN
                IF( NPHAZ(2,N).EQ.4 ) THEN
                  INDX = 28
                  N_DB = N
                  CHMSG = 'Transition from Subcritical Liquid to Gas' //
     &              ', Gas Pressure = '
                  RLMSG = PG(2,N) + PATM
                  CALL WRMSGS( INDX )
                  ICNV = 4
                ENDIF
                NPHAZ(2,N) = 3
              ELSEIF( I_VX.EQ.1 ) THEN
                IF( NPHAZ(2,N).EQ.2 ) THEN
                  INDX = 28
                  N_DB = N
                  CHMSG = 'Transition from Subcritical Gas to Liquid' //
     &              ', Gas Pressure = '
                  RLMSG = PG(2,N) + PATM
                  CALL WRMSGS( INDX )
                  ICNV = 4
                ENDIF
                NPHAZ(2,N) = 5
              ELSE
                NPHAZ(2,N) = 7
              ENDIF
!
!---        No trapped gas transition to saturated condition
!           w/o entrapped gas  ---
!
            ELSE
              PX = PL(2,N)+PATM
              XLS(2,N) = XLSX
              CALL SP_B( T(2,N),XLSX,PSBX )
              PVBX = PSBX
              CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &          XGAX,XGWX,XLA(2,N),XLS(2,N),XLWX,XMGAX,XMGWX,XMLAX,
     &          XMLSX,XMLWX )
              PG(2,N) = PL(2,N)+(ENPR/BTGLX)-EPSLX
              NPHAZ(2,N) = 1
            ENDIF
!
!---      No aqueous phase, transition
!         to fully unsaturated condition  ---
!
          ELSEIF( SLX.LT.EPSL .AND. (1.D+0-SL(1,N)).GT.EPSL ) THEN
            PGX = PG(2,N) + PATM
            CALL SOL_LS( T(2,N),XLSMX )
            IF( TMS(2,N).GT.EPSL ) THEN
              XLSX = XLSMX
            ELSE
              XLSX = 0.D+0
            ENDIF
            XLS(2,N) = XLSX
            CALL SFT_L( T(2,N),XLSX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
            SL(2,N) = 0.D+0
            SGT(2,N) = 0.D+0
            CALL CAP_CO2( ASLMIN(2,N),BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
            PL(2,N) = PG(2,N) - PCX
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVW(2,N),XLSX,IZN )
            PVA(2,N) = PGX - PVW(2,N)
!
!---        Check for transitions across the saturation line  ---
!
            CALL DENS_A( T(2,N),PVA(2,N),RHOAX,I_VX )
            IF( I_VX.EQ.0 ) THEN
              IF( NPHAZ(2,N).EQ.4 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Liquid to Gas' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 8
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.2 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Gas to Liquid' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 9
            ELSE
              NPHAZ(2,N) = 10
            ENDIF
!
!---      Gas pressure remains above gas-entry pressure, remain
!         as unsaturated condition  ---
!
          ELSE
!
!---        Check for transition to fully unsaturated conditions
!           and adjust capillary pressure  ---
!
            IF( SLX.LT.EPSL ) THEN
              SLX = 1.D-3
              SGTX = 0.D+0
              CALL CAP_CO2( ASLMIN(2,N),BTGL(2,N),PCX,SLX,SGTX,N )
              PL(2,N) = PG(2,N)-PCX
            ENDIF
!
!---        Check for transitions across the saturation line  ---
!
            CALL DENS_A( T(2,N),PGAX,RHOAX,I_VX )
            IF( I_VX.EQ.0 ) THEN
              IF( NPHAZ(2,N).EQ.4 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Liquid to Gas' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 2
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.2 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Gas to Liquid' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 4
            ELSE
              NPHAZ(2,N) = 6
            ENDIF
          ENDIF
!
!---    Saturated system w/ entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - trapped gas saturation
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.5
     &    .OR. NPHAZ(2,N).EQ.7 ) THEN
          CALL SFT_L( T(2,N),XLSX,SFTLX )
          BTGL(2,N) = 1.D+0
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &      BTGL(2,N) = SCHR(16,IZN)/SFTLX
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
          IF( NPHAZ(1,N).EQ.3 .OR. NPHAZ(1,N).EQ.5
     &      .OR. NPHAZ(1,N).EQ.7 ) THEN
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
!---      Aqueous-gas equilibrium  ---
!
          CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLA(2,N),XLS(2,N),XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
!
!---      Trapped gas saturation disappears, transition to
!         saturated condition w/o entrapped gas  ---
!
          IF( SG(2,N).LT.EPSL ) THEN
            SGT(2,N) = 0.D+0
            SG(2,N) = 0.D+0
            PG(2,N) = PL(2,N)+(ENPR/BTGLX)-EPSLX
            ASLMIN(2,N) = 1.D+0
            NPHAZ(2,N) = 1
!
!---      Trapped gas saturation exceeds previous trapped gas
!         saturation or maximum trapped gas saturation, transition 
!         to free gas w/ trapped gas phase condition  ---
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
            CALL CAP_CO2( ASLMIN(2,N),BTGL(2,N),PCX,SLX,SGT(2,N),N )
            PCX = MIN( PCX,(ENPR/BTGLX)+1.D+5 )
            PG(2,N) = PL(2,N) + PCX
            CALL DENS_A( T(2,N),PGAX,RHOAX,I_VX )
            IF( I_VX.EQ.0 ) THEN
              IF( NPHAZ(2,N).EQ.5 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Liquid to Gas' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 2
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.3 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Gas to Liquid' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 4
            ELSE
              NPHAZ(2,N) = 6
            ENDIF
!
!---      Trapped gas saturation neither disappears nor
!         exceeds previous trapped gas saturation,
!         remain as saturated w/ entrapped gas condition  ---
!
          ELSE
            IF( ISCHR(IZN).EQ.202 ) THEN
              PG(2,N) = PL(2,N) + (SCHR(5,IZN)/BTGLX)-EPSLX
            ELSE
              PG(2,N) = PL(2,N) + (ENPR/BTGLX)-EPSLX
            ENDIF
            CALL DENS_A( T(2,N),PGAX,RHOAX,I_VX )
            IF( I_VX.EQ.0 ) THEN
              IF( NPHAZ(2,N).EQ.5 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Liquid to Gas' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 3
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.3 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Gas to Liquid' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 5
            ELSE
              NPHAZ(2,N) = 7
            ENDIF
          ENDIF
!
!---    Fully unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous saturation
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.8 .OR. NPHAZ(2,N).EQ.9 .OR.
     &    NPHAZ(2,N).EQ.10 ) THEN
          PGX = PG(2,N) + PATM
          CALL SOL_LS( T(2,N),XLSMX )
          IF( TMS(2,N).GT.EPSL ) THEN
            XLSX = XLSMX
          ELSE
            XLSX = 0.D+0
          ENDIF
          XLS(2,N) = XLSX
          CALL SFT_L( T(2,N),XLSX,SFTLX )
          BTGL(2,N) = 1.D+0
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &      BTGL(2,N) = SCHR(16,IZN)/SFTLX
          SL(2,N) = 0.D+0
          SGT(2,N) = 0.D+0
          CALL CAP_CO2( ASLMIN(2,N),BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
          PL(2,N) = PG(2,N) - PCX
          CALL SP_B( T(2,N),XLSX,PSBX )
          PX = MAX( PGX,PSBX )
          PVA(2,N) = PGX - PVW(2,N)
          CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVW(2,N),
     &      XGAX,XGWX,XLAX,XLS(2,N),XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
          CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
!
!---      Aqueous phase appears, transition to
!         unsaturated conditions  ---
!
          IF( PVW(2,N).GT.PVBX ) THEN
            SL(2,N) = 1.D-4
            YLS(2,N) = TMS(2,N)/(RHOL(2,N)*SL(2,N)*PORD(2,N))
            CALL SFT_L( T(2,N),XLSMX,SFTLX )
            BTGL(2,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(2,N) = SCHR(16,IZN)/SFTLX
            SGTX = 0.D+0
            CALL CAP_CO2( ASLMIN(2,N),BTGL(2,N),PCX,SL(2,N),SGTX,N )
            PL(2,N) = PG(2,N) - PCX
!
!---        Check for transitions across the saturation line  ---
!
            CALL DENS_A( T(2,N),PGAX,RHOAX,I_VX )
            IF( I_VX.EQ.0 ) THEN
              IF( NPHAZ(2,N).EQ.9 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Liquid to Gas' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 2
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.8 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Gas to Liquid' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 4
            ELSE
              NPHAZ(2,N) = 6
            ENDIF
!
!---      No aqueous phase, no transition from
!         fully unsaturated condition  ---
!
          ELSE
!
!---        Check for transitions across the saturation line  ---
!
            CALL DENS_A( T(2,N),PGAX,RHOAX,I_VX )
            IF( I_VX.EQ.0 ) THEN
              IF( NPHAZ(2,N).EQ.9 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Liquid to Gas' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 8
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.8 ) THEN
                INDX = 28
                N_DB = N
                CHMSG = 'Transition from Subcritical Gas to Liquid' //
     &            ', Gas Pressure = '
                RLMSG = PG(2,N) + PATM
                CALL WRMSGS( INDX )
                ICNV = 4
              ENDIF
              NPHAZ(2,N) = 9
            ELSE
              NPHAZ(2,N) = 10
            ENDIF
          ENDIF
       ENDIF
!
!---  Compute increments  ---
!
!
!---    Assign aqueous-salt mass fraction increments  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          CALL SOL_LS( T(2,N),XLSMX )
          XLSX = MIN( YLS(2,N),XLSMX )
          DNR(IEQS,N) = 1.D-5*XLSMX
        ENDIF
!
!---    Surface tension and saturation  ---
!
        CALL SFT_L( T(2,N),XLSX,SFTLX )
        BTGLX = 1.D+0
        IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &    BTGLX = SCHR(16,IZN)/SFTLX
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
!---    Saturated system w/o trapped gas
!       Water mass - aqueous pressure
!       CO2 mass - aqueous-CO2 mass fraction
!       NaCl mass - total NaCl brine mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          DNR(IEQW,N) = MAX( 1.D-1,1.D-7*(PL(2,N)+PATM) )
          PX = PL(2,N)+PATM
          CALL SP_B( T(2,N),XLSX,PSBX )
          PVBX = PSBX
          CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAMX,XLSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          IF( XLA(2,N).GT.(1.D-2*XLAMX) ) THEN
            DNR(IEQA,N) = SIGN( 1.D-4*XLAMX,5.D-1*XLAMX-XLA(2,N) )
          ELSE
            DNR(IEQA,N) = SIGN( 1.D-3*XLAMX,5.D-1*XLAMX-XLA(2,N) )
          ENDIF
!
!---    Unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 .OR. NPHAZ(2,N).EQ.4 
     &    .OR. NPHAZ(2,N).EQ.6 ) THEN
!          IF( 1.D-5*ABS(PG(2,N)-PL(2,N)).GT.1.D-1 ) THEN
!            DNR(IEQW,N) = 1.D-5*ABS(PG(2,N)-PL(2,N))
!          ELSE
!            DNR(IEQW,N) = -1.D-1
!          ENDIF
          DNRX = MAX( 1.D-2,1.D-6*ABS(PG(2,N)-PL(2,N)) )
          DNR(IEQW,N) = SIGN( DNRX,(5.D-1-SL(2,N)) )
          DNR(IEQA,N) = -DNR(IEQW,N)
!
!---    Saturated system w/ entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - trapped gas saturation
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.5
     &    .OR. NPHAZ(2,N).EQ.7 ) THEN
          DNR(IEQW,N) = MAX( 1.D-1,1.D-6*(PL(2,N)+PATM) )
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
          DNR(IEQA,N) = SIGN( 1.D-7,5.D-1*SGTMX-SG(2,N) )
!
!---    Fully unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous saturation
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.8 .OR. NPHAZ(2,N).EQ.9
     &    .OR. NPHAZ(2,N).EQ.10 ) THEN
          DNR(IEQW,N) = 1.D-3
          DNR(IEQA,N) = MAX( 1.D-2,1.D-7*ABS(PG(2,N)-PL(2,N)) )
          IF( ISLC(32).EQ.0 ) DNR(IEQS,N) = 1.D-6
        ENDIF
!
!--- Increment the primary variables  ---
!
        DO 200 M = 3,ISVC+2
          T(M,N) = T(2,N)
          PL(M,N) = PL(2,N)
          PG(M,N) = PG(2,N)
          PVW(M,N) = PVW(2,N)
          XLA(M,N) = XLA(2,N)
          SG(M,N) = SG(2,N)
          SL(M,N) = SL(2,N)
          YLS(M,N) = YLS(2,N)
          TMS(M,N) = TMS(2,N)
!
!---      Saturated system w/o trapped gas
!         Water mass - aqueous pressure
!         CO2 mass - aqueous-CO2 mass fraction
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ(2,N).EQ.1 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
              PG(M,N) = PL(M,N) + (ENPR/BTGLX)-EPSLX
            ELSEIF( M.EQ.IEQA+2 ) THEN
              XLA(M,N) = XLA(M,N) + DNR(IEQA,N)
!
!---        Isobrine option  ---
!
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Unsaturated w/ or w/o trapped gas  ---
!
          ELSEIF( NPHAZ(2,N).EQ.2 .OR. NPHAZ(2,N).EQ.4 
     &      .OR. NPHAZ(2,N).EQ.6 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQA+2 ) THEN
              PG(M,N) = PG(M,N) + DNR(IEQA,N)
!
!---        Isobrine option  ---
!
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Saturated system w/ entrapped gas
!         Water mass - aqueous pressure
!         CO2 mass - trapped gas saturation
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.5
     &      .OR. NPHAZ(2,N).EQ.7 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
              IF( ISCHR(IZN).EQ.202 ) THEN
                PG(M,N) = PL(M,N) + (SCHR(5,IZN)/BTGLX)-EPSLX
              ELSE
                PG(M,N) = PL(M,N) + (ENPR/BTGLX)-EPSLX
              ENDIF
            ELSEIF( M.EQ.IEQA+2 ) THEN
              SG(M,N) = SG(M,N) + DNR(IEQA,N)
!
!---        Isobrine option  ---
!
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---    Fully unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous saturation
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.8 .OR. NPHAZ(2,N).EQ.9
     &    .OR. NPHAZ(2,N).EQ.10 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              PVW(M,N) = PVW(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQA+2 ) THEN
              PG(M,N) = PG(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              TMS(M,N) = TMS(M,N) + DNR(IEQS,N)
            ENDIF
          ENDIF
  200   CONTINUE
  300 CONTINUE
!
!---  Illegal phase transition, stop simulation  ---
!
      IF( ICNV.EQ.4 ) THEN
        TM = TM - DT
        NSTEP = NSTEP - 1
        DT = DTO
        DTI = 1.D+0/DT
        DO N = 1,NFLD+NWN_LW
          PL(2,N) = PL(1,N)
          PG(2,N) = PG(1,N)
          PVW(2,N) = PVW(1,N)
          XLA(2,N) = XLA(1,N)
          SG(2,N) = SG(1,N)
          SGT(2,N) = SGT(1,N)
          ASLMIN(2,N) = ASLMIN(1,N)
          YLS(2,N) = YLS(1,N)
          TMS(2,N) = TMS(1,N)
          NPHAZ(2,N) = NPHAZ(1,N)
        ENDDO
!
!---    Coupled-well pressure  ---
!
        IF( L_CW.EQ.1 ) THEN
          DO NCW = 1,N_CW
            P_CW(2,NCW) = P_CW(1,NCW)
          ENDDO
        ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KSP_CO2( N,PGX,PLX,SGX,SGTX,SLX,SLFX,SLMX,RKLX,RKGX,
     &  ASLX,ASLMINX,ESGTX,ESGTMX,SLRX,BTGLX,INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Aqueous saturation, gas saturation, aqueous relative permeability,
!     and gas relative permeability.
!
!     INDX = 0 : Trapped-gas saturation computed.
!     INDX = 1 : Trapped-gas saturation given.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
      SUB_LOG(ISUB_LOG) = '/KSP_CO2'
!
!---  Aqueous and gas saturation  ---
!
      CALL SP_CO2( ASL_F,ASL_M,ASLX,ASLMINX,BTGLX,ESLX,ESGTX,ESGTMX,
     &  PGX,PLX,SGX,SGTX,SLX,SLFX,SLMX,SLRX,INDX,N )
!
!---  Aqueous relative permeability  ---
!
      CALL RKLS_CO2( ASL_F,ASL_M,BTGLX,ESLX,PGX,PLX,RKLX,
     &  SLX,N )
!
!---  Gas relative permeability  ---
!
      CALL RKGS_CO2( ASL_F,ASL_M,ASLX,ASLMINX,BTGLX,ESGTX,PGX,PLX,RKGX,
     &  SGX,SGTX,SLX,N )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KSP_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Load the current time step values into the old time step
!     variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE LEAK_WELL
      USE HYST
      USE GRID
      USE FDVS
      USE FDVP
      USE FDVG
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
      SUB_LOG(ISUB_LOG) = '/LDO_CO2'
!
!---  Assign old time step values  ---
!
      DO 100 N = 1,NFLD+NWN_LW
!
!---    Saturated system w/o trapped gas
!       Water mass - aqueous pressure
!       CO2 mass - aqueous-CO2 mass fraction
!       NaCl mass - total NaCl brine mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          ASLMIN(1,N) = 1.D+0
          ASLMIN(2,N) = 1.D+0
!
!---    Unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 .OR. NPHAZ(2,N).EQ.4 
     &    .OR. NPHAZ(2,N).EQ.6 ) THEN
          ASLMIN(1,N) = MIN( ASL(N),ASLMIN(1,N),ASLMIN(2,N) )
          ASLMIN(1,N) = MAX( ASLMIN(1,N),0.D+0 )
          ASLMIN(2,N) = ASLMIN(1,N)
!
!---    Saturated system w/ entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - trapped gas saturation
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.5
     &    .OR. NPHAZ(2,N).EQ.7 ) THEN
          ASLMIN(1,N) = ASLMIN(2,N)
!
!---    Fully unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous saturation
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.8 .OR. NPHAZ(2,N).EQ.9
     &    .OR. NPHAZ(2,N).EQ.10 ) THEN
          ASLMIN(1,N) = 0.D+0
          ASLMIN(2,N) = 0.D+0
        ENDIF
        BTGL(1,N) = BTGL(2,N)
        DFGW(1,N) = DFGW(2,N)
        DFLA(1,N) = DFLA(2,N)
        NPHAZ(1,N) = NPHAZ(2,N)
        PG(1,N) = PG(2,N)
        PSO(1,N) = PSO(2,N)
        PL(1,N) = PL(2,N)
        PVW(1,N) = PVW(2,N)
        PORD(1,N) = PORD(2,N)
        PORT(1,N) = PORT(2,N)
        RHOG(1,N) = RHOG(2,N)
        RHOL(1,N) = RHOL(2,N)
        RHOSP(1,N) = RHOSP(2,N)
        RKG(1,N) = RKG(2,N)
        RKL(1,1,N) = RKL(1,2,N)
        RKL(2,1,N) = RKL(2,2,N)
        RKL(3,1,N) = RKL(3,2,N)
        TMS(1,N) = TMS(2,N)
        SG(1,N) = SG(2,N)
        SGT(1,N) = SGT(2,N)
        SL(1,N) = SL(2,N)
        SS(1,N) = SS(2,N)
        T(1,N) = T(2,N)
        VISG(1,N) = VISG(2,N)
        VISL(1,N) = VISL(2,N)
        XGA(1,N) = XGA(2,N)
        XGW(1,N) = XGW(2,N)
        XLA(1,N) = XLA(2,N)
        XLS(1,N) = XLS(2,N)
        XLW(1,N) = XLW(2,N)
        XMGA(1,N) = XMGA(2,N)
        XMGW(1,N) = XMGW(2,N)
        YLS(1,N) = YLS(2,N)
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
      DO 110 NB = 1,NBC
        PLB(1,NB) = PLB(2,NB)
        PGB(1,NB) = PGB(2,NB)
  110 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Compute hydrologic, thermodynamic and physical properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE LEAK_WELL
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROP_CO2'
!
!---  Loop over all nodes  ---
!
      DO 200 N = 1,NFLD+NWN_LW
        IZN = IZ(N)
!#ifdef ecke
!        POR0(1,N) = POR0(1,N)
!        POR0(2,N) = POR0(2,N)
!#endif
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) GOTO 200
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
        DO 100 M = 2,ISVC+2
          PLX = PL(M,N)+PATM
          PGX = PG(M,N)+PATM
!
!---      Saturated system w/o trapped gas
!         Water mass - aqueous pressure
!         CO2 mass - aqueous-CO2 mass fraction
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ(2,N).EQ.1 ) THEN
            PX = PLX
            PCX = 0.D+0
            SL(M,N) = 1.D+0
            CALL SOL_LS( T(M,N),XLSMX )
            XLSX = MIN(YLS(M,N),XLSMX)
            XLS(M,N) = XLSX
            CALL SP_B( T(M,N),XLSX,PSBX )
            CALL DENS_B( T(M,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(M,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            ISRX = 2
            CALL DENS_W( T(M,N),PVBX,RHOLWX,RHOGWX,ISRX )
            PVA(M,N) = MAX( PX-PVBX,0.D+0 )
            PVW(M,N) = PVBX
            XLSSX = XLS(M,N)
            CALL EQUIL( T(M,N),PX,PGAX,PGWX,PSBX,PVW(M,N),
     &        XGA(M,N),XGW(M,N),XLAX,XLSSX,XLWX,XMGA(M,N),
     &        XMGW(M,N),XMLAX,XMLSX,XMLWX )
            XLS(M,N) = XLS(M,N) + (XLSSX-XLS(M,N))*(XLA(M,N)/XLAX)
            XLW(M,N) = 1.D+0-XLA(M,N)-XLS(M,N)
            WTMLX = 1.D+0/(XLA(M,N)/WTMA+XLS(M,N)/WTMS+XLW(M,N)/WTMW)
            XMLA(M,N) = XLA(M,N)*WTMLX/WTMA
            XMLS(M,N) = XLS(M,N)*WTMLX/WTMS
            XMLW(M,N) = XLW(M,N)*WTMLX/WTMW
!
!---      Unsaturated system w/ or w/o entrapped gas
!         Water mass - aqueous pressure
!         CO2 mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.2 .OR. NPHAZ(2,N).EQ.4 
     &      .OR. NPHAZ(2,N).EQ.6 ) THEN
            PCX = PG(M,N)-PL(M,N)
            CALL SOL_LS( T(M,N),XLSMX )
            XLSX = MIN(YLS(M,N),XLSMX)
            XLS(M,N) = XLSX
            CALL SP_B( T(M,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(M,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(M,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            ISRX = 2
            CALL DENS_W( T(M,N),PVBX,RHOLWX,RHOGWX,ISRX )
            PVA(M,N) = MAX( PGX-PVBX,0.D+0 )
            PVW(M,N) = PVBX
            CALL EQUIL( T(M,N),PX,PGAX,PGWX,PSBX,PVW(M,N),
     &        XGA(M,N),XGW(M,N),XLA(M,N),XLS(M,N),XLW(M,N),XMGA(M,N),
     &        XMGW(M,N),XMLA(M,N),XMLS(M,N),XMLW(M,N) )
!
!---      Saturated system w/ entrapped gas
!         Water mass - aqueous pressure
!         CO2 mass - trapped gas saturation
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.5
     &      .OR. NPHAZ(2,N).EQ.7 ) THEN
            PX = PLX
            PCX = 0.D+0
            CALL SOL_LS( T(M,N),XLSMX )
            XLSX = MIN(YLS(M,N),XLSMX)
            XLS(M,N) = XLSX
            CALL SP_B( T(M,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(M,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(M,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
            ISRX = 2
            CALL DENS_W( T(M,N),PVBX,RHOLWX,RHOGWX,ISRX )
            PVA(M,N) = MAX( PX-PVBX,0.D+0 )
            PVW(M,N) = PVBX
            CALL EQUIL( T(M,N),PX,PGAX,PGWX,PSBX,PVBX,
     &        XGA(M,N),XGW(M,N),XLA(M,N),XLS(M,N),XLW(M,N),XMGA(M,N),
     &        XMGW(M,N),XMLA(M,N),XMLS(M,N),XMLW(M,N) )
!
!---      Fully unsaturated system w/ or w/o entrapped gas
!         Water mass - aqueous saturation
!         CO2 mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.8 .OR. NPHAZ(2,N).EQ.9
     &      .OR. NPHAZ(2,N).EQ.10 ) THEN
            PGX = PG(M,N) + PATM
            CALL SOL_LS( T(M,N),XLSMX )
            IF( TMS(M,N).GT.EPSL ) THEN
              XLSX = XLSMX
            ELSE
              XLSX = 0.D+0
            ENDIF
            XLS(M,N) = XLSX
            CALL SFT_L( T(M,N),XLSX,SFTLX )
            BTGL(M,N) = 1.D+0
            IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &        BTGL(M,N) = SCHR(16,IZN)/SFTLX
            SL(M,N) = 0.D+0
            SGT(M,N) = 0.D+0
            CALL CAP_CO2( ASLMIN(2,N),BTGL(M,N),PCX,SL(M,N),SGT(M,N),N )
            PL(M,N) = PG(M,N) - PCX
            CALL SP_B( T(M,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            PVBX = PVW(M,N)
            ISRX = 2
            CALL DENS_W( T(M,N),PVBX,RHOLWX,RHOGWX,ISRX )
            PVA(M,N) = MAX( PGX-PVBX,0.D+0 )
            CALL EQUIL( T(M,N),PX,PGAX,PGWX,PSBX,PVW(M,N),
     &        XGA(M,N),XGW(M,N),XLA(M,N),XLS(M,N),XLW(M,N),XMGA(M,N),
     &        XMGW(M,N),XMLA(M,N),XMLS(M,N),XMLW(M,N) )
          ENDIF
!
!---      Porous-media porosity  ---
!
          CALL PORSTY_CO2E( N,PX,PCMP(N),PORD(M,N),PORT(M,N) )
          PORD(M,N) = MAX( PORD(M,N),EPSL )
          PORT(M,N) = MAX( PORT(M,N),PORD(M,N) )
!
!---      Surface tension, saturation, relative permeability  ---
!
          CALL SFT_L( T(M,N),XLSX,SFTLX )
          BTGL(M,N) = 1.D+0
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &      BTGL(M,N) = SCHR(16,IZN)/SFTLX
          INDX = 0
          IF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.5
     &      .OR. NPHAZ(2,N).EQ.7 ) THEN
            SGT(M,N) = SG(M,N)
            INDX = 1
          ENDIF
          CALL KSP_CO2( N,PG(M,N),PL(M,N),SG(M,N),SGT(M,N),SL(M,N),
     &      SDPF(N),SDPM(N),RKL(1,M,N),RKG(M,N),ASLX,ASLMINX,
     &      ESGTX,ESGTMX,SLRX,BTGL(M,N),INDX )
          IF( M.EQ.2 ) THEN
            ASL(N) = ASLX
            ASGT(N) = ESGTX
          ENDIF
!
!---      Gas density and component fractions  ---
!
          CALL DENS_A( T(M,N),PVA(M,N),RHOGAX,I_VX )
          IF( ISLC(9).EQ.1 ) THEN
            RHOG(M,N) = RHOGI
          ELSE
            RHOG(M,N) = XGA(M,N)*RHOGAX + XGW(M,N)*RHOGWX
          ENDIF
          WTMGX = XMGA(M,N)*WTMA + XMGW(M,N)*WTMW
          RHOMG(M,N) = RHOG(M,N)/WTMGX
!
!---      Gas viscosity  ---
!
          CALL VISC_A( T(M,N),RHOGAX,VISGAX )
          CALL VISC_W( T(M,N),PVBX,RHOGWX,VISGWX )
          CALL VISC_G( VISGAX,VISGWX,XMGA(M,N),XMGW(M,N),VISG(M,N) )
!
!---      Water-vapor diffusion coefficient  ---
!
          CALL DIFC_GW( T(M,N),PX,DFGW(M,N) )
!
!---      Aqueous component fractions and density  ---
!
          CALL DENS_L( T(M,N),RHOBX,XLA(M,N),RHOL(M,N) )
          WTMLX = XMLA(M,N)*WTMA + XMLS(M,N)*WTMS + XMLW(M,N)*WTMW
          RHOML(M,N) = RHOL(M,N)/WTMLX
!
!---      Aqueous viscosity  ---
!
          ISRX = 1
          CALL DENS_W( T(M,N),PX,RHOLWX,RHOX,ISRX )
          CALL VISC_W( T(M,N),PX,RHOLWX,VISLWX )
          CALL VISC_B( T(M,N),XLS(M,N),VISLWX,VISBX )
          CALL VISC_L( XMLA(M,N),VISBX,VISGAX,VISL(M,N) )
!
!---      Aqueous-CO2 and -NaCl diffusion coefficient  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLA(M,N) = DFLAC
            DFLS(M,N) = DFLSC
          ELSEIF( ISLC(4).EQ.2 ) THEN
!            CALL DIFC_LA( VISL(M,N),VISGAX,DFLA(M,N) )
            CALL DIFC_LA( T(M,N),XLS(M,N),DFLA(M,N) )
            CALL DIFC_LS( T(M,N),XLS(M,N),VISL(M,N),DFLS(M,N) )
          ENDIF
!
!---      Precipitated NaCl density and saturation  ---
!
          CALL DENS_S( T(M,N),PX,RHOSP(M,N) )
!
!---      Fully unsaturated system w/ or w/o entrapped gas
!         Water mass - aqueous saturation
!         CO2 mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ(2,N).EQ.8 .OR. NPHAZ(2,N).EQ.9
     &      .OR. NPHAZ(2,N).EQ.10 ) THEN
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
#ifdef ecke
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
#else
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
!---        Permeability reduction factor with salt precipitation  ---
!
            PERMRF(M,N) = ((PORD(M,N)-PERM(5,IZN))/
     &        (PERM(4,IZN)-PERM(5,IZN)))**PERM(6,IZN)
!
!---      Permeability reduction factor with Verma and Preuss model  ---
!
          ELSEIF( IPRF(IZN).EQ.1 ) THEN
            CALL PERM_R( SS(M,N),PERMRF(M,N),PORD(M,N),IZN )
          ENDIF
#endif
!
!---      Aqueous and gas tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL(M,N),SG(M,N),ZERO,
     &      PORD(M,N),TORL(M,N),TORG(M,N),TORNX )
!
!---      Gas pressure for Darcy velocity calculation  ---
!
          PSO(M,N) = PG(M,N)
  100   CONTINUE
!
!---    Saturated system w/o trapped gas
!       Water mass - aqueous pressure
!       CO2 mass - aqueous-CO2 mass fraction
!       NaCl mass - total NaCl brine mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          MX = IEQW + 2
          TORL(1,N) = -RHOL(2,N)*
     &      (1.D+0/RHOL(MX,N)-1.D+0/RHOL(2,N))/(PL(MX,N)-PL(2,N)+SMALL)
          TORG(1,N) = -RHOG(2,N)*
     &      (1.D+0/RHOG(MX,N)-1.D+0/RHOG(2,N))/(PL(MX,N)-PL(2,N)+SMALL)
!
!---    Unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 .OR. NPHAZ(2,N).EQ.4 
     &    .OR. NPHAZ(2,N).EQ.6 ) THEN
          MX = IEQW + 2
          TORL(1,N) = -RHOL(2,N)*
     &      (1.D+0/RHOL(MX,N)-1.D+0/RHOL(2,N))/(PL(MX,N)-PL(2,N)+SMALL)
          MX = IEQA + 2
          TORG(1,N) = -RHOG(2,N)*
     &      (1.D+0/RHOG(MX,N)-1.D+0/RHOG(2,N))/(PG(MX,N)-PG(2,N)+SMALL)
!
!---    Saturated system w/ entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - trapped gas saturation
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.5
     &    .OR. NPHAZ(2,N).EQ.7 ) THEN
          MX = IEQW + 2
          TORL(1,N) = -RHOL(2,N)*
     &      (1.D+0/RHOL(MX,N)-1.D+0/RHOL(2,N))/(PL(MX,N)-PL(2,N)+SMALL)
          TORG(1,N) = -RHOG(2,N)*
     &      (1.D+0/RHOG(MX,N)-1.D+0/RHOG(2,N))/(PL(MX,N)-PL(2,N)+SMALL)
!
!---    Fully unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous saturation
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.8 .OR. NPHAZ(2,N).EQ.9
     &    .OR. NPHAZ(2,N).EQ.10 ) THEN
          MX = IEQA + 2
          TORL(1,N) = -RHOL(2,N)*
     &      (1.D+0/RHOL(MX,N)-1.D+0/RHOL(2,N))/(PG(MX,N)-PG(2,N)+SMALL)
          TORG(1,N) = -RHOG(2,N)*
     &      (1.D+0/RHOG(MX,N)-1.D+0/RHOG(2,N))/(PG(MX,N)-PG(2,N)+SMALL)
        ENDIF
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDAT_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Reads the CO2 transport card dispersivities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 April 2012
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
      SUB_LOG(ISUB_LOG) = '/RDAT_CO2'
!
!---  Write card information to output file  ---
!
      CARD = 'CO2 Transport Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      IDSPS = 0
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
        VARB = 'CO2 Gas Longitudinal Dispersivity: '
        IF( IJK.GT.0 ) THEN
          UNTS = 'm'
          IUNM = 1
          CALL RDIJK( ISTART,IJK,CHDUM,UNTS,DPLGA )
          IDSPS = 1
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DPLGA(IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &      DPLGA(IROCK)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,DPLGA(IROCK),INDX)
          IF( DPLGA(IROCK).GE.SMALL ) IDSPS = 1
        ENDIF
!
!---  Transverse dispersivity  ---
!
        VARB = 'CO2 Gas Transverse Dispersivity: '
        IF( IJK.GT.0 ) THEN
          UNTS = 'm'
          IUNM = 1
          CALL RDIJK( ISTART,IJK,CHDUM,UNTS,DPTGA )
          IDSPS = 1
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DPTGA(IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &      DPTGA(IROCK)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,DPTGA(IROCK),INDX)
          IF( DPTGA(IROCK).GE.SMALL ) IDSPS = 1
        ENDIF
!
!---  Longitudinal dispersivity  ---
!
        VARB = 'CO2 Aqueous Longitudinal Dispersivity: '
        IF( IJK.GT.0 ) THEN
          UNTS = 'm'
          IUNM = 1
          CALL RDIJK( ISTART,IJK,CHDUM,UNTS,DPLLA )
          IDSPS = 1
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DPLLA(IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &      DPLLA(IROCK)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,DPLLA(IROCK),INDX)
          IF( DPLLA(IROCK).GE.SMALL ) IDSPS = 1
        ENDIF
!
!---  Transverse dispersivity  ---
!
        VARB = 'CO2 Aqueous Transverse Dispersivity: '
        IF( IJK.GT.0 ) THEN
          UNTS = 'm'
          IUNM = 1
          CALL RDIJK( ISTART,IJK,CHDUM,UNTS,DPTLA )
          IDSPS = 1
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DPTLA(IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &      DPTLA(IROCK)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,DPTLA(IROCK),INDX)
          IF( DPTLA(IROCK).GE.SMALL ) IDSPS = 1
        ENDIF
!
!---    Read next rock/soil type or scaling group  ---
!
        IF( N.LT.NROCK ) WRITE(IWR,'(/)')
        GOTO 10
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDAT_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDBC_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
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
!     11 - Hydrostatic
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
!     38 - Salt Molality
!     39 - Inflow Salt Molality
!     40 - Aqu. Conc.
!     41 - Inflow Aqu. Conc.
!     42 - Dirichlet-Inflow
!     43 - Inflow-Outflow Gas
!     44 - Geothermal Gradient
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      CHARACTER*64 ADUM,BDUM(LUK+LSOLU+2),FDUM
      CHARACTER*64 UNTS
      CHARACTER*32 CHTYP(44)
      CHARACTER*512 CHDUM
      REAL*8 VAR(LBTM,LBCV)
      REAL*8 XVX(4),YVX(4),ZVX(4)
      INTEGER ITYP(LUK+LSOLU+2)
#ifdef ecke
      CHARACTER*64 SDUM
      INTEGER IBCSPX(LSPBC+1)
#endif
!
!----------------------Data Statements---------------------------------!
!
      SAVE CHTYP
      DATA CHTYP /'Dirichlet','Neumann','Zero Flux','Saturated',
     &  'Unit Gradient','Free Gradient','Outflow',
     &  'Aqueous Concentration','Gas Concentration',
     &  'NAPL Concentration','Hydrostatic',
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
     &  'Inflow Aqu. Mass Frac.','Molality','Inflow Molality',
     &  'Aqu. Conc.','Inflow Aqu. Conc.','Dirichlet-Inflow',
     &  'Inflow-Outflow Gas','Geothermal Gradient'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDBC_CO2'
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
      VARB = 'Number of Boundary Condition Cards: '
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO 400 NB = 1, NLIN
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        IF( NB.NE.1 ) WRITE(IWR, '(/)')
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
!---    Aqueous, gas, and salt boundary condition types, 
!       isothermal option  ---
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(1))
        IF( INDEX(BDUM(1)(1:),'hydrostatic').EQ.0 ) THEN
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(2))
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 )
     &      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(3))
        ELSE
          BDUM(2) = 'hydrostatic'
          BDUM(3) = 'hydrostatic'  
        ENDIF
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
        IF( ISLC(40).EQ.1 ) THEN
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
        ENDIF
#endif
!
!---    Read aqueous boundary condition type ---
!
        IF( INDEX(BDUM(1)(1:),'dirichlet').NE.0 ) THEN
           ITYP(IEQW) = 1
        ELSEIF( INDEX(BDUM(1)(1:),'neumann').NE.0 ) THEN
           ITYP(IEQW) = 2
        ELSEIF( INDEX(BDUM(1)(1:),'zero flux').NE.0 ) THEN
           ITYP(IEQW) = 3
        ELSEIF( INDEX(BDUM(1)(1:),'saturated').NE.0 ) THEN
           ITYP(IEQW) = 4
        ELSEIF( INDEX(BDUM(1)(1:),'unit gradient').NE.0 ) THEN
           ITYP(IEQW) = 5
        ELSEIF( INDEX(BDUM(1)(1:),'hydrostatic').NE.0 ) THEN
           ITYP(IEQW) = 11
        ELSEIF( INDEX(BDUM(1)(1:),'initial cond').NE.0 ) THEN
           ITYP(IEQW) = 12
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Aqueous Boundary Condition: '//BDUM(1)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read gas boundary condition type ---
!
        IF( INDEX(BDUM(1)(1:),'hydrostatic').NE.0 ) THEN
           ITYP(IEQA) = 11
        ELSEIF( INDEX(BDUM(2)(1:),'dirichlet').NE.0 ) THEN
           ITYP(IEQA) = 1
        ELSEIF( INDEX(BDUM(2)(1:),'neumann').NE.0 ) THEN
           ITYP(IEQA) = 2
        ELSEIF( INDEX(BDUM(2)(1:),'zero flux').NE.0 ) THEN
           ITYP(IEQA) = 3
        ELSEIF( INDEX(BDUM(2)(1:),'initial cond').NE.0 ) THEN
           ITYP(IEQA) = 12
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Gas Boundary Condition: '//BDUM(2)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
!
!---      Read salt boundary condition type ---
!
          IF( INDEX(BDUM(1)(1:),'hydrostatic').NE.0 ) THEN
            ITYP(IEQS) = 11
          ELSEIF( INDEX(BDUM(3)(1:),'inflow').NE.0 .AND.
     &      INDEX(BDUM(3)(1:),'molality').NE.0 ) THEN
             ITYP(IEQS) = 39
          ELSEIF( INDEX(BDUM(3)(1:),'inflow').NE.0 .AND.
     &      INDEX(BDUM(3)(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(3)(1:),'conc').NE.0 ) THEN
             ITYP(IEQS) = 41
          ELSEIF( INDEX(BDUM(3)(1:),'inflow').NE.0 .AND.
     &      (INDEX(BDUM(3)(1:),'rel').NE.0 .OR.
     &      INDEX(BDUM(3)(1:),'sat').NE.0) ) THEN
            ITYP(IEQS) = 35
          ELSEIF( INDEX(BDUM(3)(1:),'inflow').NE.0 .AND.
     &      (INDEX(BDUM(3)(1:),'mass').NE.0 .OR.
     &      INDEX(BDUM(3)(1:),'frac').NE.0) ) THEN
            ITYP(IEQS) = 37
          ELSEIF( INDEX(BDUM(3)(1:),'molality').NE.0 ) THEN
             ITYP(IEQS) = 38
          ELSEIF( INDEX(BDUM(3)(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(3)(1:),'conc').NE.0 ) THEN
             ITYP(IEQS) = 40
          ELSEIF( INDEX(BDUM(3)(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(3)(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(3)(1:),'frac').NE.0 ) THEN
             ITYP(IEQS) = 36
          ELSEIF( INDEX(BDUM(3)(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(3)(1:),'rel').NE.0 .AND.
     &      INDEX(BDUM(3)(1:),'sat').NE.0 ) THEN
            ITYP(IEQS) = 34
          ELSEIF( INDEX(BDUM(3)(1:),'zero flux').NE.0 ) THEN
             ITYP(IEQS) = 3
          ELSEIF( INDEX(BDUM(3)(1:),'outflow').NE.0 ) THEN
             ITYP(IEQS) = 7
          ELSEIF( INDEX(BDUM(3)(1:),'initial cond').NE.0 ) THEN
             ITYP(IEQS) = 12
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Salt Boundary Condition: '//BDUM(3)
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
!---      Read gas reactive species boundary condition type(s) ---
!
          IF( INDEX(BDUM(NSOLU+LUK+2)(1:),
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
          ELSEIF( INDEX(BDUM(NSOLU+LUK+2)(1:),
     &      'volumetric conc').NE.0 ) THEN
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
        WRITE(IWR,'(2X,2A)') 'Gas: ',CHTYP(ITYP(IEQA))
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
!---      Read aqueous boundary condition variables ---
!
          IF( ITYP(IEQW).EQ.1 ) THEN
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
          ELSEIF( ITYP(IEQW).EQ.2 ) THEN
            VARB = 'Volumetric Aqueous Flux'
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
            VARB = 'Reference Aqueous CO2 Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,10))
            WRITE(IWR,'(2X,2A)') VARB(1:IVR),': ',VAR(NTM,10)
            VARB = 'Reference Dissolved-CO2 Z Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,11))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,11)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,11),INDX)
            VARB = 'Dissolved-CO2 Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,12))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,12)
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,VAR(NTM,12),INDX)
          ELSE
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          ENDIF
!
!---      Read aqueous-CO2 relative saturation ---
!
          IF( ITYP(IEQW).NE.11 ) THEN
            VARB = 'Aqueous-CO2 Relative Saturation, '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            ISTART = ISX
            ICOMMA = ICX
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
            IF( VAR(NTM,7).LT.0 ) THEN
              VARB = 'Aqueous-CO2 Mass Fraction, '
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',-VAR(NTM,7)
            ELSE
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,7)
            ENDIF
!
!---        Read gas boundary condition variables ---
!
            IF( ITYP(IEQA).EQ.1 ) THEN
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
            ELSEIF( ITYP(IEQA).EQ.2 ) THEN
              VARB = 'Volumetric Gas Flux'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,4)
              INDX = 0
              IUNM = 1
              IUNS = -1
              CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            ELSEIF( ITYP(IEQA).EQ.3 ) THEN
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
            ELSEIF( ITYP(IEQA).EQ.11 ) THEN
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
!---        Read water-vapor relative saturation ---
!
            VARB = 'Water-Vapor Relative Saturation, '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            ISTART = ISX
            ICOMMA = ICX
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,5)
!
!---        Isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Read salt boundary condition variables ---
!
              IF( ITYP(IEQS).EQ.38 .OR. ITYP(IEQS).EQ.39 ) THEN
                VARB = 'Aqueous-Salt Molality, '
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),
     &            UNTS(1:NCH),': ',VAR(NTM,6)
                INDX = 0
                IUNMOL = 1
                IUNKG = -1
                CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
              ELSEIF( ITYP(IEQS).EQ.40 .OR. ITYP(IEQS).EQ.41 ) THEN
                VARB = 'Aqueous-Salt Concentration, '
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),
     &            UNTS(1:NCH),': ',VAR(NTM,6)
                INDX = 0
                IUNKG = 1
                IUNM = -3
                CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
              ELSEIF( ITYP(IEQS).EQ.34 .OR. ITYP(IEQS).EQ.35 ) THEN
                VARB = 'Aqueous-Salt Relative Saturation, '
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,6)
              ELSEIF( ITYP(IEQS).EQ.36 .OR. ITYP(IEQS).EQ.37 ) THEN
                VARB = 'Aqueous-Salt Mass Fraction, '
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR(NTM,6)
              ELSE
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              ENDIF
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
#ifdef ecke
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
        DO 320 K = K1X,K2X
          DO 310 J = J1X,J2X
            DO 300 I = I1X,I2X
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
                GOTO 300
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
              IBCT(IEQA,NBC) = ITYP(IEQA)
              IF( ISLC(32).EQ.0 ) IBCT(IEQS,NBC) = ITYP(IEQS)
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
!---  End of RDBC_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIC_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Read input file for initial conditions information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE LEAK_WELL
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
      SUB_LOG(ISUB_LOG) = '/RDIC_CO2'
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
      IF( INDEX(ADUM(1:),'hydrostatic').NE.0 ) THEN
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
          GOTO 1001
        ELSE
          GOTO 1002
        ENDIF
        IF( NWN_LW.NE.0 ) THEN
          INDX = 4
          CHMSG = 'Non-Hydrostatic Initial Conditions with Leaky Wells'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF      
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
 1001 WRITE(IWR,'(/,A)') 'Initial Condition Variable(s) and Domain(s)'
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
        ELSEIF( INDEX( ADUM(1:),'temperature' ).NE.0 ) THEN
          VARB = 'Initial Temperature'
          IUNK = 1
        ELSEIF( INDEX( ADUM(1:),'eclipse' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'gas' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Eclipse Gas Saturation'
        ELSEIF( INDEX( ADUM(1:),'eclipse' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'pres' ).NE.0 ) THEN
          VARB = 'Initial Eclipse Pressure'
        ELSEIF( INDEX( ADUM(1:),'scaling' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'gas' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Gas Saturation Scaling'
        ELSEIF( INDEX( ADUM(1:),'aqueous sat' ).NE.0 ) THEN
          VARB = 'Initial Aqueous Saturation'
        ELSEIF( INDEX( ADUM(1:),'relative' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'trapped gas' ).NE.0 ) THEN
          VARB = 'Initial Relative Trapped Gas Saturation'
        ELSEIF( INDEX( ADUM(1:),'trapped gas' ).NE.0 ) THEN
          VARB = 'Initial Trapped Gas Saturation'
        ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-CO2 Relative Saturation'
        ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &    ( INDEX( ADUM(1:),'mass' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'frac' ).NE.0 ) ) THEN
          VARB = 'Initial Aqueous-CO2 Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'conc' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-CO2 Concentration'
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
        ELSEIF( INDEX(ADUM(1:),'aqueous press').NE.0 ) THEN
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
            CALL RDINFS( TMS,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
            DO 206 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 206
              ICBRN(N) = IVAR
  206       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( TMS,VAR(1),ADDER,IROCK,INDX )
            DO 208 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 208
              IF( IZ(N).EQ.IROCK ) ICBRN(N) = IVAR
  208       CONTINUE
          ELSE
            CALL RDINIS( TMS,VAR,ADDER,IDOM,INDX )
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
        ELSEIF( INDEX(ADUM(1:),'co2').NE.0 ) THEN
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
            CALL RDINFS( PVA,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
            DO 306 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 306
              ICAIR(N) = IVAR
  306       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PVA,VAR(1),ADDER,IROCK,INDX )
            DO 308 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 308
              IF( IZ(N).EQ.IROCK ) ICAIR(N) = IVAR
  308       CONTINUE
          ELSE
            CALL RDINIS( PVA,VAR,ADDER,IDOM,INDX )
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
!
!---    Read Eclipse gas saturation input  ---
!
        ELSEIF( INDEX(ADUM(1:),'eclipse').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'sat').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 1
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            WRITE(IWR,'(2X,3A)') 'Eclipse Gas Saturation ',
     &        'from External File: ',FDUM(1:NCHF)
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( SI,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( SI,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( SI,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            WRITE(IWR,'(2X,3A)') 'Eclipse Gas Saturation ',
     &        'from External File: ',FDUM(1:NCHF)
            CALL RDINZS( SI,VAR(1),ADDER,IROCK,INDX )
          ELSE
            WRITE(IWR,'(2X,A)') 'Eclipse Gas Saturation'
            CALL RDINIS( SI,VAR,ADDER,IDOM,INDX )
          ENDIF
!
!---    Read Eclipse pressure input  ---
!
        ELSEIF( INDEX(ADUM(1:),'eclipse').NE.0 .AND.
     &    INDEX(ADUM(1:),'pres').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            WRITE(IWR,'(2X,3A)') 'Eclipse Pressure ',
     &        'from External File: ',FDUM(1:NCHF)
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( SI,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( SI,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( SI,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            WRITE(IWR,'(2X,3A)') 'Eclipse Pressure ',
     &        'from External File: ',FDUM(1:NCHF)
            CALL RDINZS( SI,VAR(1),ADDER,IROCK,INDX )
          ELSE
            WRITE(IWR,'(2X,A)') 'Eclipse Pressure'
            CALL RDINIS( SI,VAR,ADDER,IDOM,INDX )
          ENDIF  
!                
!
!---    Read gas saturation scaling input  ---
!
        ELSEIF( INDEX(ADUM(1:),'scaling').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'sat').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            WRITE(IWR,'(2X,3A)') 'Gas Saturation Scaling',
     &        'from External File: ',FDUM(1:NCHF)
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( POSM,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( POSM,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( POSM,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            WRITE(IWR,'(2X,3A)') 'Gas Saturation Scaling ',
     &        'from External File: ',FDUM(1:NCHF)
            CALL RDINZS( POSM,VAR(1),ADDER,IROCK,INDX )
          ELSE
            WRITE(IWR,'(2X,A)') 'Gas Saturation Scaling'
            CALL RDINIS( POSM,VAR,ADDER,IDOM,INDX )
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
!              DO 512 NRC = 1,NRCK
!                NSPRX = IRC_K(1,NRC)
!                NSPPX = IRC_K(2,NRC)
!                NSPMX = IRC_K(3+NSPRX+NSPPX,NRC)
!                IF( NSPMX.EQ.NSPX ) THEN
!                  INDX = 4
!                  CHMSG = 'Solid-Species Mineral ' // 
!     &              '(see Lithology Card): ' // BDUM(1:NCHB)
!                  CALL WRMSGS( INDX )
!                ENDIF
!  512         CONTINUE
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
 1002 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIC_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINPT_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Read input file cards.
!     Direct control to card reader subroutines.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      SUB_LOG(ISUB_LOG) = '/RDINPT_CO2'
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
      CARD = 'Simulation Title Card'
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
      CARD = 'Solution Control Card'
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
      CARD = 'Grid Card'
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
      CARD = 'Rock/Soil Zonation Card'
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
        GOTO 550
      ELSE
        GOTO 500
      ENDIF
  510 CONTINUE
      INDX = 1
      CHMSG = 'Missing Inactive Nodes Card'
      CARD = 'Inactive Nodes Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for leaky well card --
!
  550 CONTINUE
  559 READ(IRD,'(A)', END= 560) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 559
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'leaky').NE.0 .AND.
     &  INDEX(CHDUM(2:),'well').NE.0 ) THEN
        CALL RDLEAK_WELL
        REWIND(IRD)
        GOTO 650
      ELSE
        GOTO 550
      ENDIF
  560 CONTINUE
      INDX = 1
      CHMSG = 'Missing Leaky Well Card'
      CARD = 'Leaky Well Card Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for mechanical properties card  ---
!
  650 CONTINUE
  659 READ(IRD,'(A)', END=660) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 659
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'mechanical').NE.0 ) THEN
        CALL RDMECH
        REWIND(IRD)
        GOTO 700
      ELSE
        GOTO 650
      ENDIF
  660 CONTINUE
      INDX = 4
      CHMSG = 'Missing Mechanical Properties Card'
      CARD = 'Mechanical Properties Card'
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
      CARD = 'Hydraulic Properties Card'
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
        CALL RDSP_CO2
        REWIND(IRD)
        GOTO 900
      ELSE
        GOTO 800
      ENDIF
  810 CONTINUE
      INDX = 4
      CHMSG = 'Missing Saturation Function Card'
      CARD = 'Saturation Function Card'
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
      CARD = 'Aqueous Relative Permeability Card'
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
        GOTO 2000
      ELSE
        GOTO 1100
      ENDIF
 1110 CONTINUE
      INDX = 4
      CHMSG = 'Missing Gas Relative Permeability Card'
      CARD = 'Gas Relative Permeability Card'
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
        CALL RDTF_CO2
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
        CARD = 'Solute/Fluid Interaction Card'
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
        CALL RDTP_CO2
        REWIND(IRD)
        GOTO 2200
      ELSE
        GOTO 2100
      ENDIF
 2110 CONTINUE
      IF( IEQC.EQ.0 .AND. ISLC(40).EQ.0 ) THEN
       REWIND(IRD)
      ELSE
        INDX = 4
        CHMSG = 'Missing Solute/Porous Media Interaction Card'
        CARD = 'Solute/Porous Media Interaction Card'
        CALL WRMSGS( INDX )
      ENDIF
 2200 CONTINUE
#ifdef ecke
!
!---  Read ECKEChem input  ---
!
      IF( ISLC(40).EQ.1 ) THEN
!
!---    Search input file for aqueous species card  ---
!
 2209   READ(IRD,'(A)', END=2210) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 2209
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'aqueous').NE.0 .AND.
     &    INDEX(CHDUM(2:),'specie').NE.0) THEN
          CALL RDAQSP
          REWIND(IRD)
          GOTO 2259
        ELSE
          GOTO 2209
        ENDIF
 2210   CONTINUE
        REWIND(IRD)
        IF( ISLC(40).GT.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Aqueous Species Card'
          CARD = 'Aqueous Species Card'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Search input file for gas species card  ---
!
 2259   READ(IRD,'(A)', END=2260) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 2259
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'gas').NE.0 .AND.
     &    INDEX(CHDUM(2:),'specie').NE.0) THEN
          CALL RDGSSP
          REWIND(IRD)
          GOTO 2300
        ELSE
          GOTO 2259
        ENDIF
 2260   CONTINUE
        REWIND(IRD)
        IF( ISLC(40).GT.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Gas Species Card'
          CARD = 'Gas Species Card'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Search input file for solid species card  ---
!
 2300   CONTINUE
 2309   READ(IRD,'(A)', END=2310) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 2309
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'solid').NE.0 .AND.
     &    INDEX(CHDUM(2:),'specie').NE.0) THEN
          CALL RDSDSP
          REWIND(IRD)
          GOTO 2400
        ELSE
          GOTO 2300
        ENDIF
 2310   CONTINUE
        REWIND(IRD)
        IF( ISLC(40).GT.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Solid Species Card'
          CARD = 'Solid Species Card'
          CALL WRMSGS( INDX )
        ENDIF
        REWIND(IRD)
 2400   CONTINUE
!
!---    Search input file for equilibrium reactions card  ---
!
 2409   READ(IRD,'(A)', END=2410) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 2409
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'equil').NE.0 .AND.
     &    INDEX(CHDUM(2:),'react').NE.0) THEN
          CALL RDEQRC
          REWIND(IRD)
          GOTO 2500
        ELSE
          GOTO 2400
        ENDIF
 2410   CONTINUE
        REWIND(IRD)
        IF( ISLC(40).GT.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Equilibrium Reactions Card'
          CARD = 'Equilibrium Reactions Card'
          CALL WRMSGS( INDX )
        ENDIF
        REWIND(IRD)
!
!---    Search input file for kinetic reactions card  ---
!
 2500   CONTINUE
 2509   READ(IRD,'(A)', END=2510) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 2509
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'kinetic').NE.0 .AND.
     &    INDEX(CHDUM(2:),'react').NE.0) THEN
          CALL RDKNRC
          REWIND(IRD)
          GOTO 2700
        ELSE
          GOTO 2500
        ENDIF
 2510   CONTINUE
        REWIND(IRD)
        IF( ISLC(40).GT.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Kinetic Reactions Card'
          CARD = 'Kinetic Reactions Card'
          CALL WRMSGS( INDX )
        ENDIF
        REWIND(IRD)
!
!---    Search input file for equilibrium equation card  ---
!
 2700   CONTINUE
 2709   READ(IRD,'(A)', END=2710) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 2709
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'equil').NE.0 .AND.
     &    INDEX(CHDUM(2:),'equat').NE.0) THEN
          CALL RDEQEQ
          REWIND(IRD)
          GOTO 2800
        ELSE
          GOTO 2700
        ENDIF
 2710   CONTINUE
        REWIND(IRD)
        IF( ISLC(40).GT.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Equilibrium Equations Card'
          CARD = 'Equilibrium Equations Card'
          CALL WRMSGS( INDX )
        ENDIF
        REWIND(IRD)
!
!---    Search input file for kinetic equations card  ---
!
 2800   CONTINUE
 2809   READ(IRD,'(A)', END=2810) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 2809
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'kinetic').NE.0 .AND.
     &    INDEX(CHDUM(2:),'equat').NE.0) THEN
          CALL RDKNEQ
          REWIND(IRD)
          GOTO 2900
        ELSE
          GOTO 2800
        ENDIF
 2810   CONTINUE
        REWIND(IRD)
        IF( ISLC(40).GT.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Kinetic Equations Card'
          CARD = 'Kinetic Equations Card'
          CALL WRMSGS( INDX )
        ENDIF
        REWIND(IRD)
!
!---    Search input file for conservation equations card  ---
!
 2900   CONTINUE
 2909   READ(IRD,'(A)', END=2910) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 2909
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'conservation').NE.0 .AND.
     &    INDEX(CHDUM(2:),'equat').NE.0) THEN
          CALL RDCNEQ
          REWIND(IRD)
          GOTO 3000
        ELSE
          GOTO 2900
        ENDIF
 2910   CONTINUE
        REWIND(IRD)
        IF( ISLC(40).GT.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Conservation Equations Card'
          CARD = 'Conservation Equations Card'
          CALL WRMSGS( INDX )
        ENDIF
        REWIND(IRD)
!
!---    Search input file for lithology card  ---
!
 3000   CONTINUE
 3009   READ(IRD,'(A)', END=3010) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 3009
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'lithol').NE.0) THEN
          CALL RDLITH
          REWIND(IRD)
          GOTO 3100
        ELSE
          GOTO 3000
        ENDIF
 3010   CONTINUE
        REWIND(IRD)
        IF( ISLC(40).GT.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Lithology Card'
          CARD = 'Lithology Card'
          CALL WRMSGS( INDX )
        ENDIF
        REWIND(IRD)
!
!---    Search input file for reactive species link card  ---
!
 3100   CONTINUE
 3109   READ(IRD,'(A)', END=3110) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 3109
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'link').NE.0 .AND.
     &    INDEX(CHDUM(2:),'specie').NE.0) THEN
          CALL RDSPLK
          REWIND(IRD)
          GOTO 3150
        ELSE
          GOTO 3100
        ENDIF
 3110   CONTINUE
        REWIND(IRD)
        IF( ISLC(40).GT.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Reactive Species Link Card'
          CARD = 'Reactive Species Link Card'
          CALL WRMSGS( INDX )
        ENDIF
        REWIND(IRD)
      ENDIF
#endif  
!
!---    Search input file for CO2 transport card --
!
 3150 CONTINUE
 3159 READ(IRD,'(A)', END=3160) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3159
      CALL LCASE( CHDUM )
      IF( CHDUM(1:15).EQ.'~co2 transport' ) THEN
        CALL RDAT_CO2
        REWIND(IRD)
        GOTO 3200
      ELSE
        GOTO 3150
      ENDIF
 3160 CONTINUE
      INDX = 1
      CHMSG = 'Missing CO2 Transport Card'
      CARD = 'CO2 Transport Card Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for salt transport card --
!
 3200 CONTINUE
 3209 READ(IRD,'(A)', END=3210) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3209
      CALL LCASE( CHDUM )
      IF( CHDUM(1:15).EQ.'~salt transport' ) THEN
        CALL RDST_CO2
        REWIND(IRD)
        GOTO 3250
      ELSE
        GOTO 3200
      ENDIF
 3210 CONTINUE
      IF( ISLC(32).EQ.0 ) THEN
        INDX = 4
        CHMSG = 'Missing Salt Transport Card'
        CARD = 'Salt Transport Card Card'
        CALL WRMSGS( INDX )
      ELSE
        INDX = 1
        CHMSG = 'Missing Salt Transport Card'
        CARD = 'Salt Transport Card Card'
        CALL WRMSGS( INDX )
        REWIND(IRD)
      ENDIF
!
!---  Search input file for coupled well card --
!
 3250 CONTINUE
 3259 READ(IRD,'(A)', END= 3260) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3259
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'coupled').NE.0 .AND.
     &  INDEX(CHDUM(2:),'well').NE.0 ) THEN
        CALL RDCOUP_WELL
        REWIND(IRD)
        GOTO 3300
      ELSE
        GOTO 3250
      ENDIF
 3260 CONTINUE
      INDX = 1
      CHMSG = 'Missing Coupled Well Card'
      CARD = 'Coupled Well Card Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for initial conditions card --
!
 3300 CONTINUE
 3309 READ(IRD,'(A)', END=3310) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3309
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'initial').NE.0 ) THEN
        CALL RDIC_CO2
        REWIND(IRD)
        GOTO 3400
      ELSE
        GOTO 3300
      ENDIF
 3310 CONTINUE
      IF( IEO.EQ.2 ) THEN
        INDX = 1
        CHMSG = 'Missing Initial Conditions Card'
        CARD = 'Initial Conditions Card'
        CALL WRMSGS( INDX )
        INDX = 2
        CALL RDRST(INDX)
        ISIC = 3
        REWIND(IRD)
      ELSE
        INDX = 4
        CHMSG = 'Missing Initial Conditions Card'
        CARD = 'Initial Conditions Card'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Search input file for boundary conditions card --
!
 3400 CONTINUE
 3409 READ(IRD,'(A)', END=3410) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3409
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'boundary').NE.0 .AND.
     &  INDEX(CHDUM(2:),'geomech').EQ.0 ) THEN
        CALL RDBC_CO2
        REWIND(IRD)
        GOTO 3420
      ELSE
        GOTO 3400
      ENDIF
 3410 CONTINUE
      INDX = 1
      CHMSG = 'Missing Boundary Conditions Card'
      CARD = 'Boundary Conditions Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
 3420 CONTINUE
!
!---  Search input file for source card --
!
 3500 CONTINUE
 3509 READ(IRD,'(A)', END= 3510) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3509
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'source').NE.0 ) THEN
        CALL RDSR_CO2
        REWIND(IRD)
        GOTO 3600
      ELSE
        GOTO 3500
      ENDIF
 3510 CONTINUE
      INDX = 1
      CHMSG = 'Missing Source Card'
      CARD = 'Source Card Card'
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
        CALL RDOU_CO2
        REWIND(IRD)
        GOTO 3700
      ELSE
        GOTO 3600
      ENDIF
 3610 CONTINUE
      INDX = 1
      CHMSG = 'Missing Output Control Card'
      CARD = 'Output Control Card'
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
        CALL RDSF_CO2
        REWIND(IRD)
        GOTO 3800
      ELSE
        GOTO 3700
      ENDIF
 3710 CONTINUE
      INDX = 1
      CHMSG = 'Missing Surface Flux Card'
      CARD = 'Surface Flux Card'
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
!---  End of RDINPT_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDOU_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Read input file for output information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
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
      USE FDVP
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
#ifdef ecke
      CHARACTER*64 SPNM
#endif
      CHARACTER*64 ADUM,UNTS,SOLNM
      CHARACTER*512 CHDUM
      CHARACTER*6 FORM,FORM1
!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM,FORM1
      DATA FORM / '(I6,$)' /
      DATA FORM1 / '(A,I1)' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDOU_CO2'
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
!
!---    Check for leaky well node  ---
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX( ADUM(1:),'leaky' ).NE.0 ) THEN
          VARB = 'Leaky Well Node Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NWN)
          NDREF(N) = -NWN
          WRITE(FORM(3:3),'(I1)') ICOUNT(NWN)
          WRITE(IWR,'(2X,A,$)') 'Reference Node No. '
          WRITE(IWR,FORM) NWN
          WRITE(IWR,'(2X,A)' ) 'Leaky Well'
          CYCLE
        ENDIF
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
          NVC = NVC - 1
          GOTO 200
  120     CONTINUE
        ENDIF
#ifdef ecke
        IF( INDEX( ADUM(1:),'species' ).NE.0 ) THEN
          IF( ISLC(40).EQ.0 ) THEN
            NVC = NVC - 1
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
          NVC = NVC - 1
          GOTO 200
  150     CONTINUE
        ENDIF
#endif
        IF( INDEX(ADUM(1:),'integr').NE.0 .AND.
     &    INDEX(ADUM(1:),'precip').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 367
        ELSEIF( INDEX(ADUM(1:),'aqueous pressure').NE.0 ) THEN
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
        ELSEIF( INDEX(ADUM(1:),'eclipse').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'sat').NE.0 ) THEN
          IREF(NV) = 83
          CHREF(83) = 'ESG'
        ELSEIF( INDEX(ADUM(1:),'eclipse').NE.0 .AND.
     &    INDEX(ADUM(1:),'pres').NE.0 ) THEN
          IREF(NV) = 84
          CHREF(84) = 'ESG'
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
        ELSEIF( INDEX(ADUM(1:),'aqueous moisture cont').NE.0 ) THEN
          IREF(NV) = 15
        ELSEIF( INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 ) THEN
          CHREF(204) = 'XMLA'
          IREF(NV) = 204
        ELSEIF( INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 ) THEN
          CHREF(205) = 'XMLS'
          IREF(NV) = 205
        ELSEIF( INDEX(ADUM(1:),'effective').NE.0 .AND.
     &    INDEX(ADUM(1:),'trap').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IREF(NV) = 19
        ELSEIF( INDEX(ADUM(1:),'diffusive porosity').NE.0 ) THEN
          IREF(NV) = 20
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 21
        ELSEIF( INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 22
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 24
        ELSEIF( INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 25
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
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 38
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'int').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 40
        ELSEIF( INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'int').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
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
        ELSEIF( INDEX(ADUM(1:),'x aqueous vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'x-aqueous vol').NE.0 ) THEN
          IREF(NV) = 51
        ELSEIF( INDEX(ADUM(1:),'y aqueous vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'y-aqueous vol').NE.0 ) THEN
          IREF(NV) = 52
        ELSEIF( INDEX(ADUM(1:),'z aqueous vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'y-aqueous vol').NE.0 ) THEN
          IREF(NV) = 53
        ELSEIF( INDEX(ADUM(1:),'x gas vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'x-gas vol').NE.0 ) THEN
          IREF(NV) = 54
        ELSEIF( INDEX(ADUM(1:),'y gas vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'y-gas vol').NE.0 ) THEN
          IREF(NV) = 55
        ELSEIF( INDEX(ADUM(1:),'z gas vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'z-gas vol').NE.0 ) THEN
          IREF(NV) = 56
        ELSEIF( INDEX(ADUM(1:),'webb').NE.0 .AND.
     &    INDEX(ADUM(1:),'matching').NE.0 .AND. 
     &    INDEX(ADUM(1:),'point').NE.0 .AND. 
     &    INDEX(ADUM(1:),'head').NE.0 ) THEN
          IREF(NV) = 63
          CHREF(63) = 'WMPH'
        ELSEIF( INDEX(ADUM(1:),'webb').NE.0 .AND.
     &    INDEX(ADUM(1:),'matching').NE.0 .AND. 
     &    INDEX(ADUM(1:),'point').NE.0 .AND. 
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IREF(NV) = 150
          CHREF(150) = 'WMPS'
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
        ELSEIF( INDEX(ADUM(1:),'h2o gas mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'water gas mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas water mole frac').NE.0 ) THEN
          IREF(NV) = 70
        ELSEIF( INDEX(ADUM(1:),'co2 gas mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas co2 mole').NE.0 ) THEN
          IREF(NV) = 71
        ELSEIF( INDEX(ADUM(1:),'h2o gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'water gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas water conc').NE.0 ) THEN
          IREF(NV) = 73
        ELSEIF( INDEX(ADUM(1:),'co2 gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas co2 conc').NE.0 ) THEN
          IREF(NV) = 74
        ELSEIF( INDEX(ADUM(1:),'h2o aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous h2o conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'water aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous water conc').NE.0 ) THEN
          IREF(NV) = 76
        ELSEIF( INDEX(ADUM(1:),'co2 aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous co2 conc').NE.0 ) THEN
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
        ELSEIF( INDEX(ADUM(1:),'node').NE.0 .AND.
     &    INDEX(ADUM(1:),'number').NE.0 ) THEN
          IREF(NV) = 100
        ELSEIF( INDEX(ADUM(1:),'osmotic pressure').NE.0 ) THEN
          IREF(NV) = 101
        ELSEIF( INDEX(ADUM(1:),'osmotic eff').NE.0 ) THEN
          IREF(NV) = 102
        ELSEIF( INDEX(ADUM(1:),'aqueous compress').NE.0 ) THEN
          IREF(NV) = 103
        ELSEIF( INDEX(ADUM(1:),'gas compress').NE.0 ) THEN
          IREF(NV) = 104
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mass').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mass').NE.0 ) THEN
          IREF(NV) = 110
        ELSEIF( INDEX(ADUM(1:),'water').NE.0 .AND.
     &    INDEX(ADUM(1:),'vapor').NE.0 .AND.
     &    INDEX(ADUM(1:),'part').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IREF(NV) = 128
        ELSEIF( INDEX(ADUM(1:),'gas-aqueous scaling').NE.0 ) THEN
          IREF(NV) = 131
        ELSEIF( INDEX(ADUM(1:),'coupled').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 .AND. 
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IREF(NV) = 138
          VARB = 'Reference Node Variable: Coupled-Well Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,IREF_CW(NV))
          IF( IREF_CW(NV).LT.1 .OR. IREF_CW(NV).GT.N_CW ) THEN
            INDX = 7
            CHMSG = 'Unrecognized Well Number: '
            IMSG = IREF_CW(NV)
            CALL WRMSGS( INDX )
          ENDIF
          CHREF(138) = 'PWCW'
          UNREF(138) = 'pa'
          WRITE(ADUM(NCH+1:NCH+11),'(A,I3)') ': well #',IREF_CW(NV)
          NCH = NCH+11
        ELSEIF( INDEX(ADUM(1:),'coupled').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 .AND. 
     &    INDEX(ADUM(1:),'mass').NE.0 .AND. 
     &    INDEX(ADUM(1:),'nodal').NE.0 .AND. 
     &    INDEX(ADUM(1:),'co2').NE.0 .AND. 
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IREF(NV) = 142
          VARB = 'Reference Node Variable: Coupled-Well Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,IREF_CW(NV))
          IF( IREF_CW(NV).LT.1 .OR. IREF_CW(NV).GT.N_CW ) THEN
            INDX = 7
            CHMSG = 'Unrecognized Well Number: '
            IMSG = IREF_CW(NV)
            CALL WRMSGS( INDX )
          ENDIF
          CHREF(142) = 'QNRA'
          UNREF(142) = 'kg/s'
          WRITE(ADUM(NCH+1:NCH+11),'(A,I3)') ': well #',IREF_CW(NV)
          NCH = NCH+11
        ELSEIF( INDEX(ADUM(1:),'coupled').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 .AND. 
     &    INDEX(ADUM(1:),'mass').NE.0 .AND. 
     &    INDEX(ADUM(1:),'nodal').NE.0 .AND. 
     &    INDEX(ADUM(1:),'water').NE.0 .AND. 
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IREF(NV) = 145
          VARB = 'Reference Node Variable: Coupled-Well Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,IREF_CW(NV))
          IF( IREF_CW(NV).LT.1 .OR. IREF_CW(NV).GT.N_CW ) THEN
            INDX = 7
            CHMSG = 'Unrecognized Well Number: '
            IMSG = IREF_CW(NV)
            CALL WRMSGS( INDX )
          ENDIF
          CHREF(145) = 'QNRW'
          UNREF(145) = 'kg/s'
          WRITE(ADUM(NCH+1:NCH+11),'(A,I3)') ': well #',IREF_CW(NV)
          NCH = NCH+11
        ELSEIF( INDEX(ADUM(1:),'coupled').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 .AND. 
     &    INDEX(ADUM(1:),'mass').NE.0 .AND. 
     &    INDEX(ADUM(1:),'co2').NE.0 .AND. 
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IREF(NV) = 349
          VARB = 'Reference Node Variable: Coupled-Well Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,IREF_CW(NV))
          IF( IREF_CW(NV).LT.1 .OR. IREF_CW(NV).GT.N_CW ) THEN
            INDX = 7
            CHMSG = 'Unrecognized Well Number: '
            IMSG = IREF_CW(NV)
            CALL WRMSGS( INDX )
          ENDIF
          CHREF(349) = 'QMRA'
          UNREF(349) = 'kg/s'
          WRITE(ADUM(NCH+1:NCH+11),'(A,I3)') ': well #',IREF_CW(NV)
          NCH = NCH+11
        ELSEIF( INDEX(ADUM(1:),'coupled').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 .AND. 
     &    INDEX(ADUM(1:),'mass').NE.0 .AND. 
     &    INDEX(ADUM(1:),'co2').NE.0 .AND. 
     &    INDEX(ADUM(1:),'integral').NE.0 ) THEN
          IREF(NV) = 350
          VARB = 'Reference Node Variable: Coupled-Well Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,IREF_CW(NV))
          IF( IREF_CW(NV).LT.1 .OR. IREF_CW(NV).GT.N_CW ) THEN
            INDX = 7
            CHMSG = 'Unrecognized Well Number: '
            IMSG = IREF_CW(NV)
            CALL WRMSGS( INDX )
          ENDIF
          CHREF(350) = 'QMIA'
          UNREF(350) = 'kg'
          WRITE(ADUM(NCH+1:NCH+11),'(A,I3)') ': well #',IREF_CW(NV)
          NCH = NCH+11
        ELSEIF( INDEX(ADUM(1:),'coupled').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 .AND. 
     &    INDEX(ADUM(1:),'mass').NE.0 .AND. 
     &    INDEX(ADUM(1:),'water').NE.0 .AND. 
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IREF(NV) = 351
          VARB = 'Reference Node Variable: Coupled-Well Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,IREF_CW(NV))
          IF( IREF_CW(NV).LT.1 .OR. IREF_CW(NV).GT.N_CW ) THEN
            INDX = 7
            CHMSG = 'Unrecognized Well Number: '
            IMSG = IREF_CW(NV)
            CALL WRMSGS( INDX )
          ENDIF
          CHREF(351) = 'QMRW'
          UNREF(351) = 'kg/s'
          WRITE(ADUM(NCH+1:NCH+11),'(A,I3)') ': well #',IREF_CW(NV)
          NCH = NCH+11
        ELSEIF( INDEX(ADUM(1:),'coupled').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 .AND. 
     &    INDEX(ADUM(1:),'mass').NE.0 .AND. 
     &    INDEX(ADUM(1:),'water').NE.0 .AND. 
     &    INDEX(ADUM(1:),'integral').NE.0 ) THEN
          IREF(NV) = 352
          VARB = 'Reference Node Variable: Coupled-Well Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,IREF_CW(NV))
          IF( IREF_CW(NV).LT.1 .OR. IREF_CW(NV).GT.N_CW ) THEN
            INDX = 7
            CHMSG = 'Unrecognized Well Number: '
            IMSG = IREF_CW(NV)
            CALL WRMSGS( INDX )
          ENDIF
          CHREF(352) = 'QMIW'
          UNREF(352) = 'kg'
          WRITE(ADUM(NCH+1:NCH+11),'(A,I3)') ': well #',IREF_CW(NV)
          NCH = NCH+11
        ELSEIF( INDEX(ADUM(1:),'h2o mass source rate').NE.0 .OR.
     &    INDEX(ADUM(1:),'water mass source rate').NE.0 ) THEN
          IREF(NV) = 140
        ELSEIF( INDEX(ADUM(1:),'co2 mass source rate').NE.0 ) THEN
          IREF(NV) = 141
        ELSEIF( INDEX(ADUM(1:),'salt mass source rate').NE.0 ) THEN
          IREF(NV) = 147
        ELSEIF( INDEX(ADUM(1:),'salt mass source int').NE.0 ) THEN
          IREF(NV) = 148
        ELSEIF( INDEX(ADUM(1:),'aqueous viscosity').NE.0 ) THEN
          IREF(NV) = 176
        ELSEIF( INDEX(ADUM(1:),'disp').NE.0 .AND.
     &    INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'vari').NE.0 ) THEN
          IREF(NV) = 132
          CHREF(132) = 'SIGXX'
        ELSEIF( INDEX(ADUM(1:),'disp').NE.0 .AND.
     &    INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'vari').NE.0 ) THEN
          IREF(NV) = 133
          CHREF(133) = 'SIGYY'
        ELSEIF( INDEX(ADUM(1:),'disp').NE.0 .AND.
     &    INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'vari').NE.0 ) THEN
          IREF(NV) = 134
          CHREF(134) = 'SIGZZ'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'trapped').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IREF(NV) = 190
          CHREF(190) = 'IMTGA'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'water').NE.0 .OR.
     &    INDEX(ADUM(1:),'h2o').NE.0 ) ) THEN
          IREF(NV) = 194
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 ) THEN
          IREF(NV) = 195
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'water').NE.0 .OR.
     &    INDEX(ADUM(1:),'h2o').NE.0 ) ) THEN
          IREF(NV) = 197
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 ) THEN
          IREF(NV) = 198
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'water').NE.0 .OR.
     &    INDEX(ADUM(1:),'h2o').NE.0 ) ) THEN
          IREF(NV) = 191
        ELSEIF( INDEX(ADUM(1:),'integrated co2 mass sour').NE.0 ) THEN
          IREF(NV) = 262
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 192
        ELSEIF( (INDEX(ADUM(1:),'x intrinsic').NE.0 .OR.
     &    INDEX(ADUM(1:),'x-intrinsic').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'perm').NE.0 ) THEN
          IREF(NV) = 247
        ELSEIF( (INDEX(ADUM(1:),'y intrinsic').NE.0 .OR.
     &    INDEX(ADUM(1:),'y-intrinsic').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'perm').NE.0 ) THEN
          IREF(NV) = 248
        ELSEIF( (INDEX(ADUM(1:),'z intrinsic').NE.0 .OR.
     &    INDEX(ADUM(1:),'z-intrinsic').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'perm').NE.0 ) THEN
          IREF(NV) = 249
        ELSEIF( INDEX(ADUM(1:),'integrated h2o mass sour').NE.0 .OR.
     &    INDEX(ADUM(1:),'integrated water mass sour').NE.0 ) THEN
          IREF(NV) = 261
        ELSEIF( INDEX(ADUM(1:),'well-node').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IREF(NV) = 275
        ELSEIF( INDEX(ADUM(1:),'source-well pres').NE.0 ) THEN
          IREF(NV) = 285
        ELSEIF( INDEX(ADUM(1:),'gas viscosity').NE.0 ) THEN
          IREF(NV) = 289
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
        ELSEIF( ( INDEX(ADUM(1:),'similarity').NE.0 .OR.
     &    INDEX(ADUM(1:),'similitude').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'variable').NE.0) THEN
          IREF(NV) = 299
        ELSEIF( INDEX(ADUM(1:),'co2 aqueous diff').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous co2 diff').NE.0 ) THEN
          IREF(NV) = 347
        ELSEIF( INDEX(ADUM(1:),'h2o gas diff').NE.0 .OR.
     &    INDEX(ADUM(1:),'water gas diff').NE.0 ) THEN
          IREF(NV) = 348
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
        ELSEIF ( INDEX(ADUM(1:),'satur').NE.0 .AND.
     &    INDEX(ADUM(1:),'func').NE.0 .AND.
     &    INDEX(ADUM(1:),'index').NE.0 ) THEN
          IREF(NV) = 394
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
        ELSEIF( INDEX(ADUM(1:),'solute aqueous mol').NE.0 )THEN
          IREF(NV) = 400 + (NSL-1)*33 + 5
        ELSEIF( INDEX(ADUM(1:),'solute gas mol').NE.0 )THEN
          IREF(NV) = 400 + (NSL-1)*33 + 6
        ELSEIF( INDEX(ADUM(1:),'xnc solute flux').NE.0 .OR.
     &    INDEX(ADUM(1:),'xnc-solute flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 8
        ELSEIF( INDEX(ADUM(1:),'ync solute flux').NE.0 .OR.
     &    INDEX(ADUM(1:),'ync-solute flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 9
        ELSEIF( INDEX(ADUM(1:),'znc solute flux').NE.0 .OR.
     &    INDEX(ADUM(1:),'znc-solute flux').NE.0 ) THEN
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
        ELSEIF( INDEX(ADUM(1:),'mineral area').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 24
          IREF(NV) = INDX
          UNREF(INDX) = 'm^2'
        ELSEIF( INDEX(ADUM(1:),'mineral rate').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 25
          IREF(NV) = INDX
          UNREF(INDX) = 'mol/s'
        ELSEIF( INDEX(ADUM(1:),'volume fraction').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 26
          IREF(NV) = INDX
          CHREF(INDX) = 'SPVF'
        ELSEIF( INDEX(ADUM(1:),'ph').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33) + 27
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
          WRITE( IWR,'(2X,3A,2X,2A,I3,A)' ) ADUM(1:NCH),', ',
     &      UNREF(IREF(NV))(1:NCU),SOLNM(1:NCS),' Solute(',NSL,')'
#ifdef ecke
        ELSEIF( INDEX( ADUM(1:),'species' ).NE.0 ) THEN
          WRITE( IWR,'(2X,3A,2X,2A,I3,A)' ) ADUM(1:NCH),', ',
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
          NVC = NVC - 1
          GOTO 400
  320     CONTINUE
        ENDIF
#ifdef ecke
        IF( INDEX( ADUM(1:),'species' ).NE.0 ) THEN
          IF( ISLC(40).EQ.0 ) THEN
            NVC = NVC - 1
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
          NVC = NVC - 1
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
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'equil').NE.0 .AND.
     &    INDEX(ADUM(1:),'inter').NE.0 .AND.
     &    INDEX(ADUM(1:),'eleva').NE.0 ) THEN
          IPLOT(NV) = 353
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'equil').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IPLOT(NV) = 354
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'equil').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IPLOT(NV) = 355
        ELSEIF( INDEX(ADUM(1:),'vertical').NE.0 .AND.
     &    INDEX(ADUM(1:),'equilibr').NE.0 .AND.
     &    INDEX(ADUM(1:),'trap').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'sat').NE.0 ) THEN
          IPLOT(NV) = 357
        ELSEIF( INDEX(ADUM(1:),'vertical').NE.0 .AND.
     &    INDEX(ADUM(1:),'equilibr').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'sat').NE.0 ) THEN
          IPLOT(NV) = 356
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'equil').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &    INDEX(ADUM(1:),'sat').NE.0 ) THEN
          IPLOT(NV) = 358
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'equil').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'rel').NE.0 .AND.
     &    INDEX(ADUM(1:),'perm').NE.0 ) THEN
          IPLOT(NV) = 359
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'equil').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &    INDEX(ADUM(1:),'rel').NE.0 .AND.
     &    INDEX(ADUM(1:),'perm').NE.0 ) THEN
          IPLOT(NV) = 360
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'int').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'area').NE.0 ) THEN
          IPLOT(NV) = 364
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'int').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IPLOT(NV) = 363
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'int').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'area').NE.0 ) THEN
          IPLOT(NV) = 366
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'int').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IPLOT(NV) = 365
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'int').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'area').NE.0 ) THEN
          IPLOT(NV) = 362
        ELSEIF( INDEX(ADUM(1:),'vert').NE.0 .AND.
     &    INDEX(ADUM(1:),'int').NE.0 .AND.
     &    INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IPLOT(NV) = 361
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
        ELSEIF( INDEX(ADUM(1:),'eclipse').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'sat').NE.0 ) THEN
          IPLOT(NV) = 83
        ELSEIF( INDEX(ADUM(1:),'eclipse').NE.0 .AND.
     &    INDEX(ADUM(1:),'pres').NE.0 ) THEN
          IPLOT(NV) = 84
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
     &    INDEX(ADUM(1:),'gas h2o mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'water gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas water mass frac').NE.0 ) THEN
          IPLOT(NV) = 21
        ELSEIF( INDEX(ADUM(1:),'co2 gas mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas co2 mass frac').NE.0 ) THEN
          IPLOT(NV) = 22
        ELSEIF( INDEX(ADUM(1:),'h2o aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous h2o mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'water aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous water mass frac').NE.0 ) THEN
          IPLOT(NV) = 24
        ELSEIF( INDEX(ADUM(1:),'co2 aqueous mass frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous co2 mass frac').NE.0 ) THEN
          IPLOT(NV) = 25
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
        ELSEIF( INDEX(ADUM(1:),'total h2o mass').NE.0 .OR.
     &    INDEX(ADUM(1:),'total water mass').NE.0 ) THEN
          IPLOT(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'total co2 mass').NE.0 ) THEN
          IPLOT(NV) = 38
        ELSEIF( INDEX(ADUM(1:),'h2o mass source int').NE.0 .OR.
     &    INDEX(ADUM(1:),'water mass source int').NE.0 ) THEN
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
        ELSEIF( INDEX(ADUM(1:),'x aqueous vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'x-aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 51
        ELSEIF( INDEX(ADUM(1:),'y aqueous vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'y-aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 52
        ELSEIF( INDEX(ADUM(1:),'z aqueous vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'z-aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 53
        ELSEIF( INDEX(ADUM(1:),'x gas vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'x-gas vol').NE.0 ) THEN
          IPLOT(NV) = 54
        ELSEIF( INDEX(ADUM(1:),'y gas vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'y-gas vol').NE.0 ) THEN
          IPLOT(NV) = 55
        ELSEIF( INDEX(ADUM(1:),'z gas vol').NE.0 .OR.
     &    INDEX(ADUM(1:),'z-gas vol').NE.0 ) THEN
          IPLOT(NV) = 56
        ELSEIF( INDEX(ADUM(1:),'webb').NE.0 .AND.
     &    INDEX(ADUM(1:),'matching').NE.0 .AND. 
     &    INDEX(ADUM(1:),'point').NE.0 .AND. 
     &    INDEX(ADUM(1:),'head').NE.0 ) THEN
          IPLOT(NV) = 63
        ELSEIF( INDEX(ADUM(1:),'webb').NE.0 .AND.
     &    INDEX(ADUM(1:),'matching').NE.0 .AND. 
     &    INDEX(ADUM(1:),'point').NE.0 .AND. 
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IPLOT(NV) = 150
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
        ELSEIF( INDEX(ADUM(1:),'h2o gas mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'water gas mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas water mole').NE.0 ) THEN
          IPLOT(NV) = 70
        ELSEIF( INDEX(ADUM(1:),'co2 gas mole').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas co2 mole').NE.0 ) THEN
          IPLOT(NV) = 71
        ELSEIF( INDEX(ADUM(1:),'h2o gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas h2o conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'water gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas water conc').NE.0 ) THEN
          IPLOT(NV) = 73
        ELSEIF( INDEX(ADUM(1:),'co2 gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'gas co2 conc').NE.0 ) THEN
          IPLOT(NV) = 74
        ELSEIF( INDEX(ADUM(1:),'h2o aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous h2o conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'water aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous water conc').NE.0 ) THEN
          IPLOT(NV) = 76
        ELSEIF( INDEX(ADUM(1:),'co2 aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous co2 conc').NE.0 ) THEN
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
        ELSEIF( INDEX(ADUM(1:),'node').NE.0 .AND.
     &    INDEX(ADUM(1:),'number').NE.0 ) THEN
          IPLOT(NV) = 100
        ELSEIF( INDEX(ADUM(1:),'osmotic press').NE.0 ) THEN
          IPLOT(NV) = 101
        ELSEIF( INDEX(ADUM(1:),'osmotic eff').NE.0 ) THEN
          IPLOT(NV) = 102
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mass').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mass').NE.0 ) THEN
          IPLOT(NV) = 110
        ELSEIF( INDEX(ADUM(1:),'water').NE.0 .AND.
     &    INDEX(ADUM(1:),'vapor').NE.0 .AND.
     &    INDEX(ADUM(1:),'part').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IPLOT(NV) = 128
        ELSEIF( INDEX(ADUM(1:),'gas-aqueous scaling').NE.0 ) THEN
          IPLOT(NV) = 131
        ELSEIF( INDEX(ADUM(1:),'h2o mass source rate').NE.0 .OR.
     &    INDEX(ADUM(1:),'water mass source rate').NE.0 ) THEN
          IPLOT(NV) = 140
        ELSEIF( INDEX(ADUM(1:),'co2 mass source rate').NE.0 ) THEN
          IPLOT(NV) = 141
        ELSEIF( INDEX(ADUM(1:),'salt mass source rate').NE.0 ) THEN
          IPLOT(NV) = 147
        ELSEIF( INDEX(ADUM(1:),'salt mass source int').NE.0 ) THEN
          IPLOT(NV) = 148
        ELSEIF( INDEX(ADUM(1:),'aqueous viscosity').NE.0 ) THEN
          IPLOT(NV) = 176
        ELSEIF( INDEX(ADUM(1:),'x intrinsic perm').NE.0  .OR.
     &    INDEX(ADUM(1:),'x-intrinsic perm').NE.0 ) THEN
          IPLOT(NV) = 247
        ELSEIF( INDEX(ADUM(1:),'y intrinsic perm').NE.0  .OR.
     &    INDEX(ADUM(1:),'y-intrinsic perm').NE.0 ) THEN
          IPLOT(NV) = 248
        ELSEIF( INDEX(ADUM(1:),'z intrinsic perm').NE.0  .OR.
     &    INDEX(ADUM(1:),'z-intrinsic perm').NE.0 ) THEN
          IPLOT(NV) = 249
        ELSEIF( INDEX(ADUM(1:),'gas viscosity').NE.0 ) THEN
          IPLOT(NV) = 289
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
        ELSEIF( ( INDEX(ADUM(1:),'similarity').NE.0 .OR.
     &    INDEX(ADUM(1:),'similitude').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'variable').NE.0) THEN
          IPLOT(NV) = 299
        ELSEIF( INDEX(ADUM(1:),'co2 aqueous diff').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous co2 diff').NE.0 ) THEN
          IPLOT(NV) = 347
        ELSEIF( INDEX(ADUM(1:),'h2o gas diff').NE.0 .OR.
     &    INDEX(ADUM(1:),'water gas diff').NE.0 ) THEN
          IPLOT(NV) = 348
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
        ELSEIF( INDEX(ADUM(1:),'x-dir').NE.0 .AND.
     &    INDEX(ADUM(1:),'index').NE.0 ) THEN
          IPLOT(NV) = 388
        ELSEIF( INDEX(ADUM(1:),'y-dir').NE.0 .AND.
     &    INDEX(ADUM(1:),'index').NE.0 ) THEN
          IPLOT(NV) = 389
        ELSEIF( INDEX(ADUM(1:),'z-dir').NE.0 .AND.
     &    INDEX(ADUM(1:),'index').NE.0 ) THEN
          IPLOT(NV) = 390
        ELSEIF( INDEX(ADUM(1:),'xnc').NE.0 .AND.
     &    INDEX(ADUM(1:),'surface').NE.0 .AND.
     &    INDEX(ADUM(1:),'area').NE.0 ) THEN
          IPLOT(NV) = 391
        ELSEIF( INDEX(ADUM(1:),'ync').NE.0 .AND.
     &    INDEX(ADUM(1:),'surface').NE.0 .AND.
     &    INDEX(ADUM(1:),'area').NE.0 ) THEN
          IPLOT(NV) = 392
        ELSEIF( INDEX(ADUM(1:),'znc').NE.0 .AND.
     &    INDEX(ADUM(1:),'surface').NE.0 .AND.
     &    INDEX(ADUM(1:),'area').NE.0 ) THEN
          IPLOT(NV) = 393
        ELSEIF ( INDEX(ADUM(1:),'satur').NE.0 .AND.
     &    INDEX(ADUM(1:),'func').NE.0 .AND.
     &    INDEX(ADUM(1:),'index').NE.0 ) THEN
          IPLOT(NV) = 394
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
        ELSEIF( INDEX(ADUM(1:),'solute volumetric conc').NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 1
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
        ELSEIF( INDEX(ADUM(1:),'solute aqueous conc').NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 2
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
     &    INDEX( SPNM(1:),'total_' ).NE.0 .AND.
     &    INDEX( SPNM(1:),'solid' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 28
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 28
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
        ELSEIF( ( INDEX(ADUM(1:),'solute gas conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'species gas conc').NE.0 ) .AND.
     &    INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 3
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 3
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF
#else
        ELSEIF( INDEX(ADUM(1:),'solute gas conc').NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 3
#endif
        ELSEIF( INDEX(ADUM(1:),'solute aqueous mol').NE.0 )THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 5
        ELSEIF( INDEX(ADUM(1:),'solute gas mol').NE.0 )THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 6
        ELSEIF( INDEX(ADUM(1:),'x solute flux').NE.0  .OR.
     &    INDEX(ADUM(1:),'x-solute flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 8
        ELSEIF( INDEX(ADUM(1:),'y solute flux').NE.0  .OR.
     &    INDEX(ADUM(1:),'y-solute flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 9
        ELSEIF( INDEX(ADUM(1:),'z solute flux').NE.0  .OR.
     &    INDEX(ADUM(1:),'z-solute flux').NE.0 ) THEN
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
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33) + 27
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
!---  End of RDOU_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSF_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Read input file surface flux information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
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
      CHARACTER*64 ADUM,BDUM,FDUM,GDUM
      CHARACTER*512 CHDUM,CHDUMX
      LOGICAL FLG_CHK
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSF_CO2'
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
          IF( NSFGP.GT.LSF ) THEN
            INDX = 4
            CHMSG = 'Number of Surface Flux Files > LSF'
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
        WRITE(IWR,'(A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'aqueous').NE.0) THEN
          IF( INDEX(ADUM(1:),'co2').NE.0) THEN
            ISFT(NS) = 29
            WRITE(IWR,'(A)') 'Aqueous-Phase CO2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'volum').NE.0) THEN
            ISFT(NS) = 2
            WRITE(IWR,'(A)') 'Aqueous-Phase Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ELSE
            ISFT(NS) = 5
            WRITE(IWR,'(A)') 'Aqueous-Phase Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'co2').NE.0) THEN
            ISFT(NS) = 28
            WRITE(IWR,'(A)') 'Gas-Phase CO2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'volum').NE.0) THEN
            ISFT(NS) = 3
            WRITE(IWR,'(A)') 'Gas-Phase Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ELSE
            ISFT(NS) = 6
            WRITE(IWR,'(A)') 'Gas-Phase Mass Flux Surface'
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
            CHMSG = 'Unrecognized Solute Name: '//BDUM
            CALL WRMSGS( INDX )
   12     CONTINUE
#ifdef ecke
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
        ELSEIF( INDEX(ADUM(1:),'file').NE.0 ) THEN
          CALL RDCHR(ISTART,ICOMMA,NCHF,CHDUM,FDUM)
          NCHF = INDEX(FDUM,'  ')-1
          INQUIRE( FILE=FDUM(1:NCHF), FORM=GDUM, EXIST=FLG_CHK )
          IF( .NOT.FLG_CHK ) THEN
            INDX = 4
            CHMSG = 'Surface-Flux-Domain File: '
     &        // FDUM(1:NCHF)
            CALL WRMSGS( INDX )
          ELSEIF( GDUM.EQ.'UNFORMATTED' ) THEN
            INDX = 4
            CHMSG = 'Unformatted Surface-Flux-Domain File: '
     &        // FDUM(1:NCHF)
            CALL WRMSGS( INDX )
          ENDIF
          OPEN(UNIT=26,FILE=FDUM(1:NCHF),STATUS='OLD',FORM='FORMATTED')
          WRITE(IWR,'(/,2A)') 'Surface-Flux-Domain File: ',
     &      FDUM(1:NCHF)
          ISFD(NS) = 4
          ISFC(1,NS) = 1
          ISFC(2,NS) = 1
          ISFC(3,NS) = 1
          ISFC(4,NS) = 1
          ISFC(5,NS) = 1
          ISFC(6,NS) = 1
          NC = 0
   30     CONTINUE
          READ(26,*,END=40) IX,JX,KX,ISFDX
          NC = NC + 1
          IF( NC.GT.LSFDOM ) THEN
            INDX = 5
            CHMSG = 'Number of Surface-Flux-Domain Surfaces ' //
     &        '> Parameter LSFDOM'
            CALL WRMSGS( INDX )
          ENDIF
          ISFDOM(1,NC,NS) = IX
          ISFDOM(2,NC,NS) = JX
          ISFDOM(3,NC,NS) = KX
          ISFDOM(4,NC,NS) = ISFDX
          GOTO 30
   40     CONTINUE
          NSFDOM(NS) = NC
          CLOSE(26)
        ENDIF
!
!---    Check surface flux domain  ---
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
   50     CONTINUE  

        ELSE
!
!---      Read and check surface flux domain, checking for leaky
!         well node first  ---
!
          ISTX = ISTART
          ICMX = ICOMMA
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
          IF( INDEX( ADUM(1:),'leaky' ).NE.0 ) THEN
            IF( ABS(ISFD(NS)).EQ.2 ) THEN
              INDX = 4
              CHMSG = 'Unrecognized Leaky Well Surface Direction'
              CALL WRMSGS( INDX )
            ENDIF
            VARB = 'Starting Leaky Well Node Number'
            CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(1,NS))
            VARB = 'Ending Leaky Well Node Number'
            CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(2,NS))
            IF( ISFC(1,NS).LT.1 .OR. ISFC(1,NS).GT.LWN_LW .OR.
     &        ISFC(2,NS).LT.1 .OR. ISFC(2,NS).GT.LWN_LW .OR.
     &        ISFC(1,NS).GT.ISFC(2,NS) ) THEN
              INDX = 4
              CHMSG = 'Illegal Surface Flux Domain: Leaky Well'
              CALL WRMSGS( INDX )
            ENDIF
            ISFC(1,NS) = MAX( 1,ISFC(1,NS) )
            ISFC(1,NS) = MIN( LWN_LW,ISFC(1,NS),ISFC(2,NS) )
            ISFC(2,NS) = MAX( 1,ISFC(1,NS),ISFC(2,NS) )
            ISFC(2,NS) = MIN( LWN_LW,ISFC(2,NS) )
            WRITE(IWR,'(/,A)') 'Leaky Well Node Range'
            WRITE(IWR, '(2X,2(A,I6))') 'Leak Well Node = ',ISFC(1,NS),
     &        ' to ',ISFC(2,NS)
            ISFC(1,NS) = -ISFC(1,NS)
            ISFC(2,NS) = -ISFC(2,NS)
            ISFC(3,NS) = 1
            ISFC(4,NS) = 1
            ISFC(5,NS) = 1
            ISFC(6,NS) = 1
          ELSE
            ISTART = ISTX
            ICOMMA = ICMX
            VARB = 'Surface Flux Domain: '
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
            WRITE(IWR,'(/,A)') VARB(1:IVR)
            WRITE(IWR, '(2X,2(A,I6))') 'I = ',ISFC(1,NS),' to ',
     &        ISFC(2,NS)
            WRITE(IWR, '(2X,2(A,I6))') 'J = ',ISFC(3,NS),' to ',
     &        ISFC(4,NS)
            WRITE(IWR, '(2X,2(A,I6))') 'K = ',ISFC(5,NS),' to ',
     &        ISFC(6,NS)
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSF_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSP_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Read input file for rock/soil saturation function information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      SUB_LOG(ISUB_LOG) = '/RDSP_CO2'
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
        ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &    INDEX(ADUM(1:),'corey').NE.0 ) THEN
          ISCHRX = 2
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
        IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
          VARB = 'van Genuchten (Reference Interfacial Tension)'
          TX = 20.D+0
          XLSX = 0.D+0
          IF( IJK.GT.0 ) THEN
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            INDX = 16
            LNDX = LSCHR
            UNTS = 'n/m'
            IUNKG = 1
            IUNS = -2
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL SFT_L( TX,XLSX,SCHRX(16) )
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(16))
            UNTS = 'n/m'
            IDFLT = 1
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',SCHRX(16)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHRX(16),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(16),', N/m)'
          ENDIF
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
            DO 242 N = 1,NFLD
              SCHR(12,N) = HDOD
  242       CONTINUE
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
        IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
          VARB = 'Brooks and Corey (Reference Interfacial Tension)'
          TX = 20.D+0
          XLSX = 0.D+0
          IF( IJK.GT.0 ) THEN
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            INDX = 16
            LNDX = LSCHR
            UNTS = 'n/m'
            IUNKG = 1
            IUNS = -2
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL SFT_L( TX,XLSX,SCHRX(16) )
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(16))
            UNTS = 'n/m'
            IDFLT = 1
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',SCHRX(16)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHRX(16),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(16),', N/m)'
          ENDIF
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
            DO 246 N = 1,NFLD
              SCHR(12,N) = HDOD
  246       CONTINUE
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
        IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
          VARB = 'Dual Porosity van Genuchten ' //
     &      '(Reference Interfacial Tension)'
          TX = 20.D+0
          XLSX = 0.D+0
          IF( IJK.GT.0 ) THEN
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            INDX = 16
            LNDX = LSCHR
            UNTS = 'n/m'
            IUNKG = 1
            IUNS = -2
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL SFT_L( TX,XLSX,SCHRX(16) )
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(16))
            UNTS = 'n/m'
            IDFLT = 1
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',SCHRX(16)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHRX(16),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(16),', N/m)'
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
        IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
          VARB = 'Dual Porosity Brooks and Corey ' //
     &      '(Reference Interfacial Tension)'
          TX = 20.D+0
          XLSX = 0.D+0
          IF( IJK.GT.0 ) THEN
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            INDX = 16
            LNDX = LSCHR
            UNTS = 'n/m'
            IUNKG = 1
            IUNS = -2
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL SFT_L( TX,XLSX,SCHRX(16) )
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(16))
            UNTS = 'n/m'
            IDFLT = 1
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',SCHRX(16)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHRX(16),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(16),', N/m)'
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
        IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
          VARB = 'Haverkamp (Reference Interfacial Tension)'
          TX = 20.D+0
          XLSX = 0.D+0
          IF( IJK.GT.0 ) THEN
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            INDX = 16
            LNDX = LSCHR
            UNTS = 'n/m'
            IUNKG = 1
            IUNS = -2
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL SFT_L( TX,XLSX,SCHRX(16) )
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(16))
            UNTS = 'n/m'
            IDFLT = 1
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',SCHRX(16)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHRX(16),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(16),', N/m)'
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
        IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
          VARB = 'Russo (Reference Interfacial Tension)'
          TX = 20.D+0
          XLSX = 0.D+0
          IF( IJK.GT.0 ) THEN
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            INDX = 16
            LNDX = LSCHR
            UNTS = 'n/m'
            IUNKG = 1
            IUNS = -2
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL SFT_L( TX,XLSX,SCHRX(16) )
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(16))
            UNTS = 'n/m'
            IDFLT = 1
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',SCHRX(16)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHRX(16),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(16),', N/m)'
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
!---      Reference interfacial tension option  ---
!
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
            VARB = 'Tabular (Reference Interfacial Tension)'
            TX = 20.D+0
            XLSX = 0.D+0
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            INDX = 16
            LNDX = LSCHR
            UNTS = 'n/m'
            IUNKG = 1
            IUNS = -2
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ENDIF
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
!!
!!---        Build cubic splines  ---
!!
!            IF( ISCHRX.EQ.11 ) THEN
!              CALL SPLINY( NTBLX,NTBL )
!              CALL SPLINX( NTBLX,NTBL )
!            ENDIF
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
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
            VARB = 'Tabular (Reference Interfacial Tension)'
            TX = 20.D+0
            XLSX = 0.D+0
            CALL SFT_L( TX,XLSX,SCHRX(16) )
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(16))
            UNTS = 'n/m'
            IDFLT = 1
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',SCHRX(16)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHRX(16),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(16),', N/m)'
          ENDIF
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
!          IF( ISCHRX.EQ.11 ) THEN
!            CALL SPLINY( ISLTBL(1,IROCK),ISLTBL(2,IROCK) )
!            CALL SPLINX( ISLTBL(1,IROCK),ISLTBL(2,IROCK) )
!          ENDIF
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
!---      Reference interfacial tension option  ---
!
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
            VARB = 'Tabular (Reference Interfacial Tension)'
            TX = 20.D+0
            XLSX = 0.D+0
            CALL SFT_L( TX,XLSX,SCHR(16,IROCK) )
            INDX = 16
            LNDX = LSCHR
            UNTS = 'n/m'
            IUNKG = 1
            IUNS = -2
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ENDIF
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
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
            VARB = 'Tabular (Reference Interfacial Tension)'
            TX = 20.D+0
            XLSX = 0.D+0
            CALL SFT_L( TX,XLSX,SCHRX(16) )
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(16))
            UNTS = 'n/m'
            IDFLT = 1
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',SCHRX(16)
            INDX = 0
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SCHRX(16),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(16),', N/m)'
          ENDIF
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
!!
!!---        Build cubic splines  ---
!!
!          IF( ISCHRX.EQ.13 ) THEN
!            CALL SPLINY( ISLTBL(1,IROCK),ISLTBL(2,IROCK) )
!            CALL SPLINX( ISLTBL(1,IROCK),ISLTBL(2,IROCK) )
!          ENDIF
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
!!
!!---        Build cubic splines  ---
!!
!          IF( ISCHRX.EQ.13 ) THEN
!            CALL SPLINY( ISLTBL(3,IROCK),ISLTBL(4,IROCK) )
!            CALL SPLINX( ISLTBL(3,IROCK),ISLTBL(4,IROCK) )
!          ENDIF
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSP_CO2 group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSR_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Read input file for source information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      SUB_LOG(ISUB_LOG) = '/RDSR_CO2'
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
        IF( INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &    INDEX(ADUM(1:),'vol').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Volumetric Source'
          ISRTX = 3
          VARB = 'Aqueous-Salt Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(1:),'conc').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-Salt Concentration'
            ISRTX = ISRTX + 100
          ELSEIF( INDEX(BDUM(1:),'rel').NE.0 .AND.
     &      INDEX(BDUM(1:),'sat').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-Salt Relative Saturation'
            ISRTX = ISRTX + 200
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-Salt Mass Fraction'
            ISRTX = ISRTX + 300
          ELSEIF( INDEX(BDUM(1:),'molality').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-Salt Molality'
            ISRTX = ISRTX + 400
          ELSE
            WRITE(IWR,'(2X,A)') 'No Aqueous Salt'
          ENDIF
          VARB = 'Aqueous-CO2 Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHC,CHDUM,CDUM)
          IF( INDEX(CDUM(1:),'aqu').NE.0 .AND.
     &      INDEX(CDUM(1:),'conc').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-CO2 Concentration'
            ISRTX = ISRTX + 1000
          ELSEIF( INDEX(CDUM(1:),'mass').NE.0 .AND.
     &      INDEX(CDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-CO2 Mass Fraction'
            ISRTX = ISRTX + 3000
          ELSEIF( INDEX(CDUM(1:),'rel').NE.0 .AND.
     &      INDEX(CDUM(1:),'sat').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-CO2 Rel. Saturation'
            ISRTX = ISRTX + 2000
          ELSE
            WRITE(IWR,'(2X,A)') 'No Aqueous CO2'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'vol').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Gas Volumetric Source'
          VARB = 'Water-Vapor Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
            ISRTX = 4
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
            ISRTX = 5
          ELSE
            ISRTX = 6
            WRITE(IWR,'(2X,A)') 'No Water-Vapor'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Mass Source'
          ISRTX = 7
          VARB = 'Aqueous-Salt Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'aqu').NE.0 .AND.
     &      INDEX(BDUM(1:),'conc').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-Salt Concentration'
            ISRTX = ISRTX + 100
          ELSEIF( INDEX(BDUM(1:),'rel').NE.0 .AND.
     &      INDEX(BDUM(1:),'sat').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-Salt Rel. Saturation'
            ISRTX = ISRTX + 200
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-Salt Mass Fraction'
            ISRTX = ISRTX + 300
          ELSEIF( INDEX(BDUM(1:),'molality').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-Salt Molality'
            ISRTX = ISRTX + 400
          ELSE
            WRITE(IWR,'(2X,A)') 'No Aqueous Salt'
          ENDIF
          VARB = 'Aqueous-CO2 Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHC,CHDUM,CDUM)
          IF( INDEX(CDUM(1:),'aqu').NE.0 .AND.
     &      INDEX(CDUM(1:),'conc').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-CO2 Concentration'
            ISRTX = ISRTX + 1000
          ELSEIF( INDEX(CDUM(1:),'mass').NE.0 .AND.
     &      INDEX(CDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-CO2 Mass Fraction'
            ISRTX = ISRTX + 3000
          ELSEIF( INDEX(CDUM(1:),'rel').NE.0 .AND.
     &      INDEX(CDUM(1:),'sat').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous-CO2 Rel. Saturation'
            ISRTX = ISRTX + 2000
          ELSE
            WRITE(IWR,'(2X,A)') 'No Aqueous CO2'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Gas Mass Source'
          VARB = 'Water-Vapor Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
            ISRTX = 8
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
            ISRTX = 9
          ELSE
            ISRTX = 10
            WRITE(IWR,'(2X,A)') 'No Water-Vapor'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'salt').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'density').NE.0 ) THEN
            ISRTX = 11
            WRITE(IWR,'(2X,A)') 'Salt Density Source'
          ELSE
            ISRTX = 12
            WRITE(IWR,'(2X,A)') 'Salt Source'
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
        DO K = K1X,K2X
        DO J = J1X,J2X
        DO I = I1X,I2X
          IF( IXP(ND(I,J,K)).EQ.0 ) THEN
            INDX = 7
            IMSG = ND(I,J,K)
            CHMSG = 'Source Applied to an Inactive Node: '
            CALL WRMSGS( INDX )
          ENDIF
        ENDDO
        ENDDO
        ENDDO
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
!---      Aqueous Volumetric Source  ---
!
          IF( MOD(ISRTX,100).EQ.3 ) THEN
            VARB = 'Source Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Source Aqueous Volumetric Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            IF( MOD(ISRTX,1000)/100.EQ.1 ) THEN
              VARB = 'Source Aqueous-Salt Concentration'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,8)
              INDX = 0
              IUNM = -3
              IUNKG = 1
              CALL RDUNIT(UNTS,VAR(NTM,8),INDX)
              WRITE(IWR,'(2X,A,1PE11.4,A)') ' (',VAR(NTM,8),', kg/m^3)'
            ELSEIF( MOD(ISRTX,1000)/100.EQ.2 ) THEN
              VARB = 'Source Aqueous-Salt Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
              WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,8)
            ELSEIF( MOD(ISRTX,1000)/100.EQ.3 ) THEN
              VARB = 'Source Aqueous-Salt Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
              WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,8)
            ELSEIF( MOD(ISRTX,1000)/100.EQ.4 ) THEN
              VARB = 'Source Aqueous-Salt Molality'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
              WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,8)
            ENDIF
            IF( (ISRTX/1000).EQ.1 ) THEN
              VARB = 'Source Aqueous-CO2 Concentration'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,5)
              INDX = 0
              IUNM = -3
              IUNKG = 1
              CALL RDUNIT(UNTS,VAR(NTM,5),INDX)
              WRITE(IWR,'(2X,A,1PE11.4,A)') ' (',VAR(NTM,5),', kg/m^3)'
            ELSEIF( (ISRTX/1000).EQ.2 ) THEN
              VARB = 'Source Aqueous-CO2 Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( (ISRTX/1000).EQ.3 ) THEN
              VARB = 'Source Aqueous-CO2 Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,5)
            ENDIF
!
!---      Gas Volumetric Source  ---
!
          ELSEIF( ISRTX.EQ.4 .OR. ISRTX.EQ.5 .OR. ISRTX.EQ.6 ) THEN
            VARB = 'Source Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Source Gas Volumetric Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            IF( ISRTX.EQ.4 ) THEN
              VARB = 'Source Water Vapor Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( ISRTX.EQ.5 ) THEN
              VARB = 'Source Water Vapor Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( ISRTX.EQ.6 ) THEN
              ISRTX = 5
              VAR(NTM,5) = 0.D+0
            ENDIF
!
!---      Aqueous Mass Source  ---
!
          ELSEIF( MOD(ISRTX,10).EQ.7 ) THEN
            VARB = 'Source Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Source Aqueous Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            IF( MOD(ISRTX,100)/10.EQ.1 ) THEN
              VARB = 'Source Aqueous-Salt Concentration'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,8)
              INDX = 0
              IUNM = -3
              IUNKG = 1
              CALL RDUNIT(UNTS,VAR(NTM,8),INDX)
              WRITE(IWR,'(2X,A,1PE11.4,A)') ' (',VAR(NTM,8),', kg/m^3)'
            ELSEIF( MOD(ISRTX,1000)/100.EQ.2 ) THEN
              VARB = 'Source Aqueous-Salt Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
              WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,8)
            ELSEIF( MOD(ISRTX,1000)/100.EQ.3 ) THEN
              VARB = 'Source Aqueous-Salt Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
              WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,8)
            ELSEIF( MOD(ISRTX,1000)/100.EQ.4 ) THEN
              VARB = 'Source Aqueous-Salt Molality'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
              WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,8)
            ENDIF
            IF( (ISRTX/1000).EQ.1 ) THEN
              VARB = 'Source Aqueous-CO2 Concentration'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,5)
              INDX = 0
              IUNM = -3
              IUNKG = 1
              CALL RDUNIT(UNTS,VAR(NTM,5),INDX)
              WRITE(IWR,'(2X,A,1PE11.4,A)') ' (',VAR(NTM,5),', kg/m^3)'
            ELSEIF( (ISRTX/1000).EQ.2 ) THEN
              VARB = 'Source Aqueous-CO2 Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( (ISRTX/1000).EQ.3 ) THEN
              VARB = 'Source Aqueous-CO2 Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,5)
            ENDIF
!
!---      Gas Mass Source  ---
!
          ELSEIF( ISRTX.EQ.8 .OR. ISRTX.EQ.9 .OR. ISRTX.EQ.10 ) THEN
            VARB = 'Source Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Source Gas Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            IF( ISRTX.EQ.8 ) THEN
              VARB = 'Source Water Vapor Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( ISRTX.EQ.9 ) THEN
              VARB = 'Source Water Vapor Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( ISRTX.EQ.10 ) THEN
              ISRTX = 9
              VAR(NTM,5) = 0.D+0
            ENDIF
!
!---      Salt Density Source  ---
!
          ELSEIF( ISRTX.EQ.11 ) THEN
            VARB = 'Source Salt Density Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
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
            WRITE(IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
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
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
          ELSEIF( ISRTX.LT.-NSOLU .AND. ISRTX.GE.-2*NSOLU ) THEN
            VARB = 'Source Solute Density Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
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
!---  End of RDSR_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDST_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Reads the salt transport card dispersivities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      SUB_LOG(ISUB_LOG) = '/RDST_CO2'
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDST_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTF_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Reads solute/fluid interaction card for diffusion and partition
!     coefficients, and internodal diffusion term averaging scheme for
!     single phase (aqueous) solute transport equation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      SUB_LOG(ISUB_LOG) = '/RDTF_CO2'
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
!---    Gas-phase molecular diffusion coefficient  ---
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
!---    Gas-aqueous partition coefficient option  ---
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
!---  End of RDTF_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTP_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Reads the solute/porous media interaction card for the
!     dispersivities, half-lives, and partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      SUB_LOG(ISUB_LOG) = '/RDTP_CO2'
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
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Solute Name'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---    Search known solutes for matching name  ---
!
        DO 300 NSL = 1,NSOLU
          IF( ADUM.EQ.SOLUT(NSL)) GOTO 400
  300   CONTINUE
        INDX = 4
        CHMSG = 'Unrecognized Solute Name: '//ADUM
        CALL WRMSGS( INDX )
  400   CONTINUE
        WRITE(IWR,'(/,2A)') 'Solute Name:',SOLUT(NSL)
!
!---    Solid-aqueous partition coefficient  ---
!
        VARB = 'Solid-aqueous Adsorption Function'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'linear kd').NE.0 ) THEN
          VARB = 'Kd parameter'
          WRITE(IWR,'(A)') 'Linear Isotherm w/ Kd'
          UNTS = 'm^3/kg'
          IUNKG = -1
          IUNM = 3
          IF( IJK.GT.0 ) THEN
            INDX = 1
            LNDX = 5
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,PCSL(1,1,NSL),INDX,LNDX )
            DO 402 IROCK = 1,NFLD
              PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
              IPCSL(IROCK,NSL) = 1
  402       CONTINUE
          ELSE
            WRITE(IWR,'(A)') 'Linear Isotherm w/ Kd'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PCSL(1,IROCK,NSL)
            INDX = 0
            CALL RDUNIT(UNTS,PCSL(1,IROCK,NSL),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PCSL(1,IROCK,NSL),
     &        ', m^3/kg)'
            PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
            IPCSL(IROCK,NSL) = 1
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'linear').NE.0 ) THEN
          WRITE(IWR,'(A)') 'Linear Isotherm'
          VARB = 'k parameter'
          UNTS = 'null'
          IF( IJK.GT.0 ) THEN
            INDX = 1
            LNDX = 5
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,PCSL(1,1,NSL),INDX,LNDX )
            DO 404 IROCK = 1,NFLD
              PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
              IPCSL(IROCK,NSL) = 0
  404       CONTINUE
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
            WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCSL(1,IROCK,NSL)
            PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
            IPCSL(IROCK,NSL) = 0
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'freundlich').NE.0 ) THEN
          WRITE(IWR,'(A)') 'Freundlich Isotherm'
          VARB = 'k parameter'
          UNTS = 'null'
          IF( IJK.GT.0 ) THEN
            INDX = 1
            LNDX = 5
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,PCSL(1,1,NSL),INDX,LNDX )
            DO 406 IROCK = 1,NFLD
              PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
              IPCSL(IROCK,NSL) = 2
  406       CONTINUE
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
            WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCSL(1,IROCK,NSL)
            PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
            IPCSL(IROCK,NSL) = 2
          ENDIF
          VARB = 'n Parameter'
          UNTS = 'null'
          IF( IJK.GT.0 ) THEN
            INDX = 2
            LNDX = 5
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,PCSL(1,1,NSL),INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(2,IROCK,NSL))
            WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCSL(2,IROCK,NSL)
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'langmuir').NE.0 ) THEN
          WRITE(IWR,'(A)') 'Langmuir Isotherm'
          VARB = 'a parameter'
          UNTS = 'null'
          IF( IJK.GT.0 ) THEN
            INDX = 1
            LNDX = 5
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,PCSL(1,1,NSL),INDX,LNDX )
            DO 408 IROCK = 1,NFLD
              PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
              IPCSL(IROCK,NSL) = 3
  408       CONTINUE
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
            WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCSL(1,IROCK,NSL)
            PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
            IPCSL(IROCK,NSL) = 3
          ENDIF
          VARB = 'b Parameter'
          UNTS = 'm^3'
          IUNM = 3
          IF( IJK.GT.0 ) THEN
            INDX = 2
            LNDX = 5
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,PCSL(1,1,NSL),INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(2,IROCK,NSL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PCSL(2,IROCK,NSL)
            INDX = 0
            CALL RDUNIT(UNTS,PCSL(2,IROCK,NSL),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PCSL(2,IROCK,NSL),', m^3)'
          ENDIF
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
!---  End of RDTP_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKGS_CO2( ASL_F,ASL_M,ASLX,ASLMINX,BTGLX,ESGTX,PGX,PLX,
     &  RKGX,SGX,SGTX,SLX,N )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Gas relative permeability.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/RKGS_CO2'
      IZN = IZ(N)
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
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 ) THEN
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
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 ) THEN
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
        SGZ = (1.D+0-ASLZ)*(1.D+0-SLRX)
        SGPX = MIN( MAX( (SGZ-SGRX)/(1.D+0-SGRX),0.D+0 ),1.D+0 )
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
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKGX = FNTBLY( HDGLX,IRGTBL(1,IZN),IRGTBL(2,IZN),ITBX )
!
!---  Tabular gas relative permeability versus ln(capillary head)
!     w/ cubic spline interpolation  ---
!
      ELSEIF( IRPG(IZN).EQ.13 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKGX = FSPLNY( HDGLX,IRGTBL(1,IZN),IRGTBL(2,IZN) )
!
!---  Tabular gas relative permeability versus capillary head
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( IRPG(IZN).EQ.14 ) THEN
        ITBX = 0
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKGX = FNTBLY( HDGL,IRGTBL(1,IZN),IRGTBL(2,IZN),ITBX )
!
!---  Tabular gas relative permeability versus capillary head
!     w/ cubic spline interpolation  ---
!
      ELSEIF( IRPG(IZN).EQ.15 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
!
!---  Mualem porosity distribution function w/ van Genuchten 
!     saturation function  ---
!
      ELSEIF( IRPG(IZN).EQ.41 ) THEN
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        SGP = 1.D+0 - ASLZ
        RKGX = SQRT(SGP)*((1.D+0-ASLZ**(1.D+0/RPGC(3,IZN)))
     &    **RPGC(3,IZN))**2
!
!---  Mualem porosity distribution function w/ Brooks and Corey 
!     saturation function  ---
!
      ELSEIF( IRPG(IZN).EQ.42 ) THEN
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        SGP = 1.D+0 - ASLZ
        RKGX = SQRT(SGP)*(1.D+0-ASLZ**(1.D+0+1.D+0/RPGC(3,IZN)))**2
!
!---  Burdine porosity distribution function w/ van Genuchten 
!     saturation function  ---
!
      ELSEIF( IRPG(IZN).EQ.43 ) THEN
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        SGP = 1.D+0 - ASLZ
        RKGX = (SGP**2)*((1.D+0-ASLZ**(1.D+0/RPGC(3,IZN)))
     &    **RPGC(3,IZN))
!
!---  Burdine porosity distribution function w/ Brooks and Corey 
!     saturation function  ---
!
      ELSEIF( IRPG(IZN).EQ.44 ) THEN
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        SGP = 1.D+0 - ASLZ
        RKGX = (SGP**2)*(1.D+0-ASLZ**(1.D+0 + 2.D+0/RPGC(3,IZN)))
      ENDIF
      RKGX = MIN( MAX( RKGX,0.D+0 ),1.D+0 )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RKGS_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKLS_CO2( ASL_F,ASL_M,BTGLX,ESLX,PGX,PLX,
     &  RKLX,SLX,N )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Aqueous relative permeability.
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
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RKLS_CO2'
      IZN = IZ(N)
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
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101  ) THEN
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
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 ) THEN
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
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 .OR. ISCHR(IZN).EQ.5 ) THEN
          HDEN = SCHRV(1,N)
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
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKLX(1) = FNTBLY( HDGLX,IRLTBL(1,IZN),IRLTBL(2,IZN),ITBX )
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Tabular aqueous relative permeability versus ln(capillary head)
!     w/ cubic spline interpolation  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.13 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKLX(1) = FNTBLY( HDGL,IRLTBL(1,IZN),IRLTBL(2,IZN),ITBX )
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Tabular aqueous relative permeability versus capillary head
!     w/ cubic spline interpolation  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.15 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
!
!---  Mualem porosity distribution function w/ van Genuchten 
!     saturation function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.41 ) THEN
        RKLX(1) = SQRT(ESLZ)*((1.D+0-(1.D+0-ESLZ**(1.D+0/RPLC(2,IZN)))
     &    **RPLC(2,IZN))**2)
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Mualem porosity distribution function w/ Brooks and Corey 
!     saturation function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.42 ) THEN
        RKLX(1) = ESLZ**(2.5D+0 + 2.0D+0/RPLC(2,IZN))
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Burdine porosity distribution function w/ van Genuchten 
!     saturation function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.43 ) THEN
        RKLX(1) = (ESLZ**2)*(1.D+0-(1.D+0-ESLZ**(1.D+0/RPLC(2,IZN)))
     &      **RPLC(2,IZN))
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Burdine porosity distribution function w/ Brooks and Corey 
!     saturation function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.44 ) THEN
        RKLX(1) = ESLZ**(3.0D+0 + 2.0D+0/RPLC(2,IZN))
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
      ENDIF
      RKLX(1) = MIN( MAX( RKLX(1),1.D-24 ),1.D+0 )
      RKLX(2) = MIN( MAX( RKLX(2),1.D-24 ),1.D+0 )
      RKLX(3) = MIN( MAX( RKLX(3),1.D-24 ),1.D+0 )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RKLS_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Compute the maximum relative residuals
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE OUTPU
      USE LEAK_WELL
      USE JACOB
      USE HYST
      USE GRID
      USE FILES
      USE FDVS
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
      CHARACTER*64 PH_CND(10)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Saturated w/ Aqueous CO2',
     &   'Unsaturated Subcritical Pressure Gas', 
     &   'Saturated w/ Trapped Subcritical Pressure Gas',
     &   'Unsaturated Subcritical Pressure Liquid', 
     &   'Saturated w/ Trapped Subcritical Presssure Liquid',
     &   'Unsaturated Supercritical Pressure Gas/Liquid', 
     &   'Saturated w/ Trapped Supercritical Pressure Gas/Liquid',
     &   'Fully Unsaturated Subcritical Pressure Gas',
     &   'Fully Unsaturated Subcritical Pressure Liquid',
     &   'Fully Unsaturated Supercritical Pressure Gas/Liquid'/
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL_CO2'
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
      DO 200 N = 1,NFLD+NWN_LW
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) GOTO 200
        IZN = IZ(N)
        N_DB = N
        NMD = IXP(N)
        MPL = IM(IEQW,NMD)
        MPG = IM(IEQA,NMD)
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) MPS = IM(IEQS,NMD)
!
!---    Skip selected nodes in the residual calculation  ---
!
        IF( ISKP(IZN).EQ.1 ) GOTO 200
!
!---    Saturated system w/o trapped gas
!       Water mass - aqueous pressure
!       CO2 mass - aqueous-CO2 mass fraction
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
!---      CO2 mass equation, ignore residual for small aqueous-CO2  ---
!
          CALL SOL_LS( T(2,N),XLSMX )
          XLSX = MIN( YLS(2,N),XLSMX )
          XLS(2,N) = XLSX
          CALL SP_B( T(2,N),XLSX,PSBX )
          PVBX = PSBX
          PX = PL(2,N)+PATM
          XLSSX = XLS(2,N)
          CALL EQUIL( T(2,N),PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          XLS(2,N) = XLS(2,N) + (XLSSX-XLS(2,N))*(XLA(2,N)/XLAX)
          IF( XLA(2,N).GT.(1.D-6*XLAX) ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLA(2,N))*DTI*VOL(N)
            CALL HC_LA( T(2,N),XLS(2,N),HCX )
            RSDX = MIN( ABS(BLU(MPG))/
     &        MAX( (PG(2,N)+PATM)/HCX,PATM/HCX ),
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
            ACP = TMS(2,N)*DTI*VOL(N)
            CALL SOL_LS( T(2,N),XLSMX )
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
!       Water mass - aqueous pressure
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 .OR. NPHAZ(2,N).EQ.4 
     &    .OR. NPHAZ(2,N).EQ.6 ) THEN
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
!---      CO2 mass equation  ---
!
          IF( SG(2,N).GT.1.D-3 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLA(2,N))*DTI*VOL(N)
            RSDX = MIN( ABS(BLU(MPG))/(ABS(PG(2,N))+PATM),
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
            ACP = TMS(2,N)*DTI*VOL(N)
            CALL SOL_LS( T(2,N),XLSMX )
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
!       Water mass - aqueous pressure
!       CO2 mass - trapped gas saturation
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.5
     &    .OR. NPHAZ(2,N).EQ.7 ) THEN
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
!---      CO2 mass equation  ---
!
          IF( SG(2,N).GT.1.D-5 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLA(2,N))*DTI*VOL(N)
            RSDX = MIN( ABS(BLU(MPG)/1.D+1),
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
            ACP = TMS(2,N)*DTI*VOL(N)
            CALL SOL_LS( T(2,N),XLSMX )
            RSDX = MIN( (ABS(BLU(MPS))/XLSMX),
     &        ABS(RSDL(IEQS,N)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            IF( RSDX.GT.RSD(IEQS) ) THEN
              RSD(IEQS) = RSDX
              NSD(IEQS) = N
            ENDIF
          ENDIF
!
!---    Fully unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous saturation
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.8 .OR. NPHAZ(2,N).EQ.9 
     &    .OR. NPHAZ(2,N).EQ.10 ) THEN
          PGX = PG(2,N) + PATM
          CALL SP_B( T(2,N),XLS(2,N),PSBX )
          ACP = PORD(2,N)*RHOG(2,N)*SG(2,N)*XGW(2,N)*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPL))/PSBX,
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &      RHOL(2,N)*SL(2,N)*XLA(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPG))/ABS(PGX),
     &      ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQA) ) THEN
            RSD(IEQA) = RSDX
            NSD(IEQA) = N
          ENDIF
!
!---      Salt mass equation, isobrine option  ---

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
  200 CONTINUE
!
!---  Assign a convergence index  ---
!
      RSDX = 1.D-20
      DO 300 M = 1,ISVC
        IF( RSD(M).GT.RSDMX ) ICNV = 2
        RSDX = MAX( RSD(M),RSDX )
  300 CONTINUE
      IF( ICNV.EQ.2 .AND. NITER.GE.NRIMX ) ICNV = 1
!
!---  Unconverged solution Newton-Raphson iteration limit exceeded  ---
!
      IF( ICNV.EQ.1 ) THEN
        IF( RSDX.GE.1.D+2 .AND. RSD_CW.LT.1.D+2 ) THEN
          WRITE(ISC,'(10X,A)') '---  Excessive Residual  ---'
          WRITE(IWR,'(10X,A)') '---  Excessive Residual  ---'
        ELSEIF( RSD_CW.LT.RSDMX ) THEN
          WRITE(ISC,'(10X,A)') '---  Convergence Failure  ---'
          WRITE(IWR,'(10X,A)') '---  Convergence Failure  ---'
        ENDIF
        NW = NSD(IEQW)
        NA = NSD(IEQA)
        IF( NW.GT.0 ) THEN
          NPW = NPHAZ(2,NW)
          NCHW = INDEX( PH_CND(NPW),'  ') - 1
          IF( NW.GT.NFLD .AND. L_LW.EQ.1 ) THEN
            NWX = NW-NFLD
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Water Equation Maximum Residual = ',RSD(IEQW),
     &        ': Leaky Well Node = ',NWX,
     &        ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Water Equation Maximum Residual = ',RSD(IEQW),
     &        ': Leaky Well Node = ',NWX,
     &        ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
          ELSE
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &        'Water Equation Maximum Residual = ',RSD(IEQW),
     &        ': Node = ',NW,': I,J,K = ',ID(NW),',',JD(NW),',',KD(NW),
     &        ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &        'Water Equation Maximum Residual = ',RSD(IEQW),
     &        ': Node = ',NW,': I,J,K = ',ID(NW),',',JD(NW),',',KD(NW),
     &        ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
          ENDIF
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
        IF( NA.GT.0 ) THEN
          NPA = NPHAZ(2,NA)
          NCHA = INDEX( PH_CND(NPA),'  ') - 1
          IF( NA.GT.NFLD .AND. L_LW.EQ.1 ) THEN
            NWX = NA-NFLD
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'CO2 Equation Maximum Residual = ',RSD(IEQA),
     &        ': Leaky Well Node = ',NWX,
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'CO2 Equation Maximum Residual = ',RSD(IEQA),
     &        ': Leaky Well Node = ',NWX,
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
          ELSE
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &        'CO2 Equation Maximum Residual = ',RSD(IEQA),
     &        ': Node = ',NA,': I,J,K = ',ID(NA),',',JD(NA),',',KD(NA),
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &        'CO2 Equation Maximum Residual = ',RSD(IEQA),
     &        ': Node = ',NA,': I,J,K = ',ID(NA),',',JD(NA),',',KD(NA),
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
          ENDIF
!
!---      Extended output on convergence failure  ---
!
          IF( ISLC(62).EQ.1 ) THEN
            WRITE(ISC,'(4X,2A)')
     &        'Current Phase Condition = ',PH_CND(NPA)(1:NCHA)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Current Aqueous Saturation = ',SL(2,NA)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Current Trapped Gas Saturation = ',SGT(2,NA)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Current Salt Saturation = ',SS(2,NA)
            WRITE(IWR,'(4X,2A)')
     &        'Current Phase Condition = ',PH_CND(NPA)(1:NCHA)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Current Aqueous Saturation = ',SL(2,NA)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Current Trapped Gas Saturation = ',SGT(2,NA)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Current Salt Saturation = ',SS(2,NA)
            NPA = NPHAZ(1,NA)
            NCHA = INDEX( PH_CND(NPA),'  ') - 1
            WRITE(ISC,'(4X,2A)')
     &        'Previous Phase Condition = ',PH_CND(NPA)(1:NCHA)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Previous Aqueous Saturation = ',SL(1,NA)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Previous Trapped Gas Saturation = ',SGT(1,NA)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Previous Salt Saturation = ',SS(1,NA)
            WRITE(IWR,'(4X,2A)')
     &        'Previous Phase Condition = ',PH_CND(NPA)(1:NCHA)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Previous Aqueous Saturation = ',SL(1,NA)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Previous Trapped Gas Saturation = ',SGT(1,NA)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Previous Salt Saturation = ',SS(1,NA)
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
            IF( NS.GT.NFLD .AND. L_LW.EQ.1 ) THEN
              NWX = NS-NFLD
              WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &          'Salt Equation Maximum Residual = ',RSD(IEQS),
     &          ': Leaky Well Node = ',NWX,
     &          ': Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &          'Salt Equation Maximum Residual = ',RSD(IEQS),
     &          ': Leaky Well Node = ',NWX,
     &          ': Phase Condition = ',PH_CND(NPS)(1:NCHS)
            ELSE
              WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &          'Salt Equation Maximum Residual = ',RSD(IEQS),
     &          ': Node = ',NS,': I,J,K = ',ID(NS),',',JD(NS),',',
     &          KD(NS),': Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &          'Salt Equation Maximum Residual = ',RSD(IEQS),
     &          ': Node = ',NS,': I,J,K = ',ID(NS),',',JD(NS),',',
     &          KD(NS),': Phase Condition = ',PH_CND(NPS)(1:NCHS)
            ENDIF
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
          DO 400 N = 1,NFLD+NWN_LW
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PVW(2,N) = PVW(1,N)
            XLA(2,N) = XLA(1,N)
            SG(2,N) = SG(1,N)
            SGT(2,N) = SGT(1,N)
            ASLMIN(2,N) = ASLMIN(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  400     CONTINUE
!
!---      Coupled-well pressure  ---
!
          IF( L_CW.EQ.1 ) THEN
            DO 402 NCW = 1,N_CW
              P_CW(2,NCW) = P_CW(1,NCW)
  402       CONTINUE
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDL_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SBND_CO2( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Modify the Jacobian matrix for the solute transport equation
!     to incorporate boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
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
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 BCX(LSPBC+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SBND_CO2'
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
        ELSEIF( ILES.EQ.3 ) THEN
          MA = 1
          MCOL = KLUC(MP,MA)
          MA = MA + 1
#ifdef lis
        ELSEIF( ILES.EQ.4 ) THEN
          MA = 1
          MCOL = KLUC(MP,MA)
          MA = MA + 1
#endif
#ifdef petsc
        ELSEIF( ILES.EQ.5 ) THEN
          MA = 1
          MCOL = KLUC(MP,MA)
          MA = MA + 1
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
!---        Linked aqueous CO2   ---
!
            IF( ISPLK(6).EQ.NSL ) BCX(1) = 1.D+3*XLAB(2,NB)*
     &          RHOLB(2,N)*SLB(2,NB)*PORDB(2,NB)/WTMA
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
!---        Linked aqueous CO2   ---
!
            IF( ISPLK(6).EQ.NSL ) BCX(1) = 1.D+3*XLAB(2,NB)*
     &          RHOLB(2,N)*SLB(2,NB)*PORDB(2,NB)/WTMA
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
              IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
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
              IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
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
              IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
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
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
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
              IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
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
              IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
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
              IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
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
              IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
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
              IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
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
              IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
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
          IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
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
          IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
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
          IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
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
        ELSEIF( ILES.EQ.3 ) THEN
          DLU(MCOL) = DLU(MCOL) + AP
#ifdef lis
        ELSEIF( ILES.EQ.4 ) THEN
          DLU(MCOL) = DLU(MCOL) + AP
#endif
#ifdef petsc
        ELSEIF( ILES.EQ.5 ) THEN
          DLU(MCOL) = DLU(MCOL) + AP
#endif
        ENDIF
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SBND_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORC_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Compute source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE REACT
      USE LEAK_WELL
      USE JACOB
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
      REAL*8 SRX(8+LSOLU)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORC_CO2'
#ifdef ecke
!
!---  Load CO2 sources for linked aqueous CO2   ---
!
      IF( ISPLK(6).NE.0 ) THEN
!
!---    Load CO2 sources associated with chemical
!       reactions  ---
!
        DO 20 N = 1,NFLD+NWN_LW
          SRCAX = SRCA(1,N)*DTI
          DO 10 M = 2,ISVC+2
            SRCA(M,N) = SRCAX
   10     CONTINUE
   20   CONTINUE
!
!---    Zero source terms  ---
!
        IF( NSR.GT.0 ) THEN
          DO 40 N = 1,NFLD+NWN_LW
            DO 30 M = 2,ISVC+2
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
          DO 60 N = 1,NFLD+NWN_LW
            DO 50 M = 2,ISVC+2
              SRCA(M,N) = 0.D+0
              SRCW(M,N) = 0.D+0
              SRCS(M,N) = 0.D+0
   50       CONTINUE
   60     CONTINUE
        ENDIF
#ifdef ecke
      ENDIF
#endif
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
!---    Loop source domain  ---
!
        DO I = ISRDM(1,NS),ISRDM(2,NS)
        DO J = ISRDM(3,NS),ISRDM(4,NS)
        DO K = ISRDM(5,NS),ISRDM(6,NS)
          N = ND(I,J,K)
          N_DB = N
          IZN = IZ(N)
!
!---      Skip inactive nodes  ---
!
          IF( IXP(N).EQ.0 ) CYCLE
          DO 400 M = 2,ISVC+2
            PGX = PG(M,N) + PATM
            PLX = PL(M,N) + PATM
            TX = T(M,N)
!
!---        Aqueous volumetric rate w/ aqueous-salt
!           and aqueous-CO2; limit aqueous components
!           to their solubility limit  ---
!
            IF( MOD(ISRT(NS),100).EQ.3 ) THEN
              IF( SRX(4).LT.0.D+0 ) THEN
                SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOL(M,N)*XLA(M,N)
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOL(M,N)*XLW(M,N)
                SRCS(M,N) = SRCS(M,N) + SRX(4)*RHOL(M,N)*XLS(M,N)
              ELSE
                PORDX = 1.D+0
                SLX = 1.D+0
                PX = MAX( PLX,PGX,SRX(3) )
                PCX = 0.D+0
!
!---            Aqueous-salt concentration  ---
!
                IF( MOD(ISRT(NS),1000)/100.EQ.1 ) THEN
!
!---              Aqueous-CO2 concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    RHOLSX = SRX(8)
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Conc., Salt Aqu. Conc. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_11( TX,PX,PGX,RHOLSX,RHOLAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    RHOLSX = SRX(8)
                    PHILAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Rel. Sat., Salt Aqu. Conc. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_12( TX,PX,PGX,RHOLSX,PHILAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    RHOLSX = SRX(8)
                    XLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Conc. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_13( TX,PX,PGX,RHOLSX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              No aqueous-CO2  ---
!
                  ELSE
                    RHOLSX = SRX(8)
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Conc. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_13( TX,PX,PGX,RHOLSX,
     &                YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
!
!---            Aqueous-salt relative saturation  ---
!
                ELSEIF( MOD(ISRT(NS),1000)/100.EQ.2 ) THEN
!
!---              Aqueous-CO2 concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    PHILSX = SRX(8)
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Salt Aqu. Rel. Conc., CO2 Aqu. Conc.  = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_21( TX,PX,PGX,PHILSX,RHOLAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    PHILSX = SRX(8)
                    PHILAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'Salt Aqu. Rel. Sat., CO2 Aqu. Rel. Sat. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_22( TX,PX,PGX,PHILSX,PHILAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    PHILSX = SRX(8)
                    XLAX = SRX(5)
                    CHMSGX(1) = 'null'
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_23( TX,PX,PGX,PHILSX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              No aqueous-CO2  ---
!
                  ELSE
                    PHILSX = SRX(8)
                    XLAX = 0.D+0
                    CHMSGX(1) = 'null'
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_23( TX,PX,PGX,PHILSX,
     &                YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
!
!---            Aqueous-salt mass fraction  ---
!
                ELSEIF( MOD(ISRT(NS),1000)/100.EQ.3 ) THEN
!
!---              Aqueous-CO2 concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    YLSX = SRX(8)
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Conc., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_31( TX,PX,PGX,RHOLAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    YLSX = SRX(8)
                    PHILAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Rel. Sat., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_32( TX,PX,PGX,PHILAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    YLSX = SRX(8)
                    XLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---              No aqueous-CO2  ---
!
                  ELSE
                    YLSX = SRX(8)
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
!
!---            Aqueous-salt molality  ---
!
                ELSEIF( MOD(ISRT(NS),1000)/100.EQ.4 ) THEN
!
!---              Aqueous-CO2 concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    YLSX = SRX(8)
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Conc., Salt Aqu. Molality = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_41( TX,PX,PGX,RHOLAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    YLSX = SRX(8)
                    PHILAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Rel. Sat., Salt Aqu. Molality = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_42( TX,PX,PGX,PHILAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    YLSX = SRX(8)
                    XLAX = SRX(5)
                    CHMSGX(1) = 'null'
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_43( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---              No aqueous-CO2  ---
!
                  ELSE
                    YLSX = SRX(8)
                    XLAX = 0.D+0
                    CHMSGX(1) = 'null'
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_43( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
!
!---            No aqueous-salt  ---
!
                ELSE
!
!---              Aqueous-CO2 concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    YLSX = 0.D+0
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Conc., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_31( TX,PX,PGX,RHOLAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    YLSX = 0.D+0
                    PHILAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Rel. Sat., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_32( TX,PX,PGX,PHILAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    YLSX = 0.D+0
                    XLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---              No aqueous-CO2  ---
!
                  ELSE
                    YLSX = 0.D+0
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Unconverged Source: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &                'Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
                ENDIF
!
!---            Aqueous component fractions and density  ---
!
                XLWX = 1.D+0 - XLSX - XLAX
                CALL DENS_B( TX,PX,XLSX,RHOBX )
                CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
                SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOLX*XLAX
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOLX*XLWX
                SRCS(M,N) = SRCS(M,N) + SRX(4)*RHOLX*XLSX
              ENDIF
!
!---        Gas volumetric rate w/ component relative humidity ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.4 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                PX = MAX( PGX,PLX,SRX(3) )
                CALL SP_W( TX,PSWX )
                PVWX = PSWX*SRX(5)
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
                CALL DENS_A( TX,PVAX,RHOGAX,I_VX )
                SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOGAX
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOGWX
              ELSE
                SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOG(M,N)*XGA(M,N)
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOG(M,N)*XGW(M,N)
              ENDIF
!
!---        Gas volumetric rate w/ component mass fractions ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.5 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                PX = MAX( PGX,PLX,SRX(3) )
                XGWX = SRX(5)
                XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                WTMGX = 1.D+0/(XGAX/WTMA + XGWX/WTMW)
                XMGAX = XGAX*WTMGX/WTMA
                XMGWX = XGWX*WTMGX/WTMW
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
                CALL DENS_A( TX,PVAX,RHOGAX,I_VX )
                SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOGAX
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOGWX
              ELSE
                SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOG(M,N)*XGA(M,N)
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOG(M,N)*XGW(M,N)
              ENDIF
!
!---        Salt source  ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.12 ) THEN
              SRCS(M,N) = SRCS(M,N) + SRX(4)
!
!---        Aqueous mass rate w/ aqueous-CO2 relative
!           saturation, limit aqueous salt
!           sources to solubility limit  ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.7 ) THEN
              IF( SRX(4).LT.0.D+0 ) THEN
                SRCA(M,N) = SRCA(M,N) + SRX(4)*XLA(M,N)
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XLW(M,N)
                SRCS(M,N) = SRCS(M,N) + SRX(4)*XLS(M,N)
              ELSE
                PORDX = 1.D+0
                SLX = 1.D+0
                PX = MAX( PLX,PGX,SRX(3) )
                PCX = 0.D+0
!
!---            Aqueous-salt concentration  ---
!
                IF( MOD(ISRT(NS),1000)/100.EQ.1 ) THEN
!
!---              Aqueous-CO2 concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    RHOLSX = SRX(8)
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Conc., Salt Aqu. Conc. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_11( TX,PX,PGX,RHOLSX,RHOLAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    RHOLSX = SRX(8)
                    PHILAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Rel. Sat., Salt Aqu. Conc. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_12( TX,PX,PGX,RHOLSX,PHILAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    RHOLSX = SRX(8)
                    XLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Conc. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_13( TX,PX,PGX,RHOLSX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              No aqueous-CO2  ---
!
                  ELSE
                    RHOLSX = SRX(8)
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Conc. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_13( TX,PX,PGX,RHOLSX,
     &                YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
!
!---            Aqueous-salt relative saturation  ---
!
                ELSEIF( MOD(ISRT(NS),1000)/100.EQ.2 ) THEN
!
!---              Aqueous-CO2 concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    PHILSX = SRX(8)
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'Salt Aqu. Rel. Conc., CO2 Aqu. Conc.  = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_21( TX,PX,PGX,PHILSX,RHOLAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    PHILSX = SRX(8)
                    PHILAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'Salt Aqu. Rel. Sat., CO2 Aqu. Rel. Sat. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_22( TX,PX,PGX,PHILSX,PHILAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    PHILSX = SRX(8)
                    XLAX = SRX(5)
                    CHMSGX(1) = 'null'
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_23( TX,PX,PGX,PHILSX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              No aqueous-CO2  ---
!
                  ELSE
                    PHILSX = SRX(8)
                    XLAX = 0.D+0
                    CHMSGX(1) = 'null'
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_23( TX,PX,PGX,PHILSX,
     &                YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
!
!---            Aqueous-salt mass fraction  ---
!
                ELSEIF( MOD(ISRT(NS),1000)/100.EQ.3 ) THEN
!
!---              Aqueous-CO2 concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    YLSX = SRX(8)
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Conc., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_31( TX,PX,PGX,RHOLAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    YLSX = SRX(8)
                    PHILAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Rel. Sat., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_32( TX,PX,PGX,PHILAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    YLSX = SRX(8)
                    XLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---              No aqueous-CO2  ---
!
                  ELSE
                    YLSX = SRX(8)
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
!
!---            Aqueous-salt mass fraction  ---
!
                ELSEIF( MOD(ISRT(NS),1000)/100.EQ.4 ) THEN
!
!---              Aqueous-CO2 concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    YLSX = SRX(8)
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Conc., Salt Aqu. Molality = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_41( TX,PX,PGX,RHOLAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    YLSX = SRX(8)
                    PHILAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Rel. Sat., Salt Aqu. Molality = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_42( TX,PX,PGX,PHILAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    YLSX = SRX(8)
                    XLAX = SRX(5)
                    CHMSGX(1) = 'null'
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_43( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---              No aqueous-CO2  ---
!
                  ELSE
                    YLSX = SRX(8)
                    XLAX = 0.D+0
                    CHMSGX(1) = 'null'
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_43( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
!
!---            No aqueous-salt  ---
!
                ELSE
!
!---              Aqueous-CO2 concentration  ---
!
                  IF( ISRT(NS)/1000.EQ.1 ) THEN
                    YLSX = 0.D+0
                    RHOLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Conc., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_31( TX,PX,PGX,RHOLAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 relative saturation  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                    YLSX = 0.D+0
                    PHILAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Rel. Sat., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_32( TX,PX,PGX,PHILAX,
     &                YLSX,XLSX,XLAX,CHMSGX )
!
!---              Aqueous-CO2 mass fraction  ---
!
                  ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                    YLSX = 0.D+0
                    XLAX = SRX(5)
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---              No aqueous-CO2  ---
!
                  ELSE
                    YLSX = 0.D+0
                    XLAX = 0.D+0
                    CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &                'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                    CHMSGX(2) = 'Initial Condition Transition: ' //
     &                'Vapor Pressure > Gas Pressure: PVB + PVA = '
                    CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                  ENDIF
                ENDIF
!
!---            CO2 vapor pressure ---
!
                CALL SP_B( TX,XLSX,PSBX )
                PVBX = PSBX
                CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &            XGAX,XGWX,XLAZ,XLSX,XLWX,XMGAX,XMGWX,XMLAX,
     &            XMLSX,XMLWX )
                PVAX = MIN( XLAX/XLAZ,1.D+0 )*PGAX
!
!---            Aqueous component fraction  ---
!
                XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
                SRCA(M,N) = SRCA(M,N) + SRX(4)*XLAX
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XLWX
                SRCS(M,N) = SRCS(M,N) + SRX(4)*XLSX
              ENDIF
!
!---        Gas mass rate w/ component relative humidity ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.8 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                PX = MAX( PGX,PLX,SRX(3) )
                CALL SP_W( TX,PSWX )
                PVWX = PSWX*SRX(5)
                PVAX = MAX( PX-PVWX,0.D+0 )
!
!---            Gas density and component fractions  ---
!
                IF( PVWX.GT.EPSL ) THEN
                  ISRX = 2
                  CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
                ELSE
                  RHOGWX = 0.D+0
                ENDIF
                CALL DENS_A( TX,PVAX,RHOGAX,I_VX )
                XGAX = RHOGAX/(RHOGAX+RHOGWX)
                XGWX = RHOGWX/(RHOGAX+RHOGWX)
                SRCA(M,N) = SRCA(M,N) + SRX(4)*XGAX
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGWX
              ELSE
                SRCA(M,N) = SRCA(M,N) + SRX(4)*XGA(M,N)
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGW(M,N)
              ENDIF
!
!---        Gas mass rate w/ component mass fractions ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.9 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                PX = MAX( PGX,PLX,SRX(3) )
                XGWX = SRX(5)
                XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                WTMGX = 1.D+0/(XGAX/WTMA + XGWX/WTMW)
                XMGAX = XGAX*WTMGX/WTMA
                XMGWX = XGWX*WTMGX/WTMW
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
                CALL DENS_A( TX,PVAX,RHOGAX,I_VX )
                SRCA(M,N) = SRCA(M,N) + SRX(4)*XGAX
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGWX
              ELSE
                SRCA(M,N) = SRCA(M,N) + SRX(4)*XGA(M,N)
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGW(M,N)
              ENDIF
!
!---        Salt density source  ---
!
            ELSEIF( MOD(ISRT(NS),100).EQ.11 ) THEN
              SRCS(M,N) = SRCS(M,N) + SRX(4)*VOL(N)
            ENDIF
  400     CONTINUE
        ENDDO
        ENDDO
        ENDDO
  900 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORC_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORIC_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Compute coupled-equation source integrals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 October 2007
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE GRID
      USE COUP_WELL
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
      SUB_LOG(ISUB_LOG) = '/SORIC_CO2'
!
!---  Integrate coupled-well rates  ---
!
      IF( L_CW.EQ.1 ) THEN
        DO 10 NCW = 1,N_CW
          QM_CW(2,NCW) = QM_CW(2,NCW) + QM_CW(1,NCW)*DT
          QM_CW(4,NCW) = QM_CW(4,NCW) + QM_CW(3,NCW)*DT
   10   CONTINUE
      ENDIF
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
          DO 200 NC = 1,ISRDM(1,NS)
            N = IWSI(NC,NS)
            IF( IXP(N).EQ.0 ) GOTO 200
            SRCP(1,NS) = SRCP(1,NS) + SWSI(1,NC,NS)*SWSI(2,NC,NS)*DT
            SRCP(2,NS) = SRCP(2,NS) + SWSI(1,NC,NS)*SWSI(3,NC,NS)*DT
  200     CONTINUE
        ENDIF
  900 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORIC_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORIT_CO2( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Compute solute transport source integrals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      REAL*8 SRX(8+LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORIT_CO2'
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
        DO I = ISRDM(1,NS),ISRDM(2,NS)
          DO J = ISRDM(3,NS),ISRDM(4,NS)
            DO K = ISRDM(5,NS),ISRDM(6,NS)
              N = ND(I,J,K)
!
!---          Skip inactive nodes  ---
!
              IF( IXP(N).EQ.0 ) CYCLE
!
!---          Aqueous Volumetric Sink  ---
!
              IF( MOD(ISRT(NS),100).EQ.3 .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YL(N,NSL)*DT/(PORD(2,N)*SL(2,N))
!
!---          Gas Volumetric Sink  ---
!
              ELSEIF( (MOD(ISRT(NS),100).EQ.4 .OR. 
     &          MOD(ISRT(NS),100).EQ.5) .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YG(N,NSL)*DT/(PORD(2,N)*SG(2,N))
!
!---          Aqueous Mass Sink  ---
!
              ELSEIF( MOD(ISRT(NS),100).EQ.7 .AND. 
     &          SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YL(N,NSL)*DT/(RHOL(2,N)*PORD(2,N)*SL(2,N))
!
!---          Gas Mass Sink  ---
!
              ELSEIF( (MOD(ISRT(NS),100).EQ.8 .OR. 
     &          MOD(ISRT(NS),100).EQ.9) .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YG(N,NSL)*DT/(RHOG(2,N)*PORD(2,N)*SG(2,N))
!
!---          Solute source  ---
!
              ELSEIF( ISRT(NS).LT.0 .AND. ISRT(NS).GE.-NSOLU ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) + SRX(4)*DT
!
!---         Solute source  ---
!
              ELSEIF( ISRT(NS).LT.-NSOLU .AND.
     &          ISRT(NS).GE.-2*NSOLU ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) + SRX(4)*DT*VOL(N)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORIT_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT_CO2( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Compute solute transport source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
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
      SUB_LOG(ISUB_LOG) = '/SORT_CO2'
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
             DO 80 N = 2,8+NSOLU
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
        DO I = ISRDM(1,NS),ISRDM(2,NS)
        DO J = ISRDM(3,NS),ISRDM(4,NS)
        DO K = ISRDM(5,NS),ISRDM(6,NS)
          N = ND(I,J,K)
!
!---      Skip inactive nodes  ---
!
          IF( IXP(N).EQ.0 ) CYCLE
          MP = IXP(N)
          IF( ILES.EQ.1 ) THEN
            MCOL = MP
            MROW = MDT
          ELSEIF( ILES.EQ.3 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#ifdef lis
          ELSEIF( ILES.EQ.4 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#endif
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
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
          ELSEIF( ILES.EQ.3 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#ifdef lis
          ELSEIF( ILES.EQ.4 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#endif
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#endif
          ENDIF
        ENDDO
        ENDDO
        ENDDO
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORT_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SP_CO2( ASL_F,ASL_M,ASLX,ASLMINX,BTGLX,ESLX,ESGTX,
     &  ESGTMX,PGX,PLX,SGX,SGTX,SLX,SLFX,SLMX,SLRX,INDX,N )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
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
      SUB_LOG(ISUB_LOG) = '/SP_CO2'
      IZN = IZ(N)
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
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
!zfz
!            ESLXX = (1.D+0 + (SCHRV(1,N)*HDGL)**CN)**(-CM)
!            SLRX = 1.0D+0 -
!     &        (1.0D+0 + (LOG10(SCHR(12,IZN))-LOG10(HDGL))*DMPX)/
!     &        (1.0D+0 - ESLXX)
!zfz
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
            ASLX = (1.D+0/(1.D+0 + (SCHRV(1,N)*HDGL)**CN))**CM
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
          ASLX = (1.D+0/(1.D+0 + (SCHRV(1,N)*HDGL)**CN))**CM
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
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
!zfz
!            CL = MAX( SCHR(3,IZN),SMALL )
!            ESLXX = (SCHRV(1,N)/HDGL)**(CL)
!            SLRX = 1.0D+0 -
!     &        (1.0D+0 + (LOG10(SCHR(12,IZN))-LOG10(HDGL))*DMPX)/
!     &        (1.0D+0 - ESLXX)
!zfz
!
!---      Capillary head at or below the matching point head,
!         use Brooks-Corey function
!
          ELSE
            CL = MAX( SCHR(3,IZN),SMALL )
            IF( HDGL.LE.SCHRV(1,N) ) THEN
              ASLX = 1.D+0
            ELSE
              ASLX = (SCHRV(1,N)/HDGL)**CL
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
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
          CL = MAX( SCHR(3,IZN),SMALL )
          IF( HDGL.LE.SCHRV(1,N) ) THEN
            ASLX = 1.D+0
          ELSE
            ASLX = (SCHRV(1,N)/HDGL)**CL
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
        ASL_M = (1.D+0/(1.D+0 + (SCHRV(1,N)*HDGL)**CN))**CM
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,IZN)),ZERO )*REALX
        SLR_M = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        SLMX = ASL_M*(1.D+0-SLR_M) + SLR_M
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
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,IZN)),ZERO )*REALX
        SLR_F = MAX( (1.D+0-HSCL)*SCHR(7,IZN),ZERO )
        SLFX = ASL_F*(1.D+0-SLR_F) + SLR_F
        PORD_MX = (1.D+0-POR(4,IZN))*POR0(2,N)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR0(2,N) + SMALL )
        PORD_FX = POR(4,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR0(2,N) + SMALL )
        SLX = SLFX*PORD_FX + SLMX*PORD_MX
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ASLX = ASL_F*PORD_FX + ASL_M*PORD_MX
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        ESGTX = 0.D+0
        SGTX = 0.D+0
!
!---  Dual porosity Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.4 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        CL = MAX( SCHR(3,IZN),SMALL )
        IF( HDGL.LE.SCHRV(1,N) ) THEN
          ASL_M = 1.D+0
        ELSE
          ASL_M = (SCHRV(1,N)/HDGL)**CL
        ENDIF
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,IZN)),ZERO )*REALX
        SLR_M = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        SLMX = ASL_M*(1.D+0-SLR_M) + SLR_M
        CL = MAX( SCHR(6,IZN),SMALL )
        IF( HDGL.LE.SCHR(5,IZN) ) THEN
          ASL_F = 1.D+0
        ELSE
          ASL_F = (SCHR(5,IZN)/HDGL)**CL
        ENDIF
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,IZN)),ZERO )*REALX
        SLR_F = MAX( (1.D+0-HSCL)*SCHR(7,IZN),ZERO )
        SLFX = ASL_F*(1.D+0-SLR_F) + SLR_F
        PORD_MX = (1.D+0-POR(4,IZN))*POR0(2,N)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR0(2,N) + SMALL )
        PORD_FX = POR(4,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR0(2,N) + SMALL )
        SLX = SLFX*PORD_FX + SLMX*PORD_MX
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ASLX = ASL_F*PORD_FX + ASL_M*PORD_MX
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        ESGTX = 0.D+0
        SGTX = 0.D+0
!
!---  Haverkamp saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.5 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( HDGL.LE.SCHRV(1,N) ) THEN
          ASLX = 1.D+0
        ELSE
          ASLX = SCHR(2,IZN)/(SCHR(2,IZN)
     &      + (((HDGL-SCHRV(1,N))/SCHR(5,IZN))**SCHR(3,IZN)))
        ENDIF
        REALX = REAL(ISM(IZN))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,IZN)),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,IZN),ZERO )
        SLX = ASLX*(1.D+0-SLRX) + SLRX
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        ESGTX = 0.D+0
        SGTX = 0.D+0
!
!---  Linear or linear-log interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.10 .OR. ISCHR(IZN).EQ.11 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
!---    Compute trapped gas saturation, using the minimum
!       apparent aqueous saturation and the drainage saturation  ---
!
        IF( INDX.EQ.0 ) THEN
          ESGTX = SLDX-ASLX
          IF( ESGTX.LT.EPSL ) ESGTX = 0.D+0
          SGTX = ESGTX*(1.D+0-SLRX)
!
!---    Trapped-gas saturation given, compute the minimum
!       apparent aqueous saturation  ---
!
        ELSEIF( INDX.EQ.1 ) THEN
          IF( (SLDX-SLIX-ESGTX).GT.0.D+0 ) THEN
            ASLM = (SLDX-SLIX-ESGTX)*(ASLX/(SLDX-SLIX))
            ASLMINX = MAX( MIN( ASLM,ASLX ),0.D+0 )
          ELSE
            ALSMINX = 0.D+0
          ENDIF
        ENDIF
        SLX = ASLX
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
!!
!!---  Cubic-spline or cubic-spline-log interpolation function  ---
!!
!      ELSEIF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.13 ) THEN
!        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
!        IF( ISCHR(IZN).EQ.13 ) HDGL = LOG(HDGL)
!        ASLX = FSPLNY( HDGL,ISLTBL(1,IZN),ISLTBL(2,IZN) )
!        ESLX = ASLX
!        SLX = ASLX
!        SGX = 1.D+0-SLX
!        IF( SGX.LT.EPSL ) SGX = 0.D+0
!        ESGTX = 0.D+0
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
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
!zfz
!            ESLXX = (1.D+0 + (SCHRV(1,N)*HDGL)**CN)**(-CM)
!            SLRX = 1.0D+0 -
!     &        (1.0D+0 + (LOG10(SCHR(12,IZN))-LOG10(HDGL))*DMPX)/
!     &        (1.0D+0 - ESLXX)
!zfz
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
            ASLX = (1.D+0/(1.D+0 + (SCHRV(1,N)*HDGL)**CN))**CM
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
          ASLX = (1.D+0/(1.D+0 + (SCHRV(1,N)*HDGL)**CN))**CM
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
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
!zfz
!            CL = MAX( SCHR(3,IZN),SMALL )
!            ESLXX = (SCHRV(1,N)/HDGL)**(CL)
!            SLRX = 1.0D+0 -
!     &        (1.0D+0 + (LOG10(SCHR(12,IZN))-LOG10(HDGL))*DMPX)/
!     &        (1.0D+0 - ESLXX)
!zfz
!
!---      Capillary head at or below the matching point head,
!         use van Genuchten function
!
          ELSE
            CL = MAX( SCHR(3,IZN),SMALL )
            IF( HDGL.LE.SCHRV(1,N) ) THEN
              ASLX = 1.D+0
            ELSE
              ASLX = (SCHRV(1,N)/HDGL)**CL
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
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
          CL = MAX( SCHR(3,IZN),SMALL )
          IF( HDGL.LE.SCHRV(1,N) ) THEN
            ASLX = 1.D+0
          ELSE
            ASLX = (SCHRV(1,N)/HDGL)**CL
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
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
            ASLDX = (1.D+0/(1.D+0 + (SCHRV(1,N)*HDGL)**CND))**CMD
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
          ASLDX = (1.D+0/(1.D+0 + (SCHRV(1,N)*HDGL)**CND))**CMD
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
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
            IF( HDGL.LE.SCHRV(1,N) ) THEN
              ASLDX = 1.D+0
            ELSE
              ASLDX = (SCHRV(1,N)/HDGL)**CLD
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
          IF( HDGL.LE.SCHRV(1,N) ) THEN
            ASLDX = 1.D+0
          ELSE
            ASLDX = (SCHRV(1,N)/HDGL)**CLD
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SP_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPRP_CO2( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Calculates the aqueous- and gas-phase solute
!     mole fractions from user-specified partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE LEAK_WELL
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
      SUB_LOG(ISUB_LOG) = '/SPRP_CO2'
!
!---  Loop over all nodes  ---
!
      DO 900 N = 1,NFLD+NWN_LW
        N_DB = N
!
!---    Skip inactive nodes  ---
!
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
     &      + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &      + PCGL(5,NSL)*TK**2 )
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
  900 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SPRP_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TPORT_CO2( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Solute/Reactive Species Transport Shell.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 16 August 2005
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
      USE SOLTN
      USE JACOB
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TPORT_CO2'
!
!---  Zero Jacobian matrix  ---
!
#ifdef petsc
      IF( ILES.EQ.5 ) CALL MatZeroEntries( T_MAT,IERR )
#endif
      INDX = 1
      CALL JCBZ( ISVT,MUT,MLT,MKT,INDX )
!
!---  Compute solute sources ---
!
      CALL SORT_CO2( NSL )
!
!---  Compute solute sources from injection and production wells ---
!
      CALL SORT_COUP_WELL( NSL )
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
      CALL SBND_CO2( NSL )
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
        CALL STOMP_LIS_SOLVE( -1,T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC,INDX )
#endif
#ifdef petsc
      ELSEIF( ILES.EQ.5 ) THEN
        INDX = 1
        CALL STOMP_PETSC_SOLVE( T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC,INDX )
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
      CALL SORIT_CO2( NSL )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TPORT_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMP-CO2
!
!     Update the primary variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE OUTPU
      USE LEAK_WELL
      USE JACOB
      USE HYST
      USE GRID
      USE FILES
      USE FDVS
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
      CHARACTER*64 PH_CND(10)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Saturated w/o Trapped Gas',
     &  'Unsaturated w/ or w/o Entrapped Gas (Subcritical Gas)',
     &  'Saturated w/ Trapped Gas (Subcritical Gas)',
     &  'Unsaturated w/ or w/o Entrapped Gas (Subcritical Liquid)',
     &  'Saturated w/ Trapped Gas (Subcritical Liquid)',
     &  'Unsaturated w/ or w/o Entrapped Gas (Supercritical Gas)',
     &  'Saturated w/ Trapped Gas (Supercritical Gas)',
     &  'Fully Unsaturated (Subcritical Gas)',
     &  'Fully Unsaturated (Subcritical Liquid)',
     &  'Fully Unsaturated (Supercritical Gas)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_CO2'
      IF( ICNV.EQ.1 ) GOTO 300
      IERR = 0
      NSDX = NFLD
!
!---  Update primary variables
!
      DO 200 N = 1,NFLD+NWN_LW
        N_DB = N
        IZN = IZ(N)
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) GOTO 200
        NMD = IXP(N)
        MPL = IM(IEQW,NMD)
        MPG = IM(IEQA,NMD)
        DPL = BLU(MPL)
        DPG = BLU(MPG)
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          MPS = IM(IEQS,NMD)
          DPS = BLU(MPS)
        ENDIF
!
!---    Saturated system w/o trapped gas
!       Water mass - aqueous pressure
!       CO2 mass - aqueous-CO2 mass fraction
!       NaCl mass - total NaCl brine mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
!
!---      Limit changes in pressure  ---
!
!          DPX = 2.5D-1*MAX(PL(2,N)+PATM,1.D+6)
          DPX = 1.D+6
          DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
          PL(2,N) = PL(2,N) + DPL
          PL(2,N) = MIN( PL(2,N),5.D+8 )
!
!---      Zero negative corrections for zero aqueous CO2  ---
!
          IF( XLA(2,N)/EPSL.LT.EPSL .AND. BLU(MPG)/EPSL.LT.EPSL ) THEN
            BLU(MPG) = 0.D+0
            DPG = 0.D+0
          ENDIF
          XLA(2,N) = MAX( (XLA(2,N)+DPG),0.D+0 )
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
!       Water mass - aqueous pressure
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 .OR. NPHAZ(2,N).EQ.4 
     &    .OR. NPHAZ(2,N).EQ.6 ) THEN
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
!
!---      Limit changes in pressure, and trap large changes in 
!         gas pressure  ---
!
          DPX = MAX( 1.D+6,2.5D-1*(PG(2,N)-PL(2,N)) )
          DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
!          IF( ABS(DPG).GT.1.D+6 ) IERR = 1
          DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
!
!---      Relax pressure updates when transitioning to unsaturated
!         conditions  ---
!
          IF( (PG(2,N)+DPG)-(PL(2,N)+DPL).LT.ENPR ) THEN
            DPG = 6.D-1*DPG
            DPL = 6.D-1*DPL
          ENDIF
          PG(2,N) = PG(2,N) + DPG
          PG(2,N) = MIN( PG(2,N),5.D+8 )
          PL(2,N) = PL(2,N) + DPL
          PL(2,N) = MAX( PL(2,N),(PG(2,N)-SCHR(12,IZN)*RHORL*GRAV) )
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
!---        Zero negative corrections for zero aqueous salt  ---
!
            IF( YLS(2,N)/EPSL.LT.EPSL .AND. BLU(MPS)/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS(2,N) = YLS(2,N)+DPS
            IF( YLS(2,N).LT.EPSL ) YLS(2,N) = 0.D+0
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
          ENDIF
!
!---      Maintain the gas pressure above or at the water vapor
!         pressure  ---
!
          CALL SP_B( T(2,N),XLS(2,N),PSBX )
          CALL DENS_B( T(2,N),PSBX,XLS(2,N),RHOBX )
          PCX = MAX( PSBX-PL(2,N),0.D+0 )
          CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(2,N),IZN )
          PG(2,N) = MAX( PG(2,N),(PVBX-PATM) )
!
!---    Saturated system w/ entrapped gas
!       Water mass - aqueous pressure
!       CO2 mass - trapped gas saturation
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.5
     &    .OR. NPHAZ(2,N).EQ.7 ) THEN
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
          DPX = 1.D-1*SCHR(15,IZ(N))
          DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
          SG(2,N) = SG(2,N) + DPG
!
!---      Relax disappearance of trapped gas  ---
!
!          IF( (SG(2,N)+DPG).LT.0.D+0 ) THEN
!            SG(2,N) = MAX( (SG(2,N)+6.D-1*DPG),0.D+0 )
!          ELSE
!            SG(2,N) = MAX( (SG(2,N)+DPG),0.D+0 )
!          ENDIF
!          DPX = MAX( 2.5D-1*SG(2,N),1.D-2 )
!          DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
!          SG(2,N) = SG(2,N) + DPG
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
!---    Fully unsaturated system
!       Water mass - aqueous saturation
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl mass concentration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.8 .OR. NPHAZ(2,N).EQ.9
     &    .OR. NPHAZ(2,N).EQ.10 ) THEN
!
!---      Limit changes in gas pressure  ---
!
          DPX = 2.5D-1*MAX(PG(2,N)+PATM,1.D+6)
          DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
          PG(2,N) = PG(2,N) + DPG
!
!---      Limit changes in water vapor pressure  ---
!
          IF( PVW(2,N)/EPSL.LT.EPSL .AND. BLU(MPL)/EPSL.LT.EPSL ) THEN
            BLU(MPL) = 0.D+0
            DPL = 0.D+0
          ENDIF
          CALL SOL_LS( T(2,N),XLSMX )
          XLSX = MIN( XLS(2,N),XLSMX )
          CALL SP_B( T(2,N),XLSX,PSBX )
          DPX = 2.5D-1*PSBX
          DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
          PVW(2,N) = PVW(2,N) + DPL
          IF( PVW(2,N).LT.1.D-6 ) PVW(2,N) = 0.D+0
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
            IF( TMS(2,N).LT.EPSL ) TMS(2,N) = 0.D+0
          ENDIF
!
!---      Maintain the gas pressure above or at the water vapor
!         pressure  ---
!
          CALL SP_B( T(2,N),XLSX,PSBX )
          CALL DENS_B( T(2,N),PSBX,XLSX,RHOBX )
          PCX = MAX( PSBX-PL(2,N),0.D+0 )
          CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,IZN )
          PG(2,N) = MAX( PG(2,N),(PVBX-PATM) )
       ENDIF
!
!---    Check for excessive pressure or temperature   ---
!
        PGX = PG(2,N)+PATM
        PLX = PL(2,N)+PATM
        TKX = T(2,N)+TABS
        IF( PGX.GT.8.D+8 .OR. PGX.LT.0.D+0 ) IERR = 1
        IF( PLX.GT.8.D+8 ) IERR = 1
        IF( TKX.GT.TCRW .OR. TKX.LT.TABS ) IERR = 1
        IF( IERR.EQ.1 ) NSDX = MIN( NSDX,N )
  200 CONTINUE
!
!---  Reduce time step for excessive changes in primary variables   ---
!
        IF( IERR.EQ.1 ) THEN
          ICNV = 1
          N = NSDX
          WRITE(ISC,'(10X,A)') '---  Excessive Primary Variable Change
     &---'
          WRITE(IWR,'(10X,A)') '---  Excessive Primary Variable Change
     &---'
          IF( N.GT.NFLD .AND. L_LW.EQ.1 ) THEN
            NWX = N - NFLD
            WRITE(ISC,'(4X,A,I6)') 'Leaky Well Node = ',NWX
            WRITE(IWR,'(4X,A,I6)') 'Leaky Well Node = ',NWX
          ELSE
            WRITE(ISC,'(4X,A,I6)') 'Node = ',N
            WRITE(IWR,'(4X,A,I6)') 'Node = ',N
          ENDIF
          WRITE(ISC,'(4X,2A)') 'Phase Condition = ',PH_CND(NPHAZ(2,N))
          WRITE(IWR,'(4X,2A)') 'Phase Condition = ',PH_CND(NPHAZ(2,N))
          WRITE(ISC,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',PL(2,N)+PATM
          WRITE(IWR,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',PL(2,N)+PATM
          WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Pressure = ',PG(2,N)+PATM
          WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Pressure = ',PG(2,N)+PATM
          WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG(2,N)
          WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG(2,N)
          WRITE(ISC,'(4X,A,1PE12.5)')
     &      'Aqueous-CO2 Mass Fraction = ',XLA(2,N)
          WRITE(IWR,'(4X,A,1PE12.5)')
     &      'Aqueous-CO2 Mass Fraction = ',XLA(2,N)
          WRITE(ISC,'(4X,A,1PE12.5)')
     &      'Water Vapor Pressure = ',PVW(2,N)
          WRITE(IWR,'(4X,A,1PE12.5)')
     &      'Water Vapor Pressure = ',PVW(2,N)
          WRITE(ISC,'(4X,A,1PE12.5,A,I6)')
     &      'Total-Salt Aqu. Mass Fraction = ',YLS(2,N)
          WRITE(IWR,'(4X,A,1PE12.5,A,I6)')
     &      'Total-Salt Aqu. Mass Fraction = ',YLS(2,N)
          WRITE(ISC,'(4X,A,1PE12.5,A,I6)')
     &      'Salt Volumetric Concentration = ',TMS(2,N)
          WRITE(IWR,'(4X,A,1PE12.5,A,I6)')
     &      'Salt Volumetric Concentration = ',TMS(2,N)
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
          DO 400 N = 1,NFLD+NWN_LW
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PVW(2,N) = PVW(1,N)
            XLA(2,N) = XLA(1,N)
            SG(2,N) = SG(1,N)
            SGT(2,N) = SGT(1,N)
            ASLMIN(2,N) = ASLMIN(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  400     CONTINUE
!
!---      Coupled-well pressure  ---
!
          IF( L_CW.EQ.1 ) THEN
            DO 402 NCW = 1,N_CW
              P_CW(2,NCW) = P_CW(1,NCW)
  402       CONTINUE
          ENDIF
!
!---    Number of time step reductions failure: stop simulation  ---
!
        ELSE
          DO 410 N = 1,NFLD+NWN_LW
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PVW(2,N) = PVW(1,N)
            XLA(2,N) = XLA(1,N)
            SG(2,N) = SG(1,N)
            SGT(2,N) = SGT(1,N)
            ASLMIN(2,N) = ASLMIN(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  410     CONTINUE
!
!---      Coupled-well pressure  ---
!
          IF( L_CW.EQ.1 ) THEN
            DO 412 NCW = 1,N_CW
              P_CW(2,NCW) = P_CW(1,NCW)
  412       CONTINUE
          ENDIF
          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          ICNV = 4
        ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_CO2 group
!
      RETURN
      END


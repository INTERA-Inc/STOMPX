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
!     STOMPX: Subsurface Transport Over Multiple Phases Extensible
!
!     STOMPX-W (Water) Mode
!
!     This engineering program numerically simulates the transport
!     of water through multifluid subsurface environments
!     under isothermal conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PROP
      USE OUTPU
      USE JACOB
      USE FLUX
      USE FILES
      USE FDVP
      USE GRID
      USE GLB_PAR
      USE GEO_MECH
      USE COUP_WELL
      USE BCV

!
!----------------------Lis Modules-----------------------------------!
!
      USE LIS_STOMP








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




      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = 1
      SUB_LOG(1) = 'STOMPX-W'
!
!---  Initialize MPI  ---
!
      CALL MPI_INIT( IERR )
!
!---  Get the individual process ID  ---
!
      CALL MPI_COMM_RANK( MPI_COMM_WORLD,ID,IERR )
!
!---  Get the number of processes  ---
!
      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD,NP,IERR )

!
!---  Initialize Lis ---
!
      CALL lis_init_f( IERR )
!      PRINT *,'LIS_INITIALIZE: IERR = ',IERR,'ID = ',ID







!
!---  Print banner to screen ---
!
      IF( ID.EQ.0 ) THEN
      PRINT *,' Welcome to ...'
      PRINT *,' '
      PRINT *,'                       STOMPX-W'
      PRINT *,'        Subsurface Transport Over Multiple Phases'
      PRINT *,'                 Water Operational Mode'
      PRINT *,'               OpenMPI Extensible Version'
      PRINT *,' '
      PRINT *,' This screen echo was produced by STOMPX, a numerical '
      PRINT *,' simulator developed by the Pacific Northwest National '
      PRINT *,' Laboratory. The simulator additionally produces'
      PRINT *,' binary formatted output files: output.bin, plot_xxx.bin'
      PRINT *,' and surface.bin, which can be converted to text format'
      PRINT *,' with postprocessing utilities.'
      PRINT *,' '
      PRINT *,'                     Disclaimer'
      PRINT *,' This material was prepared as an account of work '
      PRINT *,' sponsored by an agency of the United States Government.'
      PRINT *,' Neither the United States Government nor the United '
      PRINT *,' States Department of Energy, nor Battelle, nor any of '
      PRINT *,' their employees, makes any warranty, express or '
      PRINT *,' implied, or assumes any legal liability or '
      PRINT *,' responsibility for the accuracy, completeness, or '
      PRINT *,' usefulness of any information, apparatus, product, '
      PRINT *,' software or process disclosed, or represents that its '
      PRINT *,' use would not infringe privately owned rights.'
      PRINT *,' '
      PRINT *,'                     Acknowledgement'
      PRINT *,' This software and its documentation were produced with '
      PRINT *,' Government support under Contract Number '
      PRINT *,' DE-AC06-76RLO-1830 awarded by the United Department of '
      PRINT *,' of Energy. The Government retains a paid-up '
      PRINT *,' non-exclusive, irrevocable worldwide license to'
      PRINT *,' implied, or assumes any legal liability or '
      PRINT *,' reproduce, prepare derivative works, perform publicly '
      PRINT *,' and display publicly by or for the Government, '
      PRINT *,' including the right to distribute to other Government '
      PRINT *,' contractors.'
      PRINT *,' '
      PRINT *,' For support:  Tel: 509.372.6070'
      PRINT *,'               E-mail:  mark.white@pnnl.gov'
      PRINT *,' '
      PRINT *,'                     Notice'
      PRINT *,' This screen echo only reports basic time stepping '
      PRINT *,' information to track the simulation progress. Expanded'
      PRINT *,' simulation results are reported in binary formatted '
      PRINT *,' output files: output.bin, plot_xxx.bin and '
      PRINT *,' surface.bin.'
      PRINT *,' '
      ENDIF
!
!---  Read binary files from preprocessor and allocate
!     memory for global arrays  ---
!
      CALL READ_BIN_W
!      PRINT *,'Post INTLZ_W: ID = ',ID
!      ISKIP = 1
!      IF( ISKIP.EQ.0 ) THEN
!
!---  Initialize global arrays  ---
!
      CALL INTLZ_W
!      PRINT *,'Post INTLZ_W: ID = ',ID
!
!---  Configure the arrays for compressed sparse row matrix storage  ---
!
      CALL JCBP_W
!      PRINT *,'Post JCBP_W: ID = ',ID
!
!---  Compute primary variable increments  ---
!
      CALL INCRM_W
!      PRINT *,'Post INCRM_W: ID = ',ID
!
!---  Hydrologic and thermodynamic properties  ---
!
      CALL PROP_W
!      PRINT *,'Post PROP_W: ID = ',ID
      CALL BCP_W
!      PRINT *,'Post BCP_W: ID = ',ID
!
!---  Compute initial solute concentrations  ---
!
      CALL CISC_W
!      PRINT *,'Post CISC_W: ID = ',ID
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
        DO NEQ = 1,NEQC
          NSL = NEQ + NSOLU
!
!---      Mobile conservation component fractions   ---
!
          CALL MOBCF( NEQ )
!
!---      Add immobile conservation component fractions   ---
!
          CALL IMOBCF( NEQ )
        ENDDO
!
!---    Loop over number of kinetic component species  ---
!
        DO NEQ = 1,NEQK
          NSL = NEQ + NEQC + NSOLU
! 
!---      Mobile kinetic component fractions   ---
!
          CALL MOBKF( NEQ )
! 
!---      Add immobile kinetic component fractions   ---
!
          CALL IMOBKF( NEQ )
        ENDDO
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
!---  Communicate the state of field nodes with coupled well nodes
!     to all processors, increment coupled-well primary variables  ---
!
      IF( L_CW.EQ.1 ) THEN
        CALL COMST_COUP_WELL
!        PRINT *,'Post COMST_COUP_WELL: ID = ',ID
        CALL INCRM_COUP_WELL
!        PRINT *,'Post INCRM_COUP_WELL: ID = ',ID
      ENDIF
!
!---  Compute initial fluxes on non-boundary and boundary surfaces  ---
!
      ISVF = 1
      CALL FLUX_W
      CALL BCF_W
      ISVF = 2*ISVC + 1
!
!---  Surface flux integrator for zero time step  ---
!
      DTX = DT
      DT = 0.D+0
      CALL SFIN_W
      DT = DTX

!
!---  Create Lis matrix, solver, and solution and problem vectors
!     for coupled flow  ---
!
!      PRINT *,'NUKFL(ID+1) = ',NUKFL(ID+1),' ID = ',ID
!      PRINT *,'NUKFO(ID+1) = ',NUKFO(ID+1),' ID = ',ID
!      PRINT *,'NUKFG = ',NUKFG,' ID = ',ID
      INDX = 0
      CALL STOMP_LIS_CREATE( F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC,
     &  NUKFL(ID+1),INDX )
!      PRINT *,'Post Flow STOMP_LIS_CREATE: ID = ',ID
!
!---    Create Lis matrix, solver, and solution and problem vectors
!       for solute transport  ---
!
      IF( IEQC.NE.0 .OR. ISLC(40).NE.0 ) THEN
!        PRINT *,'NUKTL(ID+1) = ',NUKTL(ID+1),' ID = ',ID
!        PRINT *,'NUKTO(ID+1) = ',NUKTO(ID+1),' ID = ',ID
!        PRINT *,'NUKTG = ',NUKTG,' ID = ',ID
        INDX = 1
        CALL STOMP_LIS_CREATE( T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC,
     &    NUKTL(ID+1),INDX )
!        PRINT *,'Post Transport STOMP_LIS_CREATE: ID = ',ID
      ENDIF
!
!---    Create Lis matrix, solver, and solution and problem vectors
!       for geomechanics  ---
!
      IF( ISLC(50).NE.0 ) THEN
!        PRINT *,'NUKGL(ID+1) = ',NUKGL(ID+1),' ID = ',ID
!        PRINT *,'NUKGO(ID+1) = ',NUKGO(ID+1),' ID = ',ID
!        PRINT *,'NUKGG = ',NUKGG,' ID = ',ID
        INDX = 2
        CALL STOMP_LIS_CREATE( G_KSP,G_MAT,G_RHS_VEC,G_SOL_VEC,
     &    NUKGL(ID+1),INDX )
!        PRINT *,'Post Geomechanics STOMP_LIS_CREATE: ID = ',ID
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
!---  Check for fatal execution errors and stop simulation
!     if detected  ---
!
      CALL CHK_ERROR
!
!---  New Time Step ---
!
      IF( ID.EQ.0 ) ICNO = 10
      TSLOOP: DO
!
!---    Reference node(s) output  ---
!
        IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 ) THEN
          CALL REFNOD_W
          IF( ID.EQ.0 ) THEN
            IF( ICNO.EQ.10 ) THEN
              ICNO = 0
              NCH = INDEX(UNTM(1:),'  ') - 1
              PRINT *,'       Step         Itr           Time' // 
     &          ' [',UNTM(1:NCH),']             Timestep [',
     &          UNTM(1:NCH),']'
            ENDIF
            ICNO = ICNO + 1
            PRINT *,NSTEP,NITER,TM*CNVTM,DT*CNVTM
          ENDIF
        ENDIF
!
!---    Normalize mineral species concentrations after initial
!       output for normal simulations  ---
!
        IF( ISLC(40).EQ.1 .AND. (NSTEP-NRST).EQ.0 ) CALL NMNSP
!
!---    Load old time step arrays  ---
!
        CALL LDO_W
!
!---    Load old time step arrays for the coupled-well model  ---
!
        IF( L_CW.EQ.1 ) THEN
          CALL LDO_COUP_WELL
!          PRINT *,'Post LDO_COUP_WELL: ID = ',ID
        ENDIF
!
!---    Compute trapping number  ---
!
        CALL TRPGL_W
!
!---    Stop simulation if simulation time exceeds limit  ---
!
        IF( ABS(TMMX-TM).LE.1.D-6 ) THEN
          IF( ID.EQ.0 ) PRINT *,'Simulation Stopped:  ' // 
     &      'Simulation Time Limit'
          EXIT TSLOOP
        ENDIF
!
!---    Restart and plot file outputs  ---
!
        IF( ABS(TMPR-TM).LE.1.D-6 ) THEN
          CALL WRPLOT_W
          IF( ISLC(18).LT.1 ) CALL WRRST_W
        ENDIF
!
!---    Compute the next time step and increment time step counter  ---
!
        DTSO = DT
        CALL TMSTEP
        IF( NSTEP.EQ.0 ) DTSO = DT
        NSTEP = NSTEP + 1
        IF( NSTEP-NRST.GT.MXSTEP ) THEN
          IF( ID.EQ.0 ) PRINT *,'Simulation Stopped:  Time Step Limit'
          NSTEP = NSTEP - 1
          EXIT TSLOOP
        ENDIF
!
!---    New Newton-Raphson iteration ---
!
        NITER = 0
        NRLOOP: DO
          IF( ICNV.EQ.1 ) NITER = 0
          NITER = NITER + 1
!
!---      Compute boundary saturation, relative permeability, and
!         thermodynamic properties  ---
!
          CALL BCP_W
!
!---      Compute coupled-well fluxes  ---
!
          IF( L_CW.EQ.1 ) CALL FLUX_COUP_WELL
!          PRINT *,'P_CW(2,1) = ',P_CW(2,1),'P_CW(3,1) = ',P_CW(3,1),
!     &      'ID = ',ID
!          PRINT *,'QM_CW(1,1) = ',QM_CW(1,1),'ID = ',ID
!          PRINT *,'QM_CW(2,1) = ',QM_CW(2,1),'ID = ',ID
!          PRINT *,'QM_CW(3,1) = ',QM_CW(3,1),'ID = ',ID
!          PRINT *,'Post FLUX_COUP_WELL: ID = ',ID
!
!---      Compute source contributions  ---
!
          CALL SORC_W
!
!---      Compute fluxes on non-boundary and boundary surfaces  ---
!
          CALL FLUX_W
          CALL BCF_W
!
!---      Zero Jacobian matrix  ---
!
          CALL JCBZ_W
!
!---      Load Jacobian matrix for the water, air
!         and salt mass equations (zero flux boundary)  ---
!
          CALL JCB_W
!
!---      Modify the Jacobian matrix for boundary conditions  ---
!
          CALL BCJ_W
!
!---      Modify Jacobian matrix for coupled-well equations  ---
!
          IF( L_CW.EQ.1 ) THEN
            CALL JCB_COUP_WELL
!             PRINT *,'Post JCB_COUP_WELL:ID = ',ID
          ENDIF
!
!---      Set values of the Jacobian matrix  ---
!
          CALL JCB_SV
!
!---      Check for fatal execution errors and stop simulation
!         if detected  ---
!
          CALL CHK_ERROR
!
!---      Solve the linear system A x = b for coupled flow  ---
!

          INDX = 0
          CALL STOMP_LIS_SOLVE( F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC,
     &      NUKFO(ID+1),NUKFL(ID+1),INDX )
!          PRINT *,'NUKFL(ID+1) = ',NUKFL(ID+1),'NUKFG = ',NUKFG,
!     &      'ID = ',ID
          IF( ICNV.EQ.4 ) EXIT TSLOOP 







!
!---      Check for fatal execution errors and stop simulation
!         if detected  ---
!
          CALL CHK_ERROR
!
!---      Update primary variables for coupled wells  ---
!
          IF( L_CW.EQ.1 ) THEN
            CALL UPDT_COUP_WELL
!            PRINT *,'P_CW(2,1) = ',P_CW(2,1),'ID = ',ID
!            PRINT *,'Post UPDT_COUP_WELL:ID = ',ID
          ENDIF
!
!---      Update primary variables on field nodes w/o ghost cells  ---
!
          CALL UPDT_W

!
!---      Update primary variables on ghost cells  ---
!
          CALL UPDT_GC_W

!
!---      Convergence check for couped wells  ---
!
          IF( L_CW.EQ.1 ) THEN
            CALL RSDL_COUP_WELL
!            PRINT *,'Post RSDL_COUP_WELL:ID = ',ID
          ENDIF
!
!---      Compute convergence from maximum relative residuals  ---
!
          CALL RSDL_W
!
!---      Compute primary variable increments, saturation,
!         relative permeability, porosity, tortuosity,
!         thermodynamic properties for interior nodes,
!         except immediately after a new time step  ---
!
          CALL INCRM_W
          CALL PROP_W
!
!---      Increment coupled-well primary variables  ---
!
          IF( L_CW.EQ.1 ) THEN
            CALL INCRM_COUP_WELL
!            PRINT *,'Post INCRM_COUP_WELL:ID = ',ID
          ENDIF
!
!---      Convergence check  ---
!
!         ICNV = 1 - cut time step and restart Newton-Raphson loop
!         ICNV = 2 - next Newton-Raphson loop
!         ICNV = 3 - converged solution next time step
!         ICNV = 4 - total convergence failure stop simulation
!
!          IF( ID.EQ.0 ) THEN
!            PRINT *,'NSTEP = ',NSTEP,' NITER = ',NITER
!            DO M = 1,ISVC
!              PRINT *,'RSD(',M,') = ',RSD(M),
!     &          'NSD(',M,') = ',NSD(M),' ICNV = ',ICNV
!            ENDDO
!            PRINT *,'RSD_CW = ',RSD_CW,' ID_CW(8,1) = ',ID_CW(8,1),
!     &        'P_CW(2,1) = ',P_CW(2,1),'P_CW(3,1) = ',P_CW(3,1)
!          ENDIF
          IF( ICNV.EQ.3 ) EXIT NRLOOP
          IF( ICNV.EQ.4 ) EXIT TSLOOP
!
!---    Proceed to new Newton-Raphson iteration  ---
!
        ENDDO NRLOOP
!
!---    Compute current fluxes for transport solutions or flux
!       integrations  ---
!
        ISVF = 1
        CALL FLUX_W
        CALL BCF_W
        ISVF = 2*ISVC + 1
!
!---    Compute Local Courant Numbers  ---
!
        IF( ICRNT.EQ.1 ) CALL CRNTNB
        ISVF = 2*ISVC+1
!
!---    Beginning of transport equation solution  ---
!
        IF( IEQC.NE.0 .OR. ISLC(40).NE.0 ) THEN
!
!---      Loop over number of solutes  ---
!
          DO NSL = 1,NSOLU
!
!---        Courant number limiting  ---
!
            N_CRN(NSL) = 1
            IF( ISLC(17).NE.0 ) CALL CRN_LIM( NSL )
!            PRINT *,'N_CRN(',NSL,') = ',N_CRN(NSL),' ID = ',ID
!
!---        Sub-time step loop  ---
!
            DO NC = 1,N_CRN(NSL)
              IF( ISLC(17).NE.0 ) TM = MIN( TM+DT,TM_CRN )
!
!---          Compute solute mole fractions ---
!
              CALL SPRP_W( NSL )
!
!---          Solute transport ---
!
              CALL TPORT_W( NSL )
              IF( ICNV.EQ.4 ) EXIT TSLOOP 
!
!---        Bottom of sub-time step loop  ---
!
            ENDDO
!
!---        Courant number limiting, reset time stepping  ---
!
            IF( ISLC(17).NE.0 ) THEN
              DT = DT_CRN
              DTI = DTI_CRN
              TM = TM_CRN
            ENDIF
          ENDDO
!
!---      Decay matrix, fracture, and borehole solutes via Bateman
!         chain decay solution  ---
!
          CALL CHAIN_DECAY
!
!---      Reactive transport  ---
!
          IF( ISLC(40).EQ.1 ) THEN
            N_CRN(NSOLU+1) = 1
            IF( ISLC(17).NE.0 ) CALL CRN_LIM( NSOLU+1 )
!
!---        Courant-limiting sub-time step loop  ---
!
            DO NCR = 1,N_CRN(NSOLU+1)
              IF( ISLC(17).NE.0 ) TM = MIN( TM+DT,TM_CRN )
!
!---          Temporarily store time stepping  ---
!
              DT_RST = DT
              DTI_RST = DTI
              TM_RST = TM
              TM = TM - DT
              N_RST = 1
!
!---          Top of ECKEChem time-step reduction loop  ---
!
              TRLOOP: DO
!
!---            Zero linked sources  ---
!
                CALL ZLKSRC
!
!---            Sub-time step reduction limit exceeded  ---
!
                IF( N_RST.GT.16 ) THEN
                  IF( ID.EQ.0 ) PRINT *, '          ---  ECKEChem ' // 
     &              'Sub-Time Step Reduction Limit Exceeded  ---'
                  DT = DT_RST
                  DTI = DTI_RST
                  TM = TM_RST
                  NSTEP = NSTEP-1
                  TM = TM-DT
                  DT = DTO
                  CALL BCK_STP
                  EXIT TSLOOP
                ENDIF
!
!---            Sub-time step loop  ---
!
                DO NC = 1,N_RST
                  TM = TM + DT
!
!---              Loop over number of conservation component species ---
!
                  DO NEQ = 1,NEQC
                    NSL = NEQ + NSOLU
!
!---                Skip transport for linked aqueous air   ---
!
                    IF( ISPLK(6).EQ.NSL ) THEN
                      CALL IMOBCF( NEQ )
!
!---                Transport for conservation component species   ---
!
                    ELSE
!
!---                  Mobile conservation component fractions   ---
!
                      CALL MOBCF( NEQ )
!
!---                  Solute transport ---
!
                      CALL TPORT_W( NSL )
                      IF( ICNV.EQ.4 ) EXIT TSLOOP 
!
!---                  Add immobile conservation component fractions  ---
!
                      CALL IMOBCF( NEQ )
                    ENDIF
!
!---              End of conservation component species transport  ---
!
                  ENDDO
!
!---              Loop over number of kinetic component species  ---
!
                  DO NEQ = 1,NEQK
                    NSL = NEQ + NEQC + NSOLU
!
!---                Skip transport for linked aqueous air   ---
!
                    IF( ISPLK(6).EQ.NSL ) THEN
                      CALL IMOBKF( NEQ )
!
!---                Transport for mobile kinetic component species   ---
!
                    ELSE
! 
!---                  Mobile kinetic component fractions   ---
!
                      CALL MOBKF( NEQ )
!
!---                  Solute transport ---
!
                      CALL TPORT_W( NSL )
                      IF( ICNV.EQ.4 ) EXIT TSLOOP 
! 
!---                  Add immobile kinetic component fractions   ---
!
                      CALL IMOBKF( NEQ )
                    ENDIF
!
!---              End of kinetic component species transport  ---
!
                  ENDDO
!
!---              Equilibrium-conservation-kinetic reaction 
!                 chemistry  ---
!
                  CALL ECKECHEM
                  IF( ECKE_ER ) CYCLE TRLOOP
!
!---              Load old sub-time-step reactive species concentrations 
!                 and component species concentrations  ---
!        
                  IF( ISLC(17).NE.0 ) CALL UPDTCHEM
!
!---            Bottom of sub-time step loop  ---
!
                ENDDO
!
!---             Exit time-step reduction loop  ---
!
                 EXIT TRLOOP
!
!---          Bottom of ECKEChem time-step reduction loop  ---
!
              ENDDO TRLOOP
!
!---          Reset time stepping  ---
!
              IF( N_RST.GT.1 ) THEN
                DT = DT_RST
                DTI = DTI_RST
                TM = TM_RST
              ENDIF
            ENDDO
!
!---        Courant number limiting, reset time stepping  ---
!
            IF( ISLC(17).NE.0 ) THEN
              DT = DT_CRN
              DTI = DTI_CRN
              TM = TM_CRN
            ENDIF
          ENDIF
        ENDIF
!
!---    Surface flux integrator  ---
!
        CALL SFIN_W
!
!---  Proceed to new time step  ---
!
      ENDDO TSLOOP

!
!---  Destroy Lis matrix, solver, and solution and problem vectors
!     for coupled flow  ---
!
      CALL STOMP_LIS_DESTROY( F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC )
!
!---  Destroy Lis matrix, solver, and solution and problem vectors
!     for solute transport  ---
!
      IF( IEQC.NE.0 .OR. ISLC(40).NE.0 ) THEN
        CALL STOMP_LIS_DESTROY( T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC )
      ENDIF
!
!---  Destroy Lis matrix, solver, and solution and problem vectors
!     for geomechanics  ---
!
      IF( ISLC(50).NE.0 ) THEN
        CALL STOMP_LIS_DESTROY( G_KSP,G_MAT,G_RHS_VEC,G_SOL_VEC )
      ENDIF
!
!---  Finalize Lis execution  ---
!
      CALL lis_finalize_f( IERR )

!
!---  Write plot_xxx.bin file  ---
!
      CALL WRPLOT_W
!
!---  Write restart_xxx.bin file  ---
!
      CALL WRRST_W
!
!---  Close output.bin file, putting a closing -9 at the end
!     of the file  ---
!
      NVAR = 1
      IVARX = -9
      OFFSET = IOFFSET_REF
      IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &   MPI_INTEGER,STATUS,IERR )
      CALL MPI_FILE_CLOSE( IWR,IERR )
!
!---  Close surface.bin file(s), putting a closing -9.D+20 at the end
!     of the file  ---
!
      NSTART = 1
      IF( ISFGP(1).EQ.0 ) NSTART = 2
      DO NSG = NSTART,NSFGP
        NVAR = 1
        VARX = -9.D+20
        OFFSET = IOFFSET_SF(NSG)
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,VARX,
     &     NVAR,MPI_REAL8,STATUS,IERR )
        CALL MPI_FILE_CLOSE( ISF(NSG),IERR )
      ENDDO
!      ENDIF
!      PRINT *,'Fini STOMPX-W: ID = ',ID
      CALL MPI_FINALIZE(IERR)
!
!---  End of STOMP program  ---
!
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCF_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Compute boundary surface fluxes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/BCF_W'
      DO NB = 1,NBC(ID+1)
        N = IBCN(NB)
        IF( IBCD(NB).EQ.-3 ) THEN
          DO M = 1,ISVF
            WL(M,1,N) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          DO M = 1,ISVF
            VL(M,1,N) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          DO M = 1,ISVF
            UL(M,1,N) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          DO M = 1,ISVF
            UL(M,2,N) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          DO M = 1,ISVF
            VL(M,2,N) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          DO M = 1,ISVF
            WL(M,2,N) = 0.D+0
          ENDDO
        ENDIF
      ENDDO
!
!---  Loop over boundary conditions  ---
!
      DO NB = 1,NBC(ID+1)
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) CYCLE
!
!---  Assign local boundary condition variables  ---
!
        IF( IBCM(NB).EQ.1 ) THEN
          DO N = 1,LBCV
            BCX(N) = BC(N,1,MB)
          ENDDO
        ELSE
          IFIND = 0
          DO M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
             TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
             DTBC = MIN( BC(1,M,MB)-TMZ,DT )
             TFBC = (TMZ-BC(1,M-1,MB))/TDBC
             DO N = 1,LBCV
               BCX(N) = BC(N,M-1,MB) + TFBC*(BC(N,M,MB)-BC(N,M-1,MB))
             ENDDO
             IFIND = 1
             EXIT
            ENDIF
          ENDDO
          IF( IFIND.EQ.0 ) CYCLE
        ENDIF
        N = IBCN(NB)
        N_DB = -NB
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              WL(M,1,N) = BCX(2)
            ENDDO
            CALL DFFLWB_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).EQ.24 ) THEN
            DO M = 1,ISVF
              WL(M,1,N) = BCX(2)
            ENDDO
            CALL DFFLWB_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLB_W( N,NB )
            CALL DFFLWB_W( N,NB )
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              VL(M,1,N) = BCX(2)
            ENDDO
            CALL DFFLWS_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).EQ.24 ) THEN
            DO M = 1,ISVF
              MP = MPOS(M)
              HDGL = MAX( (PG(MP,N)-PL(MP,N))/RHORL/GRAV,1.D-14 )
              PEFX = 1.D+0
              IF( HDGL.GT.BCX(3) ) THEN
                PEFX = (BCX(3)/HDGL)**4
              ENDIF
              VL(M,1,N) = -ABS(BCX(2))*PEFX
            ENDDO
            CALL DFFLWS_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLS_W( N,NB )
            CALL DFFLWS_W( N,NB )
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              UL(M,1,N) = BCX(2)
            ENDDO
            CALL DFFLWW_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).EQ.24 ) THEN
            DO M = 1,ISVF
              MP = MPOS(M)
              HDGL = MAX( (PG(MP,N)-PL(MP,N))/RHORL/GRAV,1.D-14 )
              PEFX = 1.D+0
              IF( HDGL.GT.BCX(3) ) THEN
                PEFX = (BCX(3)/HDGL)**4
              ENDIF
              UL(M,1,N) = -ABS(BCX(2))*PEFX
            ENDDO
            CALL DFFLWW_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLW_W( N,NB )
            CALL DFFLWW_W( N,NB )
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              UL(M,2,N) = BCX(2)
            ENDDO
            CALL DFFLWE_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).EQ.24 ) THEN
            DO M = 1,ISVF
              MN = MNEG(M)
              HDGL = MAX( (PG(MN,N)-PL(MN,N))/RHORL/GRAV,1.D-14 )
              PEFX = 1.D+0
              IF( HDGL.GT.BCX(3) ) THEN
                PEFX = (BCX(3)/HDGL)**4
              ENDIF
              UL(M,2,N) = ABS(BCX(2))*PEFX
            ENDDO
            CALL DFFLWE_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLE_W( N,NB )
            CALL DFFLWE_W( N,NB )
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              VL(M,2,N) = BCX(2)
            ENDDO
            CALL DFFLWN_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).EQ.24 ) THEN
            DO M = 1,ISVF
              MN = MNEG(M)
              HDGL = MAX( (PG(MN,N)-PL(MN,N))/RHORL/GRAV,1.D-14 )
              PEFX = 1.D+0
              IF( HDGL.GT.BCX(3) ) THEN
                PEFX = (BCX(3)/HDGL)**4
              ENDIF
              VL(M,2,N) = ABS(BCX(2))*PEFX
            ENDDO
            CALL DFFLWN_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLN_W( N,NB )
            CALL DFFLWN_W( N,NB )
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              WL(M,2,N) = BCX(2)
            ENDDO
            CALL DFFLWT_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).EQ.24 ) THEN
            DO M = 1,ISVF
              WL(M,2,N) = BCX(2)
            ENDDO
            CALL DFFLWT_W( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLT_W( N,NB )
            CALL DFFLWT_W( N,NB )
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCF_W group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCJ_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Modify the Jacobian matrix for boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
      USE FDVP
      USE CONST
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
      SUB_LOG(ISUB_LOG) = '/BCJ_W'
!
!---  Loop over boundary conditions  ---
!
      DO NB = 1,NBC(ID+1)
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) CYCLE
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) CYCLE
        N = IBCN(NB)
        N_DB = -NB
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWB_W( N,NB )
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWS_W( N,NB )
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWW_W( N,NB )
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWE_W( N,NB )
          ENDIF
!
!---    North boundary
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWN_W( N,NB )
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWT_W( N,NB )
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCJ_W group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCP_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Compute saturation, relative permeability and thermodynamic
!     properties for boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE HYST
      USE GRID
      USE FLUX
      USE FDVP
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
      REAL*8 BCX(LBCV),YLSX(3)
      REAL*8 GX(3),RX(2),RPX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCP_W'
!
!---  Assign values for initial condition type boundary conditions  ---
!
      IF( NSTEP-NRST.LE.1 .AND. NITER.LE.1 ) THEN
        DO NB = 1,NBC(ID+1)
          N = IBCN(NB)
          IF( IBCT(IEQW,NB).EQ.12 .OR. IBCT(IEQA,NB).EQ.12 ) THEN
            IF( IBCD(NB).EQ.-3 ) THEN
              DB = DZGP(1,N)
              GB = GRVZ(1,N)*DB
            ELSEIF( IBCD(NB).EQ.-2 ) THEN
              DB = DYGP(1,N)*RP(N)
              GB = GRVY(1,N)*DB
            ELSEIF( IBCD(NB).EQ.-1 ) THEN
              DB = DXGP(1,N)
              GB = GRVX(1,N)*DB
            ELSEIF( IBCD(NB).EQ.1 ) THEN
              DB = -DXGP(2,N)
              GB = GRVX(2,N)*DB
            ELSEIF( IBCD(NB).EQ.2 ) THEN
              DB = -DYGP(2,N)*RP(N)
              GB = GRVY(2,N)*DB
            ELSEIF( IBCD(NB).EQ.3 ) THEN
              DB = -DZGP(2,N)
              GB = GRVZ(2,N)*DB
            ENDIF
          ENDIF
          IF( IBCT(IEQW,NB).EQ.12 ) THEN
            PLB(1,NB) = PL(2,N) + RHOL(2,N)*GB
            PLB(2,NB) = PLB(1,NB)
          ENDIF
        ENDDO
      ENDIF
!
!---  Loop over boundary conditions  ---
!
      L1: DO NB = 1,NBC(ID+1)
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) CYCLE
!
!---    Assign local boundary condition variables  ---
!
        IF( IBCM(NB).EQ.1 ) THEN
          DO N = 1,LBCV
            BCX(N) = BC(N,1,MB)
          ENDDO
        ELSE
          IFIND = 0
          DO M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
             TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
             DTBC = MIN( BC(1,M,MB)-TMZ,DT )
             TFBC = (TMZ-BC(1,M-1,MB))/TDBC
             DO N = 1,LBCV
               BCX(N) = BC(N,M-1,MB) + TFBC*(BC(N,M,MB)-BC(N,M-1,MB))
             ENDDO
             IF( IBCT(IEQW,NB).EQ.2 ) THEN
               BCX(2) = BCX(2)-5.D-1*DTBC*(BC(2,M,MB)-BC(2,M-1,MB))/TDBC
             ENDIF
             IFIND = 1
             EXIT
            ENDIF
          ENDDO
          IF( IFIND.EQ.0 ) CYCLE
        ENDIF
!
!---    Initial condition boundary condition  ---
!
        IF( IBCT(IEQW,NB).EQ.12 ) BCX(2) = PLB(1,NB)
        N = IBCN(NB)
        IBD = ABS(IBCD(NB))
        N_DB = -NB
!
!---    ECKEChem  ---
!
        POR0(1,N) = POR0(1,N)
        POR0(2,N) = POR0(2,N)
!
!---    Boundary Direction  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
          DB = 0.5D+0*DZGF(N)
          GB = GRVZ(1,N)*DB
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          DB = 0.5D+0*DYGF(N)*RP(N)
          GB = GRVY(1,N)*DB
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          DB = 0.5D+0*DXGF(N)
          GB = GRVX(1,N)*DB
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          DB = -0.5D+0*DXGF(N)
          GB = GRVX(2,N)*DB
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          DB = -0.5D+0*DYGF(N)*RP(N)
          GB = GRVY(2,N)*DB
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          DB = -0.5D+0*DZGF(N)
          GB = GRVZ(2,N)*DB
        ENDIF
!
!---    Direction gradients for x-y-z hydraulic gradient
!       or seepage face boundaries.
!
        IF( ABS(IBCT(IEQW,NB)).EQ.44 .OR. ABS(IBCT(IEQW,NB)).EQ.45) THEN
          DXX = XPBC(NB) - XSBC(MB)
          DYX = YPBC(NB) - YSBC(MB)
          DZX = ZPBC(NB) - ZSBC(MB)
        ENDIF
!        IF( ABS(IBCT(IEQW,NB)).NE.11 ) PRINT *,'MB = ',MB,' NB = ',NB,
!     &    ' N = ',ND(N),' IBD = ',IBD,' ID = ',ID
!
!---    Loop over secondary variable indices  ---
!
        DO M = 2,ISVC+2
          PLX = PL(M,N)
          PGX = PG(2,N)
          TX = T(2,N)
!
!---      Aqueous Dirichlet  ---
!
          IF( IBCT(IEQW,NB).EQ.1 .OR. IBCT(IEQW,NB).EQ.12 ) THEN
            PLX = BCX(2)
!
!---      Aqueous Neumann  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.2 ) THEN
            PLX = PLX + BCX(2)*DB*VISL(M,N)/PERM(IBD,N)
     &        + RHOL(M,N)*GB
!
!---      Aqueous Zero Flux  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.3 ) THEN
            PLX = PLX + RHOL(M,N)*GB
!
!---      Aqueous Saturated  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.4 ) THEN
            PLX = PGX
!
!---      Aqueous Unit Gradient  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.5 ) THEN
            PLX = PLX
!
!---      Aqueous Free Gradient  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.6 ) THEN
            IF( IBCD(NB).EQ.-3 ) THEN
              IF( IXP(ICM(6,N)).NE.0 ) THEN
                PCT = PL(1,ICM(6,N))
                PCP = PL(1,N)
                PLX = PCP + 0.5D+0*(PCP-PCT)*DZGF(N)/DZGP(2,N)
              ELSE
                PLX = PLX + RHOL(M,N)*GB
              ENDIF
            ELSEIF( IBCD(NB).EQ.-2 ) THEN
              IF( IXP(ICM(5,N)).NE.0 ) THEN
                PCN = PL(1,ICM(5,N))
                PCP = PL(1,N)
                PLX = PCP + 0.5D+0*(PCP-PCN)*DYGF(N)/DYGP(2,N)
              ELSE
                PLX = PLX + RHOL(M,N)*GB
              ENDIF
            ELSEIF( IBCD(NB).EQ.-1 ) THEN
              IF( IXP(ICM(4,N)).NE.0 ) THEN
                PCE = PL(1,ICM(4,N))
                PCP = PL(1,N)
                PLX = PCP + 0.5D+0*(PCP-PCE)*DXGF(N)/DXGP(2,N)
              ELSE
                PLX = PLX + RHOL(M,N)*GB
              ENDIF
            ELSEIF( IBCD(NB).EQ.1 ) THEN
              IF( IXP(ICM(3,N)).NE.0 ) THEN
                PCW = PL(1,ICM(3,N))
                PCP = PL(1,N)
                PLX = PCP + 0.5D+0*(PCP-PCW)*DXGF(N)/DXGP(2,N)
              ELSE
                PLX = PLX + RHOL(M,N)*GB
              ENDIF
            ELSEIF( IBCD(NB).EQ.2 ) THEN
              IF( IXP(ICM(2,N)).NE.0 ) THEN
                PCS = PL(1,ICM(2,N))
                PCP = PL(1,N)
                PLX = PCP + 0.5D+0*(PCP-PCS)*DYGF(N)/DYGP(2,N)
              ELSE
                PLX = PLX + RHOL(M,N)*GB
              ENDIF
            ELSEIF( IBCD(NB).EQ.3 ) THEN
              IF( IXP(ICM(1,N)).NE.0 ) THEN
                PCB = PL(1,ICM(1,N))
                PCP = PL(1,N)
                PLX = PCP + 0.5D+0*(PCP-PCB)*DZGF(N)/DZGP(2,N)
              ELSE
                PLX = PLX + RHOL(M,N)*GB
              ENDIF
            ENDIF
!
!---  Aqueous Outflow  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.7 ) THEN
            PLX = BCX(2)
!
!---  Aqueous Free Boundary  ---
!
          ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
            PLX = PL(2,N)
            IF( IBCD(NB).EQ.-3 .AND. IXP(ICM(6,N)).NE.0 ) THEN
              PLX = PL(2,N) + 0.5D+0*(PL(2,N)-PL(2,ICM(6,N)))*DZGP(1,N)/
     &          DZGP(2,N)
            ELSEIF( IBCD(NB).EQ.-2 .AND. IXP(ICM(5,N)).NE.0 ) THEN
              PLX = PL(2,N) + 0.5D+0*(PL(2,N)-PL(2,ICM(5,N)))*DYGP(1,N)/
     &          DYGP(2,N)
            ELSEIF( IBCD(NB).EQ.-1 .AND. IXP(ICM(4,N)).NE.0 ) THEN
              PLX = PL(2,N) + 0.5D+0*(PL(2,N)-PL(2,ICM(6,N)))*DXGP(1,N)/
     &          DXGP(2,N)
            ELSEIF( IBCD(NB).EQ.1 .AND. IXP(ICM(3,N)).NE.0 ) THEN
              PLX = PL(2,N) + 0.5D+0*(PL(2,N)-PL(2,ICM(6,N)))*DXGP(2,N)/
     &          DXGP(1,N)
            ELSEIF( IBCD(NB).EQ.2 .AND. IXP(ICM(2,N)).NE.0 ) THEN
              PLX = PL(2,N) + 0.5D+0*(PL(2,N)-PL(2,ICM(6,N)))*DYGP(2,N)/
     &          DYGP(1,N)
            ELSEIF( IBCD(NB).EQ.3 .AND. IXP(ICM(1,N)).NE.0 ) THEN
              PLX = PL(2,N) + 0.5D+0*(PL(2,N)-PL(2,ICM(6,N)))*DZGP(2,N)/
     &          DZGP(1,N)
            ENDIF
            PLX = PGX
!
!---  Aqueous Hydraulic Gradient  ---
!
          ELSEIF( ABS(IBCT(IEQW,NB)).EQ.11 ) THEN
            PLX = BCX(2)
            CALL FNHGBL( XSBC(MB),YSBC(MB),ZSBC(MB),XPBC(NB),
     &        YPBC(NB),ZPBC(NB),TSBC(MB),TX,PLX )
!            IF( ND(N).EQ.2450 ) PRINT *,'11: PLX = ',PLX,
!     &        ' BCX(2) = ',BCX(2),
!     &        ' TSBC(MB) = ',TSBC(MB),' TX = ',TX,
!     &        ' MB = ',MB,' NB = ',NB,' XSBC(MB) = ',XSBC(MB),
!     &        ' XPBC(NB) = ',XPBC(NB),' YSBC(MB) = ',YSBC(MB),
!     &        ' YPBC(NB) = ',YPBC(NB),' ZSBC(MB) = ',ZSBC(MB),
!     &        ' ZPBC(NB) = ',ZPBC(NB),' ID = ',ID
!
!---  Aqueous Seepage Face  ---
!
          ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 ) THEN
            PLX = BCX(2)
            CALL FNHGBL( XSBC(MB),YSBC(MB),ZSBC(MB),XPBC(NB),
     &        YPBC(NB),ZPBC(NB),TSBC(MB),TX,PLX )
            PLX = MAX( PGX,PLX )
!
!---  Aqueous Potential Evaporation  ---
!
          ELSEIF( ABS(IBCT(IEQW,NB)).EQ.24 ) THEN
            IF( IBCD(NB).LT.0 ) THEN
              PLX = PLX - ABS(BCX(2))*DB*VISL(M,N)/
     &          (PERM(IBD,N)*RKL(M,N)) + RHOL(M,N)*GB
            ELSE
              PLX = PLX + ABS(BCX(2))*DB*VISL(M,N)/
     &          (PERM(IBD,N)*RKL(M,N)) + RHOL(M,N)*GB
            ENDIF
            HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
            IF( HDGL.GT.BCX(3) ) THEN
              PLX = PGX - BCX(3)*RHORL*GRAV
              IBCT(IEQW,NB) = -24
            ELSE
              IBCT(IEQW,NB) = 24
            ENDIF
!
!---      X-Y-Z Aqueous Hydraulic Gradient
!
          ELSEIF( ABS(IBCT(IEQW,NB)).EQ.44 ) THEN
            PLX = BCX(2)
            CALL FNHGBL( XSBC(MB),YSBC(MB),ZSBC(MB),XPBC(NB),
     &        YPBC(NB),ZPBC(NB),TSBC(MB),TX,PLX )
            PLX = PLX + DXX*BCX(3) + DYX*BCX(4) + DZX*BCX(5)
!
!---      X-Y-Z Aqueous Seepage Face
!
          ELSEIF( ABS(IBCT(IEQW,NB)).EQ.45 ) THEN
            PLX = BCX(2)
            CALL FNHGBL( XSBC(MB),YSBC(MB),ZSBC(MB),XPBC(NB),
     &        YPBC(NB),ZPBC(NB),TSBC(MB),TX,PLX )
            PLX = PLX + DXX*BCX(3) + DYX*BCX(4) + DZX*BCX(5)
            PLX = MAX( PGX,PLX )
          ENDIF
!
!---      Secondary and primary boundary variables  ---
!
          SGRMX = 0.D+0
          MX = 2
          INDX = 1
          CALL KSP_W( N,MX,PGX,PLX,SLB(M,NB),RKLB(M,NB),
     &      ASLX,ASLMINX,ASGTX,SGRMX,INDX,IPH(2,N),PGB(1,NB),
     &      PLB(1,NB),SLB(1,NB) )
          SGB(M,NB) = MAX( 1.D+0-SLB(M,NB),ZERO )
          PX = MAX( PGX,PLX )+PATM
          CALL PORSTY_W(N,PX,PCMP(N),PORDB(M,NB),PORTB(M,NB))
          IF( ISLC(3).EQ.1 ) CALL TORTU_W( N,SLB(M,NB),
     &      PORDB(M,NB),TORLB(M,NB) )
!
!---  Convert pressure to absolute prior to computing physical
!     properties  ---
!
          PLX = PLX + PATM
          PGX = PGX + PATM
!
!---  Compute thermodynamic properties  ---
!
          CALL WATSP( T(2,N),PVWB(M,NB) )
          PX = MAX( PLX,PGX,PVWB(M,NB) )
          CALL WATLQD( T(2,N),PX,RHOLB(M,NB) )
          CALL WATLQV( T(2,N),PX,PVWB(2,NB),VISLB(M,NB) )
!
!---  Correct aqueous liquid density and viscosity for electrolyte
!     solute concentration  ---
!
          XLWB(M,NB) = 1.D+0
          IF( ISLC(16).EQ.1 ) THEN
            XVLB = SLB(2,NB)*PORDB(2,NB)
            XVSB = SLB(2,NB)*RHOS(N)*PCSL(1,N,NSL_ELC)*
     &        (1.D+0-PORTB(2,NB))
            CLX = CB(NB,NSL_ELC)/(XVSB+XVLB)
            XLWB(M,NB) = RHOLB(M,NB)
            CALL ELC_DEN( RHOLB(M,NB),CLX,ELC_DCF )
            XLWB(M,NB) = XLWB(M,NB)/RHOLB(M,NB)
            CALL ELC_VIS( VISLB(M,NB),CLX,ELC_VCF )
          ENDIF
          PLB(M,NB) = PLX - PATM
          PGB(M,NB) = PGX - PATM
          TB(M,NB) = TX
        ENDDO
      ENDDO L1
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCP_W group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CAP_W( N,SLX,SGTX,CPGL,SLOX,CPGLO,IPHX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Compute the gas/aqueous capillary pressure from the aqueous
!     saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
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
      SUB_LOG(ISUB_LOG) = '/CAP_W'
      IF( SLX.LT.EPSL ) THEN
        CPGL = SCHR(12,N)*RHORL*GRAV/BTGLX
!
!---    Reset subroutine string sequence  ---
!
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  van Genuchten saturation function w/o gas entrapment,
!     w/ or w/o Webb extension  ---
!
      IF( ISCHR(N).EQ.1 ) THEN
!
!---    van Genuchten saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL  ---
!
        SLRX = SCHR(4,N)
        IF( ISM(N).EQ.2 ) THEN
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          CN = MAX( SCHR(3,N),SMALL )
          IF( SCHR(14,N).LE.ZERO ) THEN
            IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
              CM = 1.D+0 - 2.D+0/CN
            ELSE
              CM = 1.D+0 - 1.D+0/CN
            ENDIF
          ELSE
            CM = SCHR(14,N)
          ENDIF
          SMPX = SCHR(8,N)
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( SLX.LT.SMPX ) THEN
            HMPX = SCHR(9,N)
            DMPX = -(LOG10(SCHR(12,N))-LOG10(HMPX))/SMPX
            HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use van Genuchten function
!
          ELSE
            HDGL = (((1.D+0/ESLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))/
     &        SCHR(1,N)
          ENDIF
!
!---    van Genuchten saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ------
!
        ELSE
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          CN = MAX( SCHR(3,N),SMALL )
          IF( SCHR(14,N).LE.ZERO ) THEN
            IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
              CM = 1.D+0 - 2.D+0/CN
            ELSE
              CM = 1.D+0 - 1.D+0/CN
            ENDIF
          ELSE
            CM = SCHR(14,N)
          ENDIF
          IF( SLX.GT.SCHR(4,N) ) THEN
            HDGL = (((1.D+0/ESLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))/
     &        SCHR(1,N)
          ELSE
            HDGL = SCHR(12,N)
            SLX = SLRX + 1.D-6
          ENDIF
        ENDIF
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Brooks and Corey saturation function w/o gas entrapment,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(N).EQ.2 ) THEN
!
!---    Brooks and Corey saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL  ---
!
        SLRX = SCHR(4,N)
        IF( ISM(N).EQ.2 ) THEN
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          CL = MAX( SCHR(3,N),SMALL )
          SMPX = SCHR(8,N)
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( SLX.LT.SMPX ) THEN
            HMPX = SCHR(9,N)
            DMPX = -(LOG10(SCHR(12,N))-LOG10(HMPX))/SMPX
            HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use Brooks and Corey function
!
          ELSE
            HDGL = SCHR(1,N)*(1.D+0/ESLX)**(1.D+0/CL)
          ENDIF
!
!---    Brooks and Corey saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ---
!
        ELSE
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          IF( (1.D+0-SLX)/EPSL.LT.EPSL ) THEN
            HDGL = SCHR(1,N)
          ELSE
            CL = MAX( SCHR(3,N),SMALL )
            IF( SLX.GT.SCHR(4,N) ) THEN
              HDGL = SCHR(1,N)*(1.D+0/ESLX)**(1.D+0/CL)
            ELSE
              HDGL = SCHR(12,N)
              SLX = SLRX + 1.D-9
            ENDIF
          ENDIF
        ENDIF
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Dual porosity van Genuchten saturation function  ---
!
      ELSEIF( ISCHR(N).EQ.3 ) THEN
        CNM = MAX( SCHR(3,N),SMALL )
        IF( SCHR(14,N).LE.ZERO ) THEN
          IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
            CMM = 1.D+0 - 2.D+0/CNM
          ELSE
            CMM = 1.D+0 - 1.D+0/CNM
          ENDIF
        ELSE
          CMM = SCHR(14,N)
        ENDIF
        CNF = MAX( SCHR(6,N),SMALL )
        IF( SCHR(15,N).LE.EPSL ) THEN
          IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
            CMF = 1.D+0 - 2.D+0/CNF
          ELSE
            CMF = 1.D+0 - 1.D+0/CNF
          ENDIF
        ELSE
          CMF = SCHR(15,N)
        ENDIF
        PORD_MX = (1.D+0-POR(4,N))*POR0(2,N)/
     &    ( POR(4,N) + (1.D+0-POR(4,N))*POR0(2,N) + SMALL )
        PORD_FX = POR(4,N)/
     &    ( POR(4,N) + (1.D+0-POR(4,N))*POR0(2,N) + SMALL )
!
!---    Use matrix properties to generate a guess for 
!       capillary head  ---
!
        IF( SLX.GT.SCHR(4,N) ) THEN
          ASLX = (SLX-SCHR(4,N))/(1.D+0-SCHR(4,N))
          HDGL = (((1.D+0/ASLX)**(1.D+0/CMM)-1.D+0)**(1.D+0/CNM))/
     &      SCHR(1,N)
        ELSE
          HDGL = SCHR(12,N)
          SLX = SCHR(4,N) + 1.D-6
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
        DO
          NC = NC + 1
          REALX = REAL(ISM(N))
          HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,N)),ZERO )*REALX
!
!---      Matrix saturation and partial derivative  ---
!
          SLRX = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
          DSLRX = -SCHR(4,N)/(HDGL*LOG(SCHR(12,N)))*REALX
          ASLX = 1.D+0/((1.D+0 + (SCHR(1,N)*HDGL)**CNM)**CMM)
          DASLX = -CMM*SCHR(1,N)*CNM*((SCHR(1,N)*HDGL)**(CNM-1.D+0))
     &    /((1.D+0 + (SCHR(1,N)*HDGL)**CNM)**(CMM+1.D+0))
          SLMZ = ASLX*(1.D+0-SLRX) + SLRX
          DSLMZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
!
!---      Fracture saturation and partial derivative  ---
!
          SLRX = MAX( (1.D+0-HSCL)*SCHR(7,N),ZERO )
          DSLRX = -SCHR(7,N)/(HDGL*LOG(SCHR(12,N)))*REALX
          ASLX = 1.D+0/((1.D+0 + (SCHR(5,N)*HDGL)**CNF)**CMF)
          DASLX = -CMF*SCHR(5,N)*CNF*((SCHR(5,N)*HDGL)**(CNF-1.D+0))
     &    /((1.D+0 + (SCHR(5,N)*HDGL)**CNF)**(CMF+1.D+0))
          SLFZ = ASLX*(1.D+0-SLRX) + SLRX
          DSLFZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
          F = SLX - SLMZ*PORD_MX - SLFZ*PORD_FX
          DF = -DSLMZ*PORD_MX -DSLFZ*PORD_FX
          DH = -F/(DF+SMALL)
          HDGL = HDGL + DH
!
!---      No convergence on dual porosity van Genuchten 
!         capillary pressure  ---
!
          IF( NC.GT.32 ) THEN
            M_ERR(1) = 'Dual Porosity van Genuchten: '
     &        // 'No Convergence on Capillary Pressure: ' //
     &        'Aqueous Saturation = '
            M_ERR(2) = ' at Node: '
            CALL PATH
            R_ERR = SLX
            I_ERR(1) = ND(N)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
!
!---        Use matrix properties to generate a guess for 
!           capillary head  ---
!
            IF( SLX.GT.SCHR(4,N) ) THEN
              ASLX = (SLX-SCHR(4,N))/(1.D+0-SCHR(4,N))
              HDGL = (((1.D+0/ASLX)**(1.D+0/CMM)-1.D+0)**(1.D+0/CNM))/
     &          SCHR(1,N)
            ELSE
              HDGL = SCHR(12,N)
              SLX = SCHR(4,N) + 1.D-6
            ENDIF
            DH = 0.D+0
          ENDIF
          IF( ABS(DH).LE.1.D-7 ) EXIT
        ENDDO
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Dual porosity Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(N).EQ.4 ) THEN
        ISKIP = 0
        IF( (1.D+0-SLX)/EPSL.LT.EPSL ) THEN
          HDGL = SCHR(1,N)
          ISKIP = 1
        ENDIF
        IF( ISKIP.EQ.0 ) THEN
          CLM = MAX( SCHR(3,N),SMALL )
          CLF = MAX( SCHR(6,N),SMALL )
          PORD_MX = (1.D+0-POR(4,N))*POR0(2,N)/
     &      ( POR(4,N) + (1.D+0-POR(4,N))*POR0(2,N) + SMALL )
          PORD_FX = POR(4,N)/
     &      ( POR(4,N) + (1.D+0-POR(4,N))*POR0(2,N) + SMALL )
!
!---      Use matrix properties to generate a guess for 
!         capillary head  ---
!
          IF( SLX.GT.SCHR(4,N) ) THEN
            ASLX = (SLX-SCHR(4,N))/(1.D+0-SCHR(4,N))
            HDGL = SCHR(1,N)*(1.D+0/ASLX)**(1.D+0/CLM)
          ELSE
            HDGL = SCHR(12,N)
            SLX = SCHR(4,N) + 1.D-9
          ENDIF
!
!---      Start Newton-Raphson solution  ---
!
          NC = 0
          DO
            NC = NC + 1
            REALX = REAL(ISM(N))
            HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,N)),ZERO )*REALX
!
!---        Matrix saturation and partial derivative  ---
!
            SLRX = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
            DSLRX = -SCHR(4,N)*REALX/(HDGL*LOG(SCHR(12,N)))
            HDGLX = MAX( SCHR(1,N),HDGL )
            ASLX = (SCHR(1,N)/HDGLX)**CLM
            DASLX = -CLM*(SCHR(1,N)/(HDGLX**2))
     &        *(SCHR(1,N)/HDGLX)**(CLM-1.D+0)
            SLMZ = ASLX*(1.D+0-SLRX) + SLRX
            DSLMZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
!
!---        Fracture saturation and partial derivative  ---
!
            SLRX = MAX( (1.D+0-HSCL)*SCHR(7,N),ZERO )
            DSLRX = -SCHR(7,N)*REALX/(HDGL*LOG(SCHR(12,N)))
            HDGLX = MAX( SCHR(5,N),HDGL )
            ASLX = (SCHR(5,N)/HDGLX)**CLF
            DASLX = -CLF*(SCHR(5,N)/(HDGLX**2))
     &        *(SCHR(5,N)/HDGLX)**(CLF-1.D+0)
            SLFZ = ASLX*(1.D+0-SLRX) + SLRX
            DSLFZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
            F = SLX - SLMZ*PORD_MX - SLFZ*PORD_FX
            DF = -DSLMZ*PORD_MX -DSLFZ*PORD_FX
            DH = -F/(DF+SMALL)
            HDGL = HDGL + DH
!
!---        No convergence Brooks-Corey capillary pressure  ---
!
            IF( NC.GT.32 ) THEN
              M_ERR(1) = 'Dual Porosity Brooks and Corey: '
     &          // 'No Convergence on Capillary Pressure: ' //
     &          'Aqueous Saturation = '
              M_ERR(2) = ' at Node: '
              CALL PATH
              R_ERR = SLX
              I_ERR(1) = ND(N)
              I_ERR(2) = 1
              I_ERR(3) = 2
              I_ERR(4) = ID
!
!---          Use matrix properties to generate a guess for 
!             capillary head  ---
!
              IF( SLX.GT.SCHR(4,N) ) THEN
                ASLX = (SLX-SCHR(4,N))/(1.D+0-SCHR(4,N))
                HDGL = SCHR(1,N)*(1.D+0/ASLX)**(1.D+0/CLM)
              ELSE
                HDGL = SCHR(12,N)
                SLX = SCHR(4,N) + 1.D-9
              ENDIF
              DH = 0.D+0
            ENDIF
            IF( ABS(DH).LE.1.D-7 ) EXIT
          ENDDO
        ENDIF
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Haverkamp saturation function  ---
!
      ELSEIF( ISCHR(N).EQ.5 ) THEN
        IF( SLX.GT.SCHR(4,N) ) THEN
          ASLX = (SLX-SCHR(4,N))/(1.D+0-SCHR(4,N))
          HDGL = SCHR(1,N) + SCHR(5,N)*
     &      ((SCHR(2,N)/ASLX)-SCHR(2,N))**(1.D+0/SCHR(3,N))
        ELSE
          HDGL = SCHR(12,N)
          SLX = SCHR(4,N) + 1.D-9
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
        DO
          NC = NC + 1
          REALX = REAL(ISM(N))
          HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,N)),ZERO )*REALX
          SLRX = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
          DSLRX = -SCHR(4,N)*REALX/(HDGL*LOG(SCHR(12,N)))
          ASLX = SCHR(2,N)/(SCHR(2,N)+((HDGL-SCHR(1,N))/
     &      SCHR(5,N))**SCHR(3,N))
          DASLX = -(SCHR(2,N)*SCHR(3,N)*
     &      (((HDGL-SCHR(1,N))/SCHR(5,N))**(SCHR(3,N)-1.D+0))
     &      /SCHR(5,N))/((SCHR(2,N)+((HDGL-SCHR(1,N))/SCHR(5,N))
     &      **SCHR(3,N))**2)
          SLZ = ASLX*(1.D+0-SLRX) + SLRX
          DSLZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
          F = SLX - SLZ
          DF = -DSLZ
          DH = -F/(DF+SMALL)
          HDGL = HDGL + DH
!
!---      No convergence Haverkamp capillary pressure  ---
!
          IF( NC.GT.32 ) THEN
            M_ERR(1) = 'Haverkamp: '
     &        // 'No Convergence on Capillary Pressure: ' //
     &        'Aqueous Saturation = '
            M_ERR(2) = ' at Node: '
            CALL PATH
            R_ERR = SLX
            I_ERR(1) = ND(N)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
!
!---        Use matrix properties to generate a guess for 
!           capillary head  ---
!
            IF( SLX.GT.SCHR(4,N) ) THEN
              ASLX = (SLX-SCHR(4,N))/(1.D+0-SCHR(4,N))
              HDGL = SCHR(1,N) + SCHR(5,N)*
     &          ((SCHR(2,N)/ASLX)-SCHR(2,N))**(1.D+0/SCHR(3,N))
            ELSE
              HDGL = SCHR(12,N)
              SLX = SCHR(4,N) + 1.D-9
            ENDIF
            DH = 0.D+0
          ENDIF
          IF( ABS(DH).LE.1.D-7 ) EXIT
        ENDDO
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Linear interpolation function  ---
!
      ELSEIF( ISCHR(N).EQ.10 ) THEN
        ITBX = 0
        HDGL = FNTBLX( SLX,ISLTBL(1,N),ISLTBL(2,N),ITBX )
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Log-linear interpolation function  ---
!
      ELSEIF( ISCHR(N).EQ.11 ) THEN
        ITBX = 0
        HDGL = FNTBLX( SLX,ISLTBL(1,N),ISLTBL(2,N),ITBX )
        HDGL = EXP(HDGL)
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Hysteretic linear interpolation function, using the drainage
!     curve as an initial guess  ---
!
      ELSEIF( ISCHR(N).EQ.12 ) THEN
        SLDX = MAX( SLX-SGTX,0.D+0 )
        ITBX = 0
        HDGL = FNTBLX( SLDX,ISLTBL(1,N),ISLTBL(2,N),ITBX )
        CPGL = HDGL*RHORL*GRAV/BTGLX
        
!
!---  Hysteretic log-linear interpolation function, using the drainage
!     curve  ---
!
      ELSEIF( ISCHR(N).EQ.13 ) THEN
        SLDX = MAX( SLX-SGTX,0.D+0 )
        ITBX = 0
        HDGL = FNTBLX( SLDX,ISLTBL(1,N),ISLTBL(2,N),ITBX )
        HDGL = EXP(HDGL)
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  van Genuchten saturation function w/ gas entrapment,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(N).EQ.101 ) THEN
!
!---    van Genuchten saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        SLRX = SCHR(4,N)
        IF( ISM(N).EQ.2 ) THEN
          ESGTMX = SCHR(15,N)
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
          ESGTMX = SCHR(15,N)/(1.D+0-SLRX)
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
        CN = MAX( SCHR(3,N),SMALL )
        IF( SCHR(14,N).LE.ZERO ) THEN
          IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(14,N)
        ENDIF
        SMPX = SCHR(8,N)
!
!---    van Genuchten saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(N).EQ.2 ) THEN
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( SLX.LT.SMPX ) THEN
            HMPX = SCHR(9,N)
            DMPX = -(LOG10(SCHR(12,N))-LOG10(HMPX))/SMPX
            HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use van Genuchten function
!
          ELSE
            HDGL = (((1.D+0/ESLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))/
     &        SCHR(1,N)
          ENDIF
!
!---    van Genuchten saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          HDGL = (((1.D+0/ASLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))/
     &      SCHR(1,N)
        ENDIF
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Brooks and Corey saturation function w/ gas entrapment,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(N).EQ.102 ) THEN
!
!---    Brooks and Corey saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        SLRX = SCHR(4,N)
        IF( ISM(N).EQ.2 ) THEN
          ESGTMX = SCHR(15,N)
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
          ESGTMX = SCHR(15,N)/(1.D+0-SLRX)
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
        CL = MAX( SCHR(3,N),SMALL )
        SMPX = SCHR(8,N)
!
!---    Brooks and Corey saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(N).EQ.2 ) THEN
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( SLX.LT.SMPX ) THEN
            HMPX = SCHR(9,N)
            DMPX = -(LOG10(SCHR(12,N))-LOG10(HMPX))/SMPX
            HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use Brooks and Corey function
!
          ELSE
            HDGL = SCHR(1,N)*(1.D+0/ESLX)**(1.D+0/CL)
          ENDIF
!
!---    Brooks and Corey saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          HDGL = SCHR(1,N)*(1.D+0/ASLX)**(1.D+0/CL)
        ENDIF
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  van Genuchten drainage-imbibition saturation function,
!     w/ or w/o extensions---
!
      ELSEIF( ISCHR(N).EQ.201 ) THEN
!
!---    van Genuchten drainage-imbibition saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        SLRX = SCHR(4,N)
        IF( ISM(N).EQ.2 ) THEN
          ESGTX = SGTX
          ASLX = SLX + SGTX
          ESLX = (ASLX-SLRX)/(1.D+0-SLRX)
          ESGTMX = SCHR(15,N)
!
!---    van Genuchten drainage-imbibition saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          ESGTX = SGTX/(1.D+0-SLRX)
          ASLX = (SLX+SGTX-SLRX)/(1.D+0-SLRX)
          ESLX = ASLX
          ESGTMX = SCHR(15,N)/(1.D+0-SLRX)
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
        CND = MAX( SCHR(3,N),SMALL )
        IF( SCHR(14,N).LE.ZERO ) THEN
          IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
            CMD = 1.D+0 - 2.D+0/CND
          ELSE
            CMD = 1.D+0 - 1.D+0/CND
          ENDIF
        ELSE
          CMD = SCHR(14,N)
        ENDIF
        SMPDX = SCHR(8,N)
!
!---    Aqueous saturation below the matching point,
!       use Webb extension  ---
!
        IF( SLX.LT.SMPDX ) THEN
          HMPDX = SCHR(9,N)
          DMPDX = -(LOG10(SCHR(12,N))-LOG10(HMPDX))/SMPDX
          HDGLD = 1.D+1**(DMPDX*(ASLX-SMPDX) + LOG10(HMPDX))
!
!---    Aqueous saturation at or above the matching point,
!       use van Genuchten function
!
        ELSE
          HDGLD = (((1.D+0/ESLX)**(1.D+0/CMD)-1.D+0)**(1.D+0/CND))/
     &      SCHR(1,N)
        ENDIF
        CNI = MAX( SCHR(3,N),SMALL )
        IF( SCHR(14,N).LE.ZERO ) THEN
          IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
            CMI = 1.D+0 - 2.D+0/CNI
          ELSE
            CMI = 1.D+0 - 1.D+0/CNI
          ENDIF
        ELSE
          CMI = SCHR(14,N)
        ENDIF
        SMPIX = SCHR(8,N)
!
!---    Aqueous saturation below the matching point,
!       use Webb extension  ---
!
        IF( ASLX.LT.SMPIX ) THEN
          HMPIX = SCHR(11,N)
          DMPIX = -(LOG10(SCHR(12,N))-LOG10(HMPIX))/SMPIX
          HDGLI = 1.D+1**(DMPIX*(ASLX-SMPIX) + LOG10(HMPIX))
!
!---    Aqueous saturation at or above the matching point,
!       use van Genuchten function
!
        ELSE
          HDGLI = (((1.D+0/ESLX)**(1.D+0/CMI)-1.D+0)**(1.D+0/CNI))/
     &      SCHR(2,N)
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
        DO
          NC = NC + 1
!
!---      No convergence on van Genuchten drainage-imbibition 
!         saturation function  ---
!
          IF( NC.GT.32 ) THEN
            M_ERR(1) = 'van Genuchten Drainage-Imbibition: '
     &        // 'No Convergence on Capillary Pressure: ' //
     &        'Aqueous Saturation = '
            M_ERR(2) = ' at Node: '
            CALL PATH
            R_ERR = SLX
            I_ERR(1) = ND(N)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
            HDGL = (ASLMINX/ASLX)*HDGLD + (1.D+0-(ASLMINX/ASLX))*HDGLI
            RETURN
          ENDIF
          DO M = 1,2
            HDGLZ = HDGL
            DHDGLZ = SIGN( MAX( 1.D-6*HDGL,1.D-6 ),
     &        (5.D-1*SCHR(12,N) - HDGL) )
            IF( M.EQ.2 ) HDGLZ = HDGL + DHDGLZ
            HMPDZ = SCHR(9,N)
            HMPIZ = SCHR(11,N)
            SLRZ = SCHR(4,N)
!
!---        van Genuchten drainage-imbibition saturation function,
!           w/ Webb extension ASL = SL + SGT  ---
!
            IF( ISM(N).EQ.2 ) THEN
!
!---          Capillary head above the drainage matching point head,
!             use Webb extension  ---
!
              IF( HDGLZ.GT.HMPDZ ) THEN
                SMPDZ = SCHR(8,N)
                HDGLZ = MIN( HDGLZ,SCHR(12,N) )
                DMPDZ = SMPDZ/(LOG10(SCHR(12,N))-LOG10(HMPDZ))
                SLDZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPDZ
                ASLDZ = SLDZ
                ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
!
!---            Capillary head above the imbibition matching point head,
!               use Webb extension  ---
!
                IF( HDGLZ.GT.HMPIZ ) THEN
                  SMPIZ = SCHR(10,N)
                  DMPIZ = SMPIZ/(LOG10(SCHR(12,N))-LOG10(HMPIZ))
                  SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPIZ
                  ASLIZ = SLIZ
                  IF( ASLDZ.GE.ASLIZ ) THEN
                    ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ELSE
                    ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ENDIF
!
!---            Capillary head at or below the imbibition matching 
!               point head, use van Genuchten function
!
                ELSE
                  CNI = MAX( SCHR(5,N),SMALL )
                  IF( SCHR(13,N).LE.ZERO ) THEN
                    IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
                      CMI = 1.D+0 - 2.D+0/CNI
                    ELSE
                       CMI = 1.D+0 - 1.D+0/CNI
                    ENDIF
                  ELSE
                    CMI = SCHR(13,N)
                  ENDIF
                  ASLIZ = (1.D+0/(1.D+0 + (SCHR(2,N)*HDGLZ)**CNI))**CMI
                  ASLIZ = ASLIZ*(1.D+0-SLRZ) + SLRZ
                  IF( ASLDZ.GE.ASLIZ ) THEN
                    ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ELSE
                    ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ENDIF
                ENDIF
!
!---          Capillary head at or below the drainage matching 
!             point head, use van Genuchten function
!
              ELSE
                CND = MAX( SCHR(3,N),SMALL )
                IF( SCHR(14,N).LE.ZERO ) THEN
                  IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
                    CMD = 1.D+0 - 2.D+0/CND
                  ELSE
                    CMD = 1.D+0 - 1.D+0/CND
                  ENDIF
                ELSE
                  CMD = SCHR(14,N)
                ENDIF
                ASLDZ = (1.D+0/(1.D+0 + (SCHR(1,N)*HDGLZ)**CND))**CMD
                ASLDZ = ASLDZ*(1.D+0-SLRZ) + SLRZ
                ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
!
!---            Capillary head above the imbibition matching point head,
!               use Webb extension  ---
!
                IF( HDGLZ.GT.HMPIZ ) THEN
                  SMPIZ = SCHR(10,N)
                  DMPIZ = SMPIZ/(LOG10(SCHR(12,N))-LOG10(HMPIZ))
                  SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPIZ
                  ASLIZ = SLIZ
                  IF( ASLDZ.GE.ASLIZ ) THEN
                    ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ELSE
                    ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ENDIF
!
!---            Capillary head at or below the imbibition matching 
!               point head, use van Genuchten function
!
                ELSE
                  CNI = MAX( SCHR(5,N),SMALL )
                  IF( SCHR(13,N).LE.ZERO ) THEN
                    IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
                      CMI = 1.D+0 - 2.D+0/CNI
                    ELSE
                      CMI = 1.D+0 - 1.D+0/CNI
                    ENDIF
                  ELSE
                    CMI = SCHR(13,N)
                  ENDIF
                  ASLIZ = (1.D+0/(1.D+0 + (SCHR(2,N)*HDGLZ)**CNI))**CMI
                  ASLIZ = ASLIZ*(1.D+0-SLRZ) + SLRZ
                  IF( ASLDZ.GE.ASLIZ ) THEN
                    ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ELSE
                    ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ENDIF
                ENDIF
              ENDIF
!
!---          Compute trapped gas saturation, using the minimum
!             apparent aqueous saturation  ---
!
              ASLM = MAX( MIN( ASLZ,ASLMINX ),0.D+0 )
              IF( ESGTMX.GT.EPSL .AND. ASLZ.GT.ASLM ) THEN
                SGTMZ = ESGTMX
                R = 1.D+0/SGTMZ - 1.D+0
                ESGTZ = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &            (1.D+0-ASLZ)/(1.D+0 + R*(1.D+0-ASLZ))
                IF( ESGTZ.LT.EPSL ) ESGTZ = 0.D+0
              ELSE
                ESGTZ = 0.D+0
              ENDIF
              SLZ = ASLZ - ESGTZ
              SGTZ = ESGTZ
              SGZ = 1.D+0-SLZ
              IF( SGZ.LT.EPSL ) SGZ = 0.D+0
!
!---        van Genuchten drainage-imbibition saturation function,
!           w/o Webb extension ASL = ESL + ESGT  ---
!
            ELSE
              CND = MAX( SCHR(3,N),SMALL )
              IF( SCHR(14,N).LE.ZERO ) THEN
                IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
                  CMD = 1.D+0 - 2.D+0/CND
                ELSE
                  CMD = 1.D+0 - 1.D+0/CND
                ENDIF
              ELSE
                CMD = SCHR(14,N)
              ENDIF
              ASLDZ = (1.D+0/(1.D+0 + (SCHR(1,N)*HDGLZ)**CND))**CMD
              CNI = MAX( SCHR(5,N),SMALL )
              IF( SCHR(13,N).LE.ZERO ) THEN
                IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
                  CMI = 1.D+0 - 2.D+0/CNI
                ELSE
                  CMI = 1.D+0 - 1.D+0/CNI
                ENDIF
              ELSE
                CMI = SCHR(13,N)
              ENDIF
              ASLIZ = (1.D+0/(1.D+0 + (SCHR(2,N)*HDGLZ)**CNI))**CMI
              ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
              IF( ASLDZ.GE.ASLIZ ) THEN
                ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &             + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
              ELSE
                ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
              ENDIF
!
!---          Compute trapped gas saturation, using the minimum
!             apparent aqueous saturation  ---
!
              ASLM = MAX( MIN( ASLZ,ASLMINX ),0.D+0 )
              IF( ESGTMX.GT.EPSL .AND. ASLZ.GT.ASLM ) THEN
                R = 1.D+0/ESGTMX - 1.D+0
                ESGTZ = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &            (1.D+0-ASLZ)/(1.D+0 + R*(1.D+0-ASLZ))
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
          ENDDO
!
!---      Convergence check  ---
!
          IF( ABS(GX(1)).LT.1.D-9 ) EXIT
!
!---      Solve linear system  ---
!
          FX = GX(1)
          DFX = (GX(2)-GX(1))/DHDGLZ
          DHDGL = -FX/DFX
!
!---      Update primary unknowns  ---
!
          DHX = MAX( 1.D-3,2.D-1*HDGL )
          DHDGL = SIGN( MIN( DHX,ABS(DHDGL) ),DHDGL )
          HDGL = HDGL + DHDGL
          IF( (HDGL-HDGLMN).LT.EPSL .AND. NC.GE.3 ) EXIT
          IF( (HDGLMX-HDGL).LT.EPSL .AND. NC.GE.3 ) EXIT
          HDGL = MIN( HDGL,HDGLMX )
          HDGL = MAX( HDGL,HDGLMN )
        ENDDO
!
!---    Converged solution return capillary pressure  ---
!
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Brooks and Corey drainage-imbibition saturation function,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(N).EQ.202 ) THEN
!
!---    Brooks and Corey drainage-imbibition saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        SLRX = SCHR(4,N)
        IF( ISM(N).EQ.2 ) THEN
          ESGTX = SGTX
          ASLX = SLX + SGTX
          ESGTMX = SCHR(15,N)
!
!---    Brooks and Corey drainage-imbibition saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          ESGTX = SGTX/(1.D+0-SLRX)
          ASLX = (SLX+SGTX-SLRX)/(1.D+0-SLRX)
          ESGTMX = SCHR(15,N)/(1.D+0-SLRX)
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
        CLD = MAX( SCHR(3,N),SMALL )
        SMPDX = SCHR(8,N)
!
!---    Aqueous saturation below the matching point,
!       use Webb extension  ---
!
        IF( ASLX.LT.SMPDX ) THEN
          HMPDX = SCHR(9,N)
          DMPDX = -(LOG10(SCHR(12,N))-LOG10(HMPDX))/SMPDX
          HDGLD = 1.D+1**(DMPDX*(ASLX-SMPDX) + LOG10(HMPDX))
!
!---    Aqueous saturation at or above the matching point,
!       use Brooks and Corey function
!
        ELSE
!
!---      w/ Webb extension  ---
!
          IF( ISM(N).EQ.2 ) THEN
            ASLZ = (ASLX-SLRX)/(1.D+0-SLRX)
            HDGLD = SCHR(1,N)*(1.D+0/ASLZ)**(1.D+0/CLD)
!
!---      w/o Webb extension  ---
!
          ELSE
            HDGLD = SCHR(1,N)*(1.D+0/ASLX)**(1.D+0/CLD)
          ENDIF
        ENDIF
!
!---    Imbibition capillary head  ---
!
        CLI = MAX( SCHR(6,N),SMALL )
        SMPIX = SCHR(10,N)
!
!---    Aqueous saturation below the matching point,
!       use Webb extension  ---
!
        IF( ASLX.LT.SMPIX ) THEN
          HMPIX = SCHR(11,N)
          DMPIX = -(LOG10(SCHR(12,N))-LOG10(HMPIX))/SMPIX
          HDGLI = 1.D+1**(DMPIX*(ASLX-SMPIX) + LOG10(HMPIX))
!
!---    Aqueous saturation at or above the matching point,
!       use Brooks and Corey function
!
        ELSE
!
!---      w/ Webb extension  ---
!
          IF( ISM(N).EQ.2 ) THEN
            ASLZ = (ASLX-SLRX)/(1.D+0-SLRX)
            HDGLI = SCHR(5,N)*(1.D+0/ASLZ)**(1.D+0/CLI)
!
!---      w/o Webb extension  ---
!
          ELSE
            HDGLI = SCHR(5,N)*(1.D+0/ASLX)**(1.D+0/CLI)
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
        DO
          NC = NC + 1
!
!---      No convergence on Brooks and Corey drainage-imbibition 
!         saturation function  ---
!
          IF( NC.GT.32 ) THEN
            M_ERR(1) = 'Brooks and Corey Drainage-Imbibition: '
     &        // 'No Convergence on Capillary Pressure: ' //
     &        'Aqueous Saturation = '
            M_ERR(2) = ' at Node: '
            CALL PATH
            R_ERR = SLX
            I_ERR(1) = ND(N)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
            HDGL = (ASLMINX/ASLX)*HDGLD + (1.D+0-(ASLMINX/ASLX))*HDGLI
            RETURN
          ENDIF
          DO M = 1,2
            HDGLZ = HDGL
            DHDGLZ = SIGN( MAX( 1.D-6*HDGL,1.D-6 ),
     &        (5.D-1*SCHR(12,N) - HDGL) )
            IF( M.EQ.2 ) HDGLZ = HDGL + DHDGLZ
            HMPDZ = SCHR(9,N)
            HMPIZ = SCHR(11,N)
            SLRZ = SCHR(4,N)
!
!---        Brooks and Corey drainage-imbibition saturation function,
!           w/ Webb extension ASL = SL + SGT  ---
!
            IF( ISM(N).EQ.2 ) THEN
!
!---          Capillary head above the drainage matching point head,
!             use Webb extension  ---
!
              IF( HDGLZ.GT.HMPDZ ) THEN
                SMPDZ = SCHR(8,N)
                HDGLZ = MIN( HDGLZ,SCHR(12,N) )
                DMPDZ = SMPDZ/(LOG10(SCHR(12,N))-LOG10(HMPDZ))
                SLDZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPDZ
                ASLDZ = SLDZ
                ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
!
!---            Capillary head above the imbibition matching point head,
!               use Webb extension  ---
!
                IF( HDGLZ.GT.HMPIZ ) THEN
                  SMPIZ = SCHR(10,N)
                  DMPIZ = SMPIZ/(LOG10(SCHR(12,N))-LOG10(HMPIZ))
                  SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPIZ
                  ASLIZ = SLIZ
                  IF( ASLDZ.GE.ASLIZ ) THEN
                    ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ELSE
                    ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ENDIF
!
!---            Capillary head at or below the imbibition matching 
!               point head, use Brooks and Corey function
!
                ELSE
                  CLI = MAX( SCHR(6,N),SMALL )
                  IF( HDGLZ.LE.SCHR(5,N) ) THEN
                    ASLIZ = 1.D+0
                  ELSE
                    ASLIZ = (SCHR(5,N)/HDGLZ)**CLI
                  ENDIF
                  ASLIZ = ASLIZ*(1.D+0-SLRZ) + SLRZ
                  IF( ASLDZ.GE.ASLIZ ) THEN
                    ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ELSE
                    ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ENDIF
                ENDIF
!
!---          Capillary head at or below the drainage matching 
!             point head, use Brooks and Corey function
!
              ELSE
                CLD = MAX( SCHR(3,N),SMALL )
                IF( HDGLZ.LE.SCHR(1,N) ) THEN
                  ASLDZ = 1.D+0
                ELSE
                  ASLDZ = (SCHR(1,N)/HDGLZ)**CLD
                ENDIF
                ASLDZ = ASLDZ*(1.D+0-SLRZ) + SLRZ
                ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
!
!---            Capillary head above the imbibition matching point head,
!               use Webb extension  ---
!
                IF( HDGLZ.GT.HMPIZ ) THEN
                  SMPIZ = SCHR(10,N)
                  DMPIZ = SMPIZ/(LOG10(SCHR(12,N))-LOG10(HMPIZ))
                  SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPIZ
                  ASLIZ = SLIZ
                  IF( ASLDZ.GE.ASLIZ ) THEN
                    ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ELSE
                    ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ENDIF
!
!---            Capillary head at or below the imbibition matching 
!               point head, use Brooks and Corey function
!
                ELSE
                  CLI = MAX( SCHR(6,N),SMALL )
                  IF( HDGLZ.LE.SCHR(5,N) ) THEN
                    ASLIZ = 1.D+0
                  ELSE
                    ASLIZ = (SCHR(5,N)/HDGLZ)**CLI
                  ENDIF
                  ASLIZ = ASLIZ*(1.D+0-SLRZ) + SLRZ
                  IF( ASLDZ.GE.ASLIZ ) THEN
                    ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ELSE
                    ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                  ENDIF
                ENDIF
              ENDIF
!
!---          Compute trapped gas saturation, using the minimum
!             apparent aqueous saturation  ---
!
              ASLM = MAX( MIN( ASLZ,ASLMINX ),0.D+0 )
              IF( ESGTMX.GT.EPSL .AND. ASLZ.GT.ASLM ) THEN
                SGTMZ = ESGTMX
                R = 1.D+0/SGTMZ - 1.D+0
                ESGTZ = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &            (1.D+0-ASLZ)/(1.D+0 + R*(1.D+0-ASLZ))
                IF( ESGTZ.LT.EPSL ) ESGTZ = 0.D+0
              ELSE
                ESGTZ = 0.D+0
              ENDIF
              SLZ = ASLZ - ESGTZ
              SGTZ = ESGTZ
              SGZ = 1.D+0-SLZ
              IF( SGZ.LT.EPSL ) SGZ = 0.D+0
!
!---        Brooks and Corey drainage-imbibition saturation function,
!           w/o Webb extension ASL = ESL + ESGT  ---
!
            ELSE
              CLD = MAX( SCHR(3,N),SMALL )
              IF( HDGLZ.LE.SCHR(1,N) ) THEN
                ASLDZ = 1.D+0
              ELSE
                ASLDZ = (SCHR(1,N)/HDGLZ)**CLD
              ENDIF
              CLI = MAX( SCHR(6,N),SMALL )
              IF( HDGLZ.LE.SCHR(5,N) ) THEN
                ASLIZ = 1.D+0
              ELSE
                ASLIZ = (SCHR(5,N)/HDGLZ)**CLI
              ENDIF
              ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
              IF( ASLDZ.GE.ASLIZ ) THEN
                ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &             + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
              ELSE
                ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &                 + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
              ENDIF
!
!---          Compute trapped gas saturation, using the minimum
!             apparent aqueous saturation  ---
!
              ASLM = MAX( MIN( ASLZ,ASLMINX ),0.D+0 )
              IF( ESGTMX.GT.EPSL .AND. ASLZ.GT.ASLM ) THEN
                R = 1.D+0/ESGTMX - 1.D+0
                ESGTZ = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &            (1.D+0-ASLZ)/(1.D+0 + R*(1.D+0-ASLZ))
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
          ENDDO
!
!---      Convergence check  ---
!
          IF( ABS(GX(1)).LT.1.D-9 ) EXIT
!
!---      Solve linear system  ---
!
          FX = GX(1)
          DFX = (GX(2)-GX(1))/DHDGLZ
          DHDGL = -FX/DFX
!
!---      Update primary unknowns  ---
!
          DHX = MAX( 1.D-3,2.D-1*HDGL )
          DHDGL = SIGN( MIN( DHX,ABS(DHDGL) ),DHDGL )
          IF( (HDGL-HDGLMN).LT.EPSL .AND. NC.GE.3 ) EXIT
          IF( (HDGLMX-HDGL).LT.EPSL .AND. NC.GE.3 ) EXIT
          HDGL = HDGL + DHDGL
          HDGL = MIN( HDGL,HDGLMX )
          HDGL = MAX( HDGL,HDGLMN )
        ENDDO
!
!---    Converged solution return capillary pressure  ---
!
        CPGL = HDGL*RHORL*GRAV/BTGLX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ELC_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Correct aqueous liquid density and viscosity for electrolyte
!     solute concentration
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE JACOB
      USE HYST
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
      IF( ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ELC_W'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = ND(N)
        PGX = PG(2,N) + PATM
        DO M = 2,ISVC+2
          PLX = PL(M,N) + PATM
          CALL WATSP( T(2,N),PVW(M,N) )
          PX = MAX( PLX,PGX,PVW(M,N) )
          CALL WATLQD( T(2,N),PX,RHOL(M,N) )
          CALL WATLQV( T(2,N),PX,PVW(M,N),VISL(M,N) )
          CLX = C(N,NSL_ELC)*YL(N,NSL_ELC)/(SL(M,N)*PORD(M,N)+SMALL)
          XLW(M,N) = RHOL(M,N)
          CALL ELC_DEN( RHOL(M,N),CLX,ELC_DCF )
          XLW(M,N) = XLW(M,N)/RHOL(M,N)
          CALL ELC_VIS( VISL(M,N),CLX,ELC_VCF )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ELC_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLUX_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Compute fluxes on internal and boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE SOLTN
      USE MPI
      USE HYST
      USE GLB_PAR
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
      SUB_LOG(ISUB_LOG) = '/FLUX_W'
!
!---  Aqueous volumetric flux (non-boundary surfaces)  ---
!
      CALL DRCVL_W
!
!---  Water aqueous volumetric flux (non-boundary surfaces)  ---
!
      CALL DFFLW_W
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLUX_W group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Compute primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE SOLTN
      USE MPI
      USE PROP
      USE HYST
      USE GLB_PAR
      USE GRID
      USE FLUX
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
      IF( ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM_W'
!
!---  Determine wetting/drying direction on the first
!     iteration  ---
!
      IF( NITER.EQ.1 ) THEN
        DO N = 1,NFCGC(ID+1)
          IF( IXP(N).EQ.0 ) CYCLE
          N_DB = ND(N)
!
!---      Draining  ---
!
          IF( (PG(2,N)-PL(2,N)).GT.(PG(1,N)-PL(1,N)) ) THEN
            IPH(2,N) = -1
!
!---        Reset reversal point  ---
!
            IF( IPH(1,N).EQ.1 ) THEN
              ASLMIN(2,N) = SL(1,N)
              ASTMIN(2,N) = BTGL(1,N)*(PG(1,N)-PL(1,N))/RHORL/GRAV
            ENDIF
!
!---      Wetting  ---
!
          ELSE
            IPH(2,N) = 1
!
!---        Reset reversal point  ---
!
            IF( IPH(1,N).EQ.-1 ) THEN
              ASLMIN(2,N) = SL(1,N)
              ASTMIN(2,N) = BTGL(1,N)*(PG(1,N)-PL(1,N))/RHORL/GRAV
            ENDIF
          ENDIF
        ENDDO
      ENDIF
!
!---  Compute increments  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = ND(N)
        IF( NITER.EQ.1 ) THEN
          IF( IPH(2,N).EQ.2 ) THEN
            IF( (PG(2,N)-PL(2,N))/RHORL/GRAV.LE.HCMWE(N) ) THEN
              IF( (PG(2,N)-PL(2,N)).GT.(PG(1,N)-PL(1,N)) ) THEN
                IPH(2,N) = -1
              ELSE
                IPH(2,N) = 1
              ENDIF
            ENDIF
          ELSE
            IF( (PG(2,N)-PL(2,N)).GT.(PG(1,N)-PL(1,N)) ) THEN
              IPH(2,N) = -1
            ELSE
              IPH(2,N) = 1
            ENDIF
          ENDIF
        ENDIF
        IF( ISCHR(N).EQ.1 .OR. ISCHR(N).EQ.2 ) THEN
          DNR(IEQW,N) = -MAX( 1.D-1,1.D-6*ABS(PG(2,N)-PL(2,N)) )
        ELSE
          DNR(IEQW,N) = -MAX( 1.D-1,1.D-6*ABS(PG(2,N)-PL(2,N)) )
        ENDIF
!
!--- Increment the primary variables  ---
!
        DO M = 3,ISVC+2
          PL(M,N) = PL(2,N)
          IF( M.EQ.IEQW+2 ) THEN
            PL(M,N) = PL(M,N) + DNR(IEQW,N)
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM_W group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KSP_W( N,M,PGX,PLX,SLX,RKLX,ASLX,ASLMINX,
     &  ASGTX,SGRMX,INDX,IPHX,PGOX,PLOX,SLOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMP-W)
!
!     Compute the aqueous saturation from the gas/aqueous capillary
!     pressure, and compute the aqueous relative permeability from the
!     aqueous saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/KSP_W'
!
!---  Aqueous and gas saturation  ---
!
      CALL SP_W( N,M,PGX,PLX,SLX,RKLX,ASLX,ASLMINX,
     &  ASGTX,SGRMX,HDGL,SLP,SLPF,SLPM,INDX,IPHX,PGOX,PLOX,SLOX )
!
!---  Aqueous relative permeability  ---
!
      CALL RKLS_W( HDGL,RKLX,SLP,SLPF,SLPM,SLX,N,IPHX,M )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KSP_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Load the current time step values into the old time step
!     variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE REACT
      USE HYST
      USE GRID
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/LDO_W'
!
!---  Assign old time step values  ---
!
      DO N = 1,NFCGC(ID+1)
        PG(1,N) = PG(2,N)
        T(1,N) = T(2,N)
        PL(1,N) = PL(2,N)
        PORD(1,N) = PORD(2,N)
        PORT(1,N) = PORT(2,N)
        SL(1,N) = SL(2,N)
        SG(1,N) = SG(2,N)
        PVW(1,N) = PVW(2,N)
        XLW(1,N) = XLW(2,N)
        RHOL(1,N) = RHOL(2,N)
        VISL(1,N) = VISL(2,N)
        TORL(1,N) = TORL(2,N)
        RKL(1,N) = RKL(2,N)
        TRPGL(1,N) = TRPGL(2,N)
        ASLMIN(1,N) = MIN( ASL(N),ASLMIN(2,N) )
        ASLMIN(2,N) = ASLMIN(1,N)
        NPHAZ(1,N) = NPHAZ(2,N)
        IPH(1,N) = IPH(2,N)
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
      DO NB = 1,NBC(ID+1)
        PLB(1,NB) = PLB(2,NB)
        PGB(1,NB) = PGB(2,NB)
        SLB(1,NB) = SLB(2,NB)
      ENDDO








!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO_W group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Compute hydrologic, thermodynamic and physical properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE JACOB
      USE HYST
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
      IF( ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROP_W'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = ND(N)
!
!---    ECKEChem  ---
!
        POR0(1,N) = POR0(1,N)
        POR0(2,N) = POR0(2,N)
        PGX = PG(2,N) + PATM
        DO M = 2,ISVC+2
!
!---      Saturation, relative permeability, porosity, and
!         tortuosity calculations  ---
!
          SGRMX = SCHR(15,N)/(1.D+0+TRPGL(2,N)/SCHR(9,N))
          ASLMINX = ASLMIN(1,N)
          CALL KSP_W( N,M,PG(M,N),PL(M,N),SL(M,N),RKL(M,N),
     &      ASLX,ASLMINX,ASGTX,SGRMX,INDX,IPH(2,N),
     &      PG(1,N),PL(1,N),SL(1,N) )
          IF( M.EQ.2 ) THEN
            ASL(N) = ASLX
            ASGT(N) = ASGTX
            ASLMIN(2,N) = ASLMINX
          ENDIF
          SGT(M,N) = ASGTX*(1.D+0-SCHR(4,N))
          SG(M,N) = MAX( 1.D+0-SL(M,N),ZERO )
          PX = PL(M,N) + PATM
          CALL PORSTY_W(N,PX,PCMP(N),PORD(M,N),PORT(M,N))
          IF( ISLC(3).EQ.1 ) CALL TORTU_W( N,SL(M,N),
     &      PORD(M,N),TORL(M,N) )
!
!---      Kozeny-Carman Permeability Model ---
!
          PERMRF(M,N) = 1.D+0
          IF( IPRF(N).EQ.2 ) CALL PERM_I( PERMRF(M,N),PORD(M,N) )
!
!---      Poroelastic Permeability Model ---
!
          IF( IPRF(N)/10.EQ.1 ) THEN
            DPX = PX-PCMP(N)
            PERMRF(M,N) = PERMRF(M,N)*EXP( PERM(4,N)*DPX/PERM(5,N) )
          ENDIF
!
!---      Aqueous density and viscosity  ---
!
          PLX = PL(M,N) + PATM
          CALL WATSP( T(2,N),PVW(M,N) )
          PX = MAX( PLX,PGX,PVW(M,N) )
          CALL WATLQD( T(2,N),PX,RHOL(M,N) )
          CALL WATLQV( T(2,N),PX,PVW(M,N),VISL(M,N) )
          XLW(M,N) = 1.D+0
!
!---      Correct aqueous liquid density and viscosity for electrolyte
!         solute concentration  ---
!
          IF( ISLC(16).EQ.1 ) THEN
            CLX = C(N,NSL_ELC)*YL(N,NSL_ELC)/(SL(M,N)*PORD(M,N)+SMALL)
            XLW(M,N) = RHOL(M,N)
            CALL ELC_DEN( RHOL(M,N),CLX,ELC_DCF )
            XLW(M,N) = XLW(M,N)/RHOL(M,N)
            CALL ELC_VIS( VISL(M,N),CLX,ELC_VCF )
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_BIN_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Read binary files from processor.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE REACT
      USE GRID
      USE GLB_PAR
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
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_BIN_W'
!
!---  Read input.bin for solution control and output data  ---
!
      CALL READ_SOLU_W
!      PRINT *,'Post READ_SOLU_W: ID = ',ID
!
!---  Read input.bin for grid data  ---
!
      CALL READ_GRID_W
!      PRINT *,'Post READ_GRID_W: ID = ',ID
!
!---  Read input.bin for property data  ---
!
      CALL READ_PROP_W
!      PRINT *,'Post READ_PROP_W: ID = ',ID
!
!---  Read input.bin for state condition data  ---
!
      CALL READ_STATE_W
!      PRINT *,'Post READ_STATE_W: ID = ',ID
!
!---  Read input.bin for boundary condition data  ---
!
      CALL READ_BOCO_W
!      PRINT *,'Post READ_BOCO_W: ID = ',ID
!
!---  Read input.bin for source data  ---
!
      CALL READ_SORC_W
!      PRINT *,'Post READ_SORC_W: ID = ',ID
!
!---  Read well.bin for coupled-well data ---
!
      CALL READ_COUP_WELL_W
!      PRINT *,'Post READ_COUP_WELL_W: ID = ',ID
!
!---  Read input.bin for transport data  ---
!
      IF( IEQC.NE.0 ) CALL READ_TPORT_W
!      PRINT *,'Post READ_TPORT_W: ID = ',ID
!
!---  Read gmec.bin for geomechanical property data  ---
!
      IF( ISLC(50).NE.0 ) CALL READ_GMEC
!      PRINT *,'Post READ_GMEC: ID = ',ID
!
!---  Read gmbc.bin for geomechanical boundary condition data  ---
!
      IF( ISLC(50).NE.0 ) CALL READ_GMBC
!      PRINT *,'Post READ_GMBC: ID = ',ID
!
!---  Read input.bin for ECKEChem (i.e., reactive transport data)  ---
!
      IF( ISLC(40).EQ.1 ) CALL READ_REACT_W
!      PRINT *,'Post READ_REACT_CO2: ID = ',ID
!
!---  Read ptzrcoef.bin  ---
!
      IF( IACTV.EQ.2 ) CALL READ_PTZRCOEF
!      PRINT *,'Post READ_PTZR_CO2: IACTV = ',IACTV,'ID = ',ID
!
!---  Read restart.bin for restart simulations  ---
!
      IF( IEO.EQ.2 ) CALL RDRST_W
!      PRINT *,'Post RDRST_W: IEO = ',IEO,'ID = ',ID
!
!---  Check for fatal execution errors and stop simulation
!     if detected  ---
!
      CALL CHK_ERROR
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_BIN_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_BOCO_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Read binary input.bin file for boundary condition data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE SOLTN
      USE PROP
      USE MPI
      USE GRID
      USE GLB_PAR
      USE FILES
      USE BCVP
      USE BCV
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARX
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_BOCO_W'
!
!---  Open boco.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'boco.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Allocate memory for NBC  ---
!
      ALLOCATE( NBC(1:NP),STAT=ISTAT )
      CHMSG = 'NBC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read number of boundary condition surfaces on each
!     processor  ---
!
      NVAR = NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NBC,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Allocate memory for boundary condition arrays  ---
!
      CALL ALLOC_BCV
      CALL ALLOC_BCVP
!
!---  Initialize boundary condition variables  ---
!
      CALL INTLZ_BCV
      CALL INTLZ_BCVP
!
!---  Read boundary condition variables
!     (duplicated across processors)  ---
!
      NVAR = LBCV*LBTM*LBCIN
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,BC,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read boundary condition base temperature
!     (duplicated across processors)  ---
!
      NVAR = LBCIN
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TSBC,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read boundary condition base x-direction surface centroid
!     (duplicated across processors)  ---
!
      NVAR = LBCIN
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,XSBC,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read boundary condition base y-direction surface centroid
!     (duplicated across processors)  ---
!
      NVAR = LBCIN
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,YSBC,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read boundary condition base z-direction surface centroid
!     (duplicated across processors)  ---
!
      NVAR = LBCIN
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ZSBC,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read array of boundary condition nodes  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NBC(I)
      ENDDO
      NCP = 0
      DO N = 1,NP
        NCP = NCP + NBC(N)
      ENDDO
!
!---  Index array of boundary condition field nodes  ---
!
      NVAR = NBC(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCN,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'Post IBCN: ID = ',ID
!
!---  Index array of boundary condition directions  ---
!
      NVAR = NBC(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCD,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'Post IBCD: ID = ',ID
!
!---  Index array of boundary condition number of time points  ---
!
      NVAR = NBC(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'Post IBCM: ID = ',ID
!
!---  Index array of boundary condition input links  ---
!
      NVAR = NBC(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCIN,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'Post IBCIN: ID = ',ID
!
!---  Index array of boundary condition cycling options  ---
!
      NVAR = NBC(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCC,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'Post IBCC: ID = ',ID
!
!---  Index array of boundary condition types  ---
!
      LX = LUK+LSOLU*LC
      NVAR = NBC(ID+1)*LX
      OFFSET = IOFFSET + NC*LX*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*LX*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCT,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'Post IBCT: ID = ',ID
!
!---  Index array of boundary condition reactive species  ---
!
      NVAR = NBC(ID+1)*(LSPBC+1)
      OFFSET = IOFFSET + NC*(LSPBC+1)*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*(LSPBC+1)*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCSP,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'IBCSP = ',((IBCSP(K,L),K=1,(LSPBC+1)),L=1,NBC(ID+1)),
!     &  ' ID = ',ID
!
!---  Array of x-coordinate of boundary surface centroid  ---
!
      NVAR = NBC(ID+1)
      OFFSET = IOFFSET + NC*NBYTR + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,XPBC,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!      PRINT *,'Post XPBC: ID = ',ID
!
!---  Array of y-coordinate of boundary surface centroid  ---
!
      NVAR = NBC(ID+1)
      OFFSET = IOFFSET + NC*NBYTR + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,YPBC,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!      PRINT *,'Post YPBC: ID = ',ID
!
!---  Array of z-coordinate of boundary surface centroid  ---
!
      NVAR = NBC(ID+1)
      OFFSET = IOFFSET + NC*NBYTR + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ZPBC,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!      PRINT *,'Post ZPBC: ID = ',ID
!      DO NB = 1,NBC(ID+1)
!        PRINT *,'XPBC(',NB,') = ',XPBC(NB),' ID = ',ID
!        PRINT *,'YPBC(',NB,') = ',YPBC(NB),' ID = ',ID
!        PRINT *,'ZPBC(',NB,') = ',ZPBC(NB),' ID = ',ID
!      ENDDO
!
!---  Allocate local temporary real initial condition boundary
!     condition array memory  ---
!
      ALLOCATE( VARX(1:NBC(ID+1)),STAT=ISTAT )
      CHMSG = 'VARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Array of initial condition boundary condition aqueous 
!     pressure  ---
!
      NVAR = NBC(ID+1)
      OFFSET = IOFFSET + NC*NBYTR + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NBC(ID+1)
        PLB(1,N) = VARX(N)
        VARX(N) = 0.D+0
      ENDDO
!
!---  Deallocate local temporary real initial condition boundary
!     condition array memory  ---
!
      DEALLOCATE( VARX,STAT=ISTAT )
      CHMSG = 'VARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Close the boco.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_BOCO_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_COUP_WELL_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Read binary well.bin file for coupled-well data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 January 2023
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE SOLTN
      USE PROP
      USE OUTPU
      USE MPI
      USE GRID
      USE GLB_PAR
      USE SOURC
      USE FILES
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
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_COUP_WELL_W'
!
!---  Open well.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'well.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Allocate memory for coupled-well variables  ---
!
      CALL ALLOC_COUP_WELL
!
!---  Initialize coupled-well variables  ---
!
      CALL INTLZ_COUP_WELL
!
!---  Read number of coupled wells (duplicated across processors)  ---
!
      NVAR = 1
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,N_CW,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'N_CW = ',N_CW,' ID = ',ID
!
!---  Read coupled well parameters, indices, and variables
!     if coupled wells are modeled  ---
!
      IF( N_CW.GT.0 ) THEN
!
!---    Read well names (duplicated across processors)  ---
!
        NVAR = 64*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,WNM_CW,NVAR,MPI_CHAR,
     &    STATUS,IERR )
!        PRINT *,'WNM_CW = ',(WNM_CW(M),M=1,LN_CW),' ID = ',ID
!
!---    Read type index of coupled wells
!       (duplicated across processors)  ---
!
        NVAR = LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IT_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'IT_CW = ',(IT_CW(M),M=1,NVAR),' ID = ',ID
!
!---    Read injection well state index
!       (duplicated across processors)  ---
!
        NVAR = LWTP_CW*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ITS_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'ITS_CW = ',((ITS_CW(M,N),M=1,LWTP_CW),N=1,LN_CW),
!       &  ' ID = ',ID
!
!---    Read x-, y-, and z-direction coupled well fraction factors
!       (duplicated across processors)  ---
!
        NVAR = 3*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,FF_CW,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        PRINT *,'FF_CW = ',((FF_CW(M,N),M=1,3),N=1,LN_CW),
!       &  ' ID = ',ID
!
!---    Read coupled well indices (duplicated across processors)  ---
!
        NVAR = 10*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ID_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        DO N = 1,LN_CW
!          DO M = 1,10
!            PRINT *,'ID_CW(',M,',',N,') = ',ID_CW(M,N),'ID = ',ID
!          ENDDO
!        ENDDO
!
!---    Read coupled well total mass limit
!       (duplicated across processors)  ---
!
        NVAR = LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,TML_CW,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        DO N = 1,LN_CW
!          PRINT *,'TML_CW(',N,') = ',TML_CW(N),'ID = ',ID
!        ENDDO
!
!---    Read coupled well cyclic index (duplicated across processors)  ---
!
        NVAR = LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ICC_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        DO N = 1,LN_CW
!          PRINT *,'ICC_CW(',N,') = ',ICC_CW(N),'ID = ',ID
!        ENDDO
!
!---    Read number of coupled well time points
!       (duplicated across processors)  ---
!
        NVAR = LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IM_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'IM_CW = ',(IM_CW(N),N=1,LN_CW),
!       &  ' ID = ',ID
!
!---    Read number of coupled well time points per time period
!       (duplicated across processors)  ---
!
        NVAR = LWTP_CW*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IMP_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'IMP_CW = ',((IMP_CW(M,N),M=1,LWTP_CW),N=1,LN_CW),
!       &  ' ID = ',ID
!
!---    Read number of coupled well reference nodes
!       (duplicated across processors)  ---
!
        NVAR = LVREF
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IREF_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'IREF_CW = ',(IREF_CW(N),N=1,LVREF),
!       &  ' ID = ',ID
!
!---    Read coupled-well field to node pointers
!       (duplicated across processors)  ---
!
        NVAR = LWF_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IWFG_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'IWFG_CW = ',(IWFG_CW(N),N=1,LWF_CW),
!       &  ' ID = ',ID
!
!---    Read coupled-well field to node-equation pointers
!       (duplicated across processors)  ---
!
        NVAR = LWF_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IXPG_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'IXPG_CW = ',(IXPG_CW(N),N=1,LWF_CW),
!       &  ' ID = ',ID
!
!---    Read coupled-well field to node pointers
!       (duplicated across processors)  ---
!
        NVAR = LWN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IWNG_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'IWNG_CW = ',(IWNG_CW(N),N=1,LWN_CW),
!       &  ' ID = ',ID
!
!---    Read coupled-well field to node pointers
!       (duplicated across processors)  ---
!
        NVAR = LWN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IWP_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'IWP_CW = ',(IWP_CW(N),N=1,LWN_CW),
!       &  ' ID = ',ID
!
!---    Read type index for coupled well interval
!       (duplicated across processors)  ---
!
        NVAR = LWI_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IS_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        DO N = 1,LWI_CW
!          PRINT *,'IS_CW(',N,') = ',IS_CW(N),'ID = ',ID
!        ENDDO
!
!---    Read x-transition points for coupled well interval
!       (duplicated across processors)  ---
!
        NVAR = 2*LWI_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,XTP_CW,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        DO N = 1,LWI_CW
!          DO M = 1,2
!            PRINT *,'XTP_CW(',M,',',N,') = ',XTP_CW(M,N),'ID = ',ID
!          ENDDO
!        ENDDO
!
!---    Read y-transition points for coupled well interval
!       (duplicated across processors)  ---
!
        NVAR = 2*LWI_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,YTP_CW,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        DO N = 1,LWI_CW
!          DO M = 1,2
!            PRINT *,'YTP_CW(',M,',',N,') = ',YTP_CW(M,N),'ID = ',ID
!          ENDDO
!        ENDDO
!
!---    Read z-transition points for coupled well interval
!       (duplicated across processors)  ---
!
        NVAR = 2*LWI_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ZTP_CW,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        DO N = 1,LWI_CW
!          DO M = 1,2
!            PRINT *,'ZTP_CW(',M,',',N,') = ',ZTP_CW(M,N),'ID = ',ID
!          ENDDO
!        ENDDO
!
!---    Read parameters for coupled well interval
!       (duplicated across processors)  ---
!
        NVAR = 5*LWI_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,PAR_CW,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        DO N = 1,LWI_CW
!          DO M = 1,5
!            PRINT *,'PAR_CW(',M,',',N,') = ',PAR_CW(M,N),'ID = ',ID
!          ENDDO
!        ENDDO
!
!---    Read coupled-well interval
!       (duplicated across processors)  ---
!
        NVAR = LWN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,INV_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!
!---    Read time dependent driving parameters for the coupled well
!       (duplicated across processors)  ---
!
        NVAR = 8*LWT_CW*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,VAR_CW,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!
!---    Read time dependent solute driving parameters for the coupled
!       well  (duplicated across processors)  ---
!
        NVAR = LSOLU_CW*LWT_CW*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,VARC_CW,NVAR,MPI_REAL8,
     &      STATUS,IERR )
!
!---    Read time dependent species driving parameters for the coupled
!       well (duplicated across processors)  ---
!
        NVAR = LSPC_CW*LWT_CW*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,VARSP_CW,NVAR,MPI_REAL8,
     &      STATUS,IERR )
!
!---    Read x-direction well projections
!       (duplicated across processors)  ---
!
        NVAR = LWN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,PLX_CW,NVAR,MPI_REAL8,
     &      STATUS,IERR )
!
!---    Read y-direction well projections
!       (duplicated across processors)  ---
!
        NVAR = LWN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,PLY_CW,NVAR,MPI_REAL8,
     &      STATUS,IERR )
!
!---    Read z-direction well projections
!       (duplicated across processors)  ---
!
        NVAR = LWN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,PLZ_CW,NVAR,MPI_REAL8,
     &      STATUS,IERR )
!
!---    Read x-direction well node points
!       (duplicated across processors)  ---
!
        NVAR = 2*LWN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,XP_CW,NVAR,MPI_REAL8,
     &      STATUS,IERR )
!
!---    Read y-direction well node points
!       (duplicated across processors)  ---
!
        NVAR = 2*LWN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,YP_CW,NVAR,MPI_REAL8,
     &      STATUS,IERR )
!
!---    Read z-direction well node points
!       (duplicated across processors)  ---
!
        NVAR = 2*LWN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ZP_CW,NVAR,MPI_REAL8,
     &      STATUS,IERR )
!        DO N = 1,LWN_CW
!          DO M = 1,2
!            PRINT *,'ZP_CW(',M,',',N,') = ',ZP_CW(M,N),'ID = ',ID
!          ENDDO
!        ENDDO
!
!---    Read active solutes in coupled well
!       (duplicated across processors)  ---
!
        NVAR = LSOLU_CW*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ISOLU_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!
!---    Read active species in coupled well
!       (duplicated across processors)  ---
!
        NVAR = LSPC_CW*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ISPC_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!
!---    Read active component species in coupled well
!       (duplicated across processors)  ---
!
        NVAR = LEQC*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ISOLC_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!
!---    Read active kinetic species in coupled well
!       (duplicated across processors)  ---
!
        NVAR = LEQK*LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ISOLK_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!
!---    Read active kinetic species in coupled well
!       (duplicated across processors)  ---
!
        NVAR = LN_CW
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,NSP_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
      ENDIF
!
!---  Return for no coupled wells  ---
!
      IF( N_CW.EQ.0 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Create a local copy of the arrays connecting well nodes
!     with field nodes and field nodes with well nodes
!     putting zeros for well nodes off the local
!     processor, start by looping over coupled wells
!
      DO NCW = 1,N_CW
!
!---    Loop over coupled-well nodes  ---
!
        DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          IWN_CW(NWN) = 0
!
!---      Loop over local nodes, skipping ghost cells  ---
!
          DO N = 1,NFCGC(ID+1)
            IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
            IF( ND(N).EQ.IWNG_CW(NWN) ) THEN
              IWN_CW(NWN) = N
              EXIT
            ENDIF
          ENDDO
        ENDDO
!
!---    Loop over field nodes that contain coupled-well nodes  ---
!
        DO NWF = ID_CW(5,NCW),ID_CW(6,NCW)
          IWF_CW(NWF) = 0
!
!---      Loop over local nodes, skipping ghost cells  ---
!
          DO N = 1,NFCGC(ID+1)
            IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
            IF( ND(N).EQ.IWFG_CW(NWF) ) THEN
              IWF_CW(NWF) = N
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Close the well.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_COUP_WELL_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_GRID_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Read binary input.bin file for grid data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FILES
      USE CONST

!
!----------------------Lis Modules-----------------------------------!
!
      USE LIS_STOMP

!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(3):: VARX
      INTEGER IVARX(10)
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_GRID_W'
!
!---  Allocate memory for NFCGC, number of field cells and ghost
!     cells on each processor  ---
!
      ALLOCATE( NFCGC(1:NP),STAT=ISTAT )
      CHMSG = 'NFCGC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NFC, number of field cells on each
!     processor  ---
!
      ALLOCATE( NFC(1:NP),STAT=ISTAT )
      CHMSG = 'NFC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for IDP  ---
!
      ALLOCATE( IDP(1:2,1:NP),STAT=ISTAT )
      CHMSG = 'IDP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for JDP  ---
!
      ALLOCATE( JDP(1:2,1:NP),STAT=ISTAT )
      CHMSG = 'JDP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for KDP  ---
!
      ALLOCATE( KDP(1:2,1:NP),STAT=ISTAT )
      CHMSG = 'KDP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NUKFL ---
!
      ALLOCATE( NUKFL(1:NP),STAT=ISTAT )
      CHMSG = 'NUKFL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NUKFO ---
!
      ALLOCATE( NUKFO(1:NP),STAT=ISTAT )
      CHMSG = 'NUKFO'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NUKTL ---
!
      ALLOCATE( NUKTL(1:NP),STAT=ISTAT )
      CHMSG = 'NUKTL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NUKTO ---
!
      ALLOCATE( NUKTO(1:NP),STAT=ISTAT )
      CHMSG = 'NUKTO'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Open grid1.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'grid1.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read GRAVX, GRAVY, GRAVZ (duplicated across processors)  ---
!
      NVAR = 3
      OFFSET = IOFFSET + NBYTB
      IOFFSET =  IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      GRAVX = VARX(1)
      GRAVY = VARX(2)
      GRAVZ = VARX(3)
!
!---  Read NFCGC array (duplicated across processors)  ---
!
      NVAR = NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET =  IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NFCGC,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      NFCGC_G = 0
      DO N = 1,NP
        NFCGC_G = NFCGC_G + NFCGC(N)
      ENDDO
!
!---  Allocate memory for geometry arrays (including ghost cells)  ---
!
      CALL ALLOC_GRID
!
!---  Set local starting point for local copies of nodal variables  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NFCGC(I)
      ENDDO
!
!---  Read IDP array (duplicated across processors)  ---
!
      NVAR = 2*NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IDP,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read JDP array (duplicated across processors)  ---
!
      NVAR = 2*NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,JDP,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read KDP array (duplicated across processors)  ---
!
      NVAR = 2*NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,KDP,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read local copies of XE array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*8
      OFFSET = IOFFSET + NBYTB + NC*8*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*8*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,XE,NVAR,MPI_REAL8,STATUS,IERR )
!
!---  Close the grid1.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open grid2.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'grid2.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of YE array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*8
      OFFSET = IOFFSET + NBYTB + NC*8*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*8*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,YE,NVAR,MPI_REAL8,STATUS,IERR )
!
!---  Close the grid2.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open grid3.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'grid3.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of ZE array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*8
      OFFSET = IOFFSET + NBYTB + NC*8*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*8*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ZE,NVAR,MPI_REAL8,STATUS,IERR )
!
!---  Close the grid3.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open grid4.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'grid4.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of XP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,XP,NVAR,MPI_REAL8,STATUS,IERR )
!
!---  Read local copies of YP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,YP,NVAR,MPI_REAL8,STATUS,IERR )
!
!---  Read local copies of ZP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ZP,NVAR,MPI_REAL8,STATUS,IERR )
!
!---  Read local copies of RP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,RP,NVAR,MPI_REAL8,STATUS,IERR )
!
!---  Close the grid4.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open grid5.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'grid5.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of DXGF array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,DXGF,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of DYGF array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,DYGF,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of DZGF array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,DZGF,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Close the grid5.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open grid6.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'grid6.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of DXGP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*2
      OFFSET = IOFFSET + NBYTB + NC*2*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*2*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,DXGP,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of DYGP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*2
      OFFSET = IOFFSET + NBYTB + NC*2*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*2*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,DYGP,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of DZGP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*2
      OFFSET = IOFFSET + NBYTB + NC*2*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*2*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,DZGP,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Close the grid6.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open grid7.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'grid7.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of AFX array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*2
      OFFSET = IOFFSET + NBYTB + NC*2*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*2*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,AFX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of AFY array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*2
      OFFSET = IOFFSET + NBYTB + NC*2*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*2*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,AFY,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of AFZ array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*2
      OFFSET = IOFFSET + NBYTB + NC*2*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*2*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,AFZ,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Close the grid7.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open grid8.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'grid8.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of VOL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VOL,NVAR,MPI_REAL8,STATUS,
     &  IERR )
!
!---  Read local copies of GRVX array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*2
      OFFSET = IOFFSET + NBYTB + NC*2*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*2*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,GRVX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of GRVY array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*2
      OFFSET = IOFFSET + NBYTB + NC*2*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*2*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,GRVY,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of GRVZ array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*2
      OFFSET = IOFFSET + NBYTB + NC*2*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*2*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,GRVZ,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Close the grid8.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open grid9.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'grid9.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of IXP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IXP,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read local copies of INBS array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*6
      OFFSET = IOFFSET + NBYTB + NC*6*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*6*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,INBS,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read local copies of ND array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ND,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read NUKFL (local flow equation unknowns)
!     (duplicated across processors)  ---
!
      NVAR = NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NUKFL,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read NUKFO (local flow equation offsets)
!     (duplicated across processors)  ---
!
      NVAR = NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NUKFO,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read NUKFG (global flow equation unknowns)
!     (duplicated across processors)  ---
!
      NVAR = 1
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NUKFG,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read NUKTL (local transport equation unknowns)
!     (duplicated across processors)  ---
!
      NVAR = NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NUKTL,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read NUKTO (local transport equation offsets)
!     (duplicated across processors)  ---
!
      NVAR = NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NUKTO,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read NUKTG  (global transport equation unknowns)
!     (duplicated across processors)  ---
!
      NVAR = 1
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NUKTG,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read IFLD, JFLD, KFLD, NFLD (global field limits)
!     IPFLD, JPFLD, KPFLD, NPFLD (processor field limits)
!     (duplicated across processors)  ---
!
      NVAR = 10
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IVARX,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      IFLD = IVARX(1)
      JFLD = IVARX(2)
      KFLD = IVARX(3)
      NFLD = IVARX(4)
      NXP = IVARX(5)
      ICS = IVARX(6)
      IPFLD = IVARX(7)
      JPFLD = IVARX(8)
      KPFLD = IVARX(9)
      NPFLD = IVARX(10)
      IF( NP.NE.NPFLD ) THEN
        IF( ID.EQ.0 ) PRINT *,'Processor Count Error: Number of ' //
     &    'Processors Requested ≠ Input File Processor Count'
        CALL MPI_FINALIZE( IERR )
        STOP
      ENDIF
!
!---  Allocate memory for NGHC ---
!
      ALLOCATE( NGHC(1:NP),STAT=ISTAT )
      CHMSG = 'NGHC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read NGHC array total number of ghost cells on a processor,
!     (duplicated across processors)  ---
!
      NVAR = NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NGHC,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Allocate memory for NCGC ---
!
      ALLOCATE( NCGC(1:6,1:NP),STAT=ISTAT )
      CHMSG = 'NCGC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read NCGC array, number of ghost cells to be sent in the
!     stencil directions, (duplicated across processors)  ---
!
      NVAR = 6*NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NCGC,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Allocate memory for NPGC ---
!
      ALLOCATE( NPGC(1:6,1:NP),STAT=ISTAT )
      CHMSG = 'NPGC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read NPGC array, receiving processors in the
!     stencil directions, (duplicated across processors)  ---
!
      NVAR = 6*NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NPGC,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Allocate memory for NLSGC ---
!
      ALLOCATE( NLSGC(1:NGHC(ID+1)),STAT=ISTAT )
      CHMSG = 'NLSGC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      NCX = 0
      DO I = 1,ID
        NCX = NCX + NGHC(I)
      ENDDO
      NCSX = 0
      DO I = 1,NP
        NCSX = NCSX + NGHC(I)
      ENDDO
!
!---  Read local copies of NLSGC array  ---
!
      NVAR = NGHC(ID+1)
      OFFSET = IOFFSET + NBYTB + NCX*NBYTI
      IOFFSET = IOFFSET + NCSX*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NLSGC,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Allocate memory for NLRGC ---
!
      ALLOCATE( NLRGC(1:NGHC(ID+1)),STAT=ISTAT )
      CHMSG = 'NLRGC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read local copies of NLRGC array  ---
!
      NVAR = NGHC(ID+1)
      OFFSET = IOFFSET + NBYTB + NCX*NBYTI
      IOFFSET = IOFFSET + NCSX*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NLRGC,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Close the grid9.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Load ITAG, sending/receiving tags
!     (duplicated across processors)  ---
!
      NC = 0
      DO NPRX = 1,NP
        DO NPSX = 1,NP
          NC = NC + 1
          ITAG(NPSX,NPRX) = NC
        ENDDO
      ENDDO
!
!---  Number of ghost-cell primary variables for STOMP-CO2  ---
!
      NPVX = 12
!
!---  Allocate memory for ghost-cell send buffers in the six
!     stencil directions  ---
!
      IF( NCGC(1,ID+1).GT.0 ) THEN
!        PRINT *,'SBFB(',NPVX*NCGC(1,ID+1),' ID = ',ID
        ALLOCATE( SBFB(1:NPVX*NCGC(1,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for SBFB on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(1,ID+1)
          SBFB(M) = 0.D+0
        ENDDO
      ENDIF
      IF( NCGC(2,ID+1).GT.0 ) THEN
!        PRINT *,'SBFS(',NPVX*NCGC(2,ID+1),' ID = ',ID
        ALLOCATE( SBFS(1:NPVX*NCGC(2,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for SBFS on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(2,ID+1)
          SBFS(M) = 0.D+0
        ENDDO
      ENDIF
      IF( NCGC(3,ID+1).GT.0 ) THEN
!        PRINT *,'SBFW(',NPVX*NCGC(3,ID+1),' ID = ',ID
        ALLOCATE( SBFW(1:NPVX*NCGC(3,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for SBFW on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(3,ID+1)
          SBFW(M) = 0.D+0
        ENDDO
      ENDIF
      IF( NCGC(4,ID+1).GT.0 ) THEN
!        PRINT *,'SBFE(',NPVX*NCGC(4,ID+1),' ID = ',ID
        ALLOCATE( SBFE(1:NPVX*NCGC(4,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for SBFE on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(4,ID+1)
          SBFE(M) = 0.D+0
        ENDDO
      ENDIF
      IF( NCGC(5,ID+1).GT.0 ) THEN
!        PRINT *,'SBFN(',NPVX*NCGC(5,ID+1),' ID = ',ID
        ALLOCATE( SBFN(1:NPVX*NCGC(5,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for SBFN on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(5,ID+1)
          SBFN(M) = 0.D+0
        ENDDO
      ENDIF
      IF( NCGC(6,ID+1).GT.0 ) THEN
!        PRINT *,'SBFT(',NPVX*NCGC(6,ID+1),' ID = ',ID
        ALLOCATE( SBFT(1:NPVX*NCGC(6,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for SBFT on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(6,ID+1)
          SBFT(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Allocate memory for ghost-cell receive buffers in the six
!     stencil directions  ---
!
      IF( NCGC(6,ID+1).GT.0 ) THEN
!        PRINT *,'RBFB(',NPVX*NCGC(6,ID+1),' ID = ',ID
        ALLOCATE( RBFB(1:NPVX*NCGC(6,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for RBFB on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(6,ID+1)
          RBFB(M) = 0.D+0
        ENDDO
      ENDIF
      IF( NCGC(5,ID+1).GT.0 ) THEN
!        PRINT *,'RBFS(',NPVX*NCGC(5,ID+1),' ID = ',ID
        ALLOCATE( RBFS(1:NPVX*NCGC(5,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for RBFS on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(5,ID+1)
          RBFS(M) = 0.D+0
        ENDDO
      ENDIF
      IF( NCGC(4,ID+1).GT.0 ) THEN
!        PRINT *,'RBFW(',NPVX*NCGC(4,ID+1),' ID = ',ID
        ALLOCATE( RBFW(1:NPVX*NCGC(4,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for RBFW on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(4,ID+1)
          RBFW(M) = 0.D+0
        ENDDO
      ENDIF
      IF( NCGC(3,ID+1).GT.0 ) THEN
!        PRINT *,'RBFE(',NPVX*NCGC(3,ID+1),' ID = ',ID
        ALLOCATE( RBFE(1:NPVX*NCGC(3,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for RBFE on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(3,ID+1)
          RBFE(M) = 0.D+0
        ENDDO
      ENDIF
      IF( NCGC(2,ID+1).GT.0 ) THEN
!        PRINT *,'RBFN(',NPVX*NCGC(2,ID+1),' ID = ',ID
        ALLOCATE( RBFN(1:NPVX*NCGC(2,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for RBFN on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(2,ID+1)
          RBFN(M) = 0.D+0
        ENDDO
      ENDIF
      IF( NCGC(1,ID+1).GT.0 ) THEN
!        PRINT *,'RBFT(',NPVX*NCGC(1,ID+1),' ID = ',ID
        ALLOCATE( RBFT(1:NPVX*NCGC(1,ID+1)),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          PRINT *,'Allocation Error for RBFT on ID = ',ID
        ENDIF
        DO M = 1,NPVX*NCGC(1,ID+1)
          RBFT(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Initialize ghost cell map  ---
!
      DO L = 1,NFCGC(ID+1)
        IGHC(L) = 0
      ENDDO
!
!---  Create a local node connection map and identify ghost cells  ---
!
      KFLD_L = KDP(2,ID+1)-KDP(1,ID+1)+1
      JFLD_L = JDP(2,ID+1)-JDP(1,ID+1)+1
      IFLD_L = IDP(2,ID+1)-IDP(1,ID+1)+1
      DO K = KDP(1,ID+1),KDP(2,ID+1)
        DO J = JDP(1,ID+1),JDP(2,ID+1)
          DO I = IDP(1,ID+1),IDP(2,ID+1)
            KX = K - KDP(1,ID+1) + 1
            JX = J - JDP(1,ID+1) + 1
            IX = I - IDP(1,ID+1) + 1
            N = (KX-1)*IFLD_L*JFLD_L + (JX-1)*IFLD_L + IX
!
!---        Bottom connection  ---
!
            ICM(1,N) = 0
            IF( KX.GT.1 ) THEN
              NB = N-IFLD_L*JFLD_L
              IF( IXP(NB).GT.0 .AND. INBS(1,N).EQ.0 ) ICM(1,N) = NB
            ENDIF
!
!---        Bottom ghost cells  ---
!
            IF( KX.EQ.1 .AND. K.GT.1 ) IGHC(N) = 1
!
!---        South connection  ---
!
            ICM(2,N) = 0
            IF( JX.GT.1 ) THEN
              NS = N-IFLD_L
              IF( IXP(NS).GT.0 .AND. INBS(2,N).EQ.0 ) ICM(2,N) = NS
            ENDIF
!
!---        South ghost cells  ---
!
            IF( JX.EQ.1 .AND. J.GT.1 ) IGHC(N) = 1
!
!---        West connection  ---
!
            ICM(3,N) = 0
            IF( IX.GT.1 ) THEN
              NW = N-1
              IF( IXP(NW).GT.0 .AND. INBS(3,N).EQ.0 ) ICM(3,N) = NW
            ENDIF
!
!---        West ghost cells  ---
!
            IF( IX.EQ.1 .AND. I.GT.1 ) IGHC(N) = 1
!
!---        East connection  ---
!
            ICM(4,N) = 0
            IF( IX.LT.IFLD_L ) THEN
              NE = N+1
              IF( IXP(NE).GT.0 .AND. INBS(4,N).EQ.0 ) ICM(4,N) = NE
            ENDIF
!
!---        East ghost cells  ---
!
            IF( IX.EQ.IFLD_L .AND. I.LT.IDP(2,NP) ) IGHC(N) = 1
!
!---        North connection  ---
!
            ICM(5,N) = 0
            IF( JX.LT.JFLD_L ) THEN
              NN = N+IFLD_L
              IF( IXP(NN).GT.0 .AND. INBS(5,N).EQ.0 ) ICM(5,N) = NN
            ENDIF
!
!---        North ghost cells  ---
!
            IF( JX.EQ.JFLD_L .AND. J.LT.JDP(2,NP) ) IGHC(N) = 1
!
!---        Top connection  ---
!
            ICM(6,N) = 0
            IF( KX.LT.KFLD_L ) THEN
              NT = N+IFLD_L*JFLD_L
              IF( IXP(NT).GT.0 .AND. INBS(6,N).EQ.0 ) ICM(6,N) = NT
            ENDIF
!
!---        Top ghost cells  ---
!
            IF( KX.EQ.KFLD_L .AND. K.LT.KDP(2,NP) ) IGHC(N) = 1
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_GRID_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_PROP_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Read binary input.bin file for property data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE TABL
      USE SOLTN
      USE PROP
      USE MPI
      USE GRID
      USE GLB_PAR
      USE GEO_MECH
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
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_PROP_W'
!
!---  Allocate local memory for property arrays ---
!
      CALL ALLOC_PROP
!
!---  Allocate local memory for lookup table arrays ---
!
      CALL ALLOC_TABL
!
!---  Set local starting point for local copies of nodal variables  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NFCGC(I)
      ENDDO
!
!---  Open mech.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'mech.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of IZ array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IZ,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!      PRINT *,'IZ(',ND(2),') = ',IZ(2),
!     &  'IZ(',ND(NVAR),') = ',IZ(NVAR),'ID = ',ID
!
!---  Read local copies of RHOS array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NC*NBYTR + NBYTB
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,RHOS,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'RHOS(',ND(2),') = ',RHOS(2),
!     &  'RHOS(',ND(NVAR),') = ',RHOS(NVAR),'ID = ',ID
!
!---  Read local copies of CMP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*4
      OFFSET = IOFFSET + NBYTB + NC*4*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*4*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,CMP,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'CMP(1,',ND(2),') = ',CMP(1,2),
!     &  'CMP(1,',ND(NVAR),') = ',CMP(1,NVAR),'ID = ',ID
!
!---  Read local copies of POR array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*6
      OFFSET = IOFFSET + NBYTB + NC*6*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*6*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,POR,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!
!---  Read local copies of TOR array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*6
      OFFSET = IOFFSET + NBYTB + NC*6*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*6*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TOR,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!
!---  Read local copies of ITOR array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ITOR,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Close the mech.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open hydr.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'hydr.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of PERM array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*9
      OFFSET = IOFFSET + NBYTB + NC*9*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*9*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,PERM,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!
!---  Read local copies of IPRF array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IPRF,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Close the hydr.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open satu1.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'satu1.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of SCHR array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*LSCHR
      OFFSET = IOFFSET + NBYTB + NC*LSCHR*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*LSCHR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,SCHR,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!
!---  Close the satu1.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open satu2.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'satu2.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of ISCHR array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISCHR,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Read local copies of ISM array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISM,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Read local copies of ISLTBL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*4
      OFFSET = IOFFSET + NBYTB + NC*4*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*4*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISLTBL,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Read local copies of IRLTBL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*4
      OFFSET = IOFFSET + NBYTB + NC*4*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*4*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IRLTBL,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Close the satu2.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open rel_perm.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'perm.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of ISKP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISKP,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Read local copies of RPLC array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*LRPLC
      OFFSET = IOFFSET + NBYTB + NC*LRPLC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*LRPLC*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,RPLC,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!
!---  Read local copies of IRPL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IRPL,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Close the rel_perm.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open comp.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'comp.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of PCMP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,PCMP,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'PCMP(',ND(2),') = ',PCMP(2),
!     &  'PCMP(',ND(NVAR),') = ',PCMP(NVAR),'ID = ',ID
!
!---  Read local copies of TCMP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TCMP,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'TCMP(',ND(2),') = ',TCMP(2),
!     &  'TCMP(',ND(NVAR),') = ',TCMP(NVAR),'ID = ',ID
!
!---  Close the comp.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open tabl.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'tabl.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of NTBL  ---
!
      NVAR = 1
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NTBL,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Read local copies of TBLX array  ---
!
      NVAR = LTBL
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TBLX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!
!---  Read local copies of TBLY array  ---
!
      NVAR = LTBL
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TBLY,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!
!---  Read local copies of TBLDDX array  ---
!
      NVAR = LTBL
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TBLDDX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!
!---  Read local copies of TBLDDY array  ---
!
      NVAR = LTBL
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TBLDDY,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!
!---  Close the table.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_PROP_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_REACT_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Read binary input.bin file for ECKEChem (reactive transport) data
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PROP
      USE HYST
      USE GRID
      USE GLB_PAR
      USE FILES
      USE FDVP
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
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARX,VARZ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVARX,IVARZ
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_REACT_W'
!
!---  Allocate and initialize memory, global and local, for 
!     reactive species arrays  ---
!
      CALL ALLOC_REACT
      CALL INTLZ_REACT
!
!---  Set local starting point for local copies of nodal variables  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NFCGC(I)
      ENDDO
!
!---  Allocate local temporary real and integer nodal arrays 
!     (including ghost cells)  ---
!
      ALLOCATE( VARX(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'VARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IVARX(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read number of reactive species integers
!     (duplicated across processors)  ---
!
      NVAR = 19
      ALLOCATE( IVARZ(1:NVAR),STAT=ISTAT )
      CHMSG = 'IVARZ'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IVARZ,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
      IACTEX = IVARZ(1)
      IACTV = IVARZ(2)
      ISP_IEDL = IVARZ(3)
      NEQE = IVARZ(4)
      NEQC = IVARZ(5)
      NEQK = IVARZ(6)
      NESITE = IVARZ(7)
      NRCE = IVARZ(8)
      NRCK = IVARZ(9)
      NRTSI = IVARZ(10)
      NSPC = IVARZ(11)
      NSPE = IVARZ(12)
      NSPG = IVARZ(13)
      NSPK = IVARZ(14)
      NSPL = IVARZ(15)
      NSPLK = IVARZ(16)
      NSPN = IVARZ(17)
      NSPR = IVARZ(18)
      NSPS = IVARZ(19)
!      PRINT *,'IVARZ = ',(IVARZ(M),M=1,NVAR),' ID = ',ID
      DEALLOCATE( IVARZ,STAT=ISTAT )
      CHMSG = 'IVARZ'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read reactive transport parameters, indices, and variables
!     if reactive species are modeled  ---
!
      IF( NEQE+NEQC+NEQK.GT.0 ) THEN
!
!---    Read index for 1 equations
!       (duplicated across processors)  ---
!
        NVAR = (LSEE+2)*LEQE
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IEQ_E,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'IEQ_E = ',((IEQ_E(L,M),L=1,LSEE+2),M=1,LEQE),
!     &    ' ID = ',ID
!
!---    Read index for conservation equations
!       (duplicated across processors)  ---
!
        NVAR = (LSEC+1)*LEQC
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IEQ_C,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'IEQ_C = ',((IEQ_C(L,M),L=1,LSEC+1),M=1,LEQC),
!     &    ' ID = ',ID
!
!---    Read index for kinetic equations
!       (duplicated across processors)  ---
!
        NVAR = (LSEK+LREK+2)*LEQK
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IEQ_K,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'IEQ_K = ',((IEQ_K(L,M),L=1,LSEK+LREK+2),M=1,LEQK),
!     &    ' ID = ',ID
!
!---    Read index for reactive equation sequencing
!       (duplicated across processors)  ---
!
        NVAR = LSPR
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IEQ_S,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'IEQ_S = ',(IEQ_S(L),L=1,LSPR),
!     &    ' ID = ',ID
!
!---    Read index for reactive species
!       (duplicated across processors)  ---
!
        NVAR = LSPE
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ISP_E,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'ISP_E = ',(ISP_E(L),L=1,LSPE),
!     &    ' ID = ',ID
!
!---    Read index for reactive species
!       (duplicated across processors)  ---
!
        NVAR = LSPE
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IEL_LK,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'IEL_LK = ',(IEL_LK(L),L=1,LSPE),
!     &    ' ID = ',ID
!
!---    Read index for reactive species
!       (duplicated across processors)  ---
!
        NVAR = LSPR
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ISP_MN,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'ISP_MN = ',(ISP_MN(L),L=1,LSPR),
!     &    ' ID = ',ID
!
!---    Read index for reactive species
!       (duplicated across processors)  ---
!
        NVAR = LEQE+LEQC+LEQK
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ISP_S,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'ISP_S = ',(ISP_S(L),L=1,LEQE+LEQC+LEQK),
!     &    ' ID = ',ID
!
!---    Read index for reactive species
!       (duplicated across processors)  ---
!
        NVAR = LSOLU+LSPT
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IMMB,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'IMMB = ',(IMMB(L),L=1,LSOLU+LSPT),
!     &    ' ID = ',ID
!
!---    Read index for reactive species
!       (duplicated across processors)  ---
!
        NVAR = (LSPK+3)*LRCK
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IRC_K,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'IRC_K = ',((IRC_K(L,M),L=1,LSPK+3),M=1,LRCK),
!     &    ' ID = ',ID
!
!---    Read index for reactive species
!       (duplicated across processors)  ---
!
        NVAR = LSPK+11
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IRCKN,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'IRCKN = ',(IRCKN(L),L=1,LSPK+11),
!     &    ' ID = ',ID
!
!---    Read index for reactive species
!       (duplicated across processors)  ---
!
        NVAR = LRCK
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IRCKT,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'IRCKT = ',(IRCKT(L),L=1,LRCK),
!     &    ' ID = ',ID
!
!---    Read index for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 2*LNGC+LSPR+14
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ISPLK,NVAR,MPI_INTEGER,
     &    STATUS,IERR)
!        PRINT *,'ISPLK = ',(ISPLK(L),L=1,2*LNGC+LSPR+14),
!     &    ' ID = ',ID
!
!---    Read index for reactive species  ---
!
        DO J = 1,LSPS
          NVAR = NFCGC(ID+1)
          OFFSET = IOFFSET + NC*NBYTI + NBYTB
          IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
          CALL MPI_FILE_READ_AT( IRD,OFFSET,IVARX,NVAR,MPI_INTEGER,
     &      STATUS,IERR)
          DO N = 1,NFCGC(ID+1)
            ISP_OW(J,N) = IVARX(N)
          ENDDO
!          PRINT *,'ISP_OW(',J,',1) = ',ISP_OW(J,1),
!     &      'ISP_OW(',J,',',ND(NVAR),') = ',ISP_OW(J,NVAR),' ID = ',ID
        ENDDO
!
!---    Read index for reactive species  ---
!
        DO M = 1,LSPR
          NVAR = NFCGC(ID+1)
          OFFSET = IOFFSET + NC*NBYTI + NBYTB
          IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
          CALL MPI_FILE_READ_AT( IRD,OFFSET,IVARX,NVAR,MPI_INTEGER,
     &      STATUS,IERR)
          DO N = 1,NFCGC(ID+1)
            IC_SP(N,M) = IVARX(N)
          ENDDO
!          PRINT *,'IC_SP(1,',M,') = ',IC_SP(1,M),
!     &      'IC_SP(',ND(NVAR),',',M,') = ',IC_SP(NVAR,M),' ID = ',ID
        ENDDO
!
!---    Read parameters for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 4
        ALLOCATE( VARZ(1:NVAR),STAT=ISTAT )
        CHMSG = 'VARZ'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,VARZ,NVAR,MPI_REAL8,
     &    STATUS,IERR)
        ACTVC = VARZ(1)
        CMIN = VARZ(2)
        SP_MDG = VARZ(3)
        SP_MDL = VARZ(4)
        DEALLOCATE( VARZ,STAT=ISTAT )
        CHMSG = 'VARZ'
        CALL DEALLOC_ERROR( CHMSG,ISTAT )
!        PRINT *,'ACTVC = ',ACTVC,' CMIN = ',CMIN,
!     &    ' SP_MDG = ',SP_MDG,' SP_MDL = ',SP_MDL,' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 3
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SP_SDCL,NVAR,MPI_REAL8,
     &    STATUS,IERR)
!        PRINT *,'SP_SDCL = ',(SP_SDCL(L),L=1,3),
!     &    ' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = LMC
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,CFMX,NVAR,MPI_REAL8,
     &    STATUS,IERR)
!        PRINT *,'LMC = ',LMC,' CFMX = ',(CFMX(L),L=1,LMC),
!     &    ' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = LSPR
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,CHARG,NVAR,MPI_REAL8,
     &    STATUS,IERR)
!        PRINT *,'LSPR = ',LSPR,' CHARG = ',(CHARG(L),L=1,LSPR),
!     &    ' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = LSEC*LEQC
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,EQ_C,NVAR,MPI_REAL8,
     &    STATUS,IERR)
!        PRINT *,'EQ_C = ',((EQ_C(L,M),L=1,LSEC),M=1,LEQC),
!     &    ' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = LSEE*LEQE
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,EQ_E,NVAR,MPI_REAL8,
     &    STATUS,IERR)
!        PRINT *,'EQ_E = ',((EQ_E(L,M),L=1,LSEE),M=1,LEQE),
!     &    ' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = (LSEK+LREK)*LEQK
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,EQ_K,NVAR,MPI_REAL8,
     &    STATUS,IERR)
!        PRINT *,'EQ_K = ',((EQ_K(L,M),L=1,LSEK+LREK),M=1,LEQK),
!     &    ' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 5*LRCE
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,RC_E,NVAR,MPI_REAL8,
     &    STATUS,IERR)
!        PRINT *,'RC_E = ',((RC_E(L,M),L=1,5),M=1,LRCE),
!     &    ' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = (LSPK+11)*LCKN*LRCK
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,RC_K,NVAR,MPI_REAL8,
     &    STATUS,IERR)
!        PRINT *,'RC_K = ',(RC_K(K,1,1),K=1,LSPK+11),' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 3*LSPL
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SP_L,NVAR,MPI_REAL8,
     &    STATUS,IERR)
!        PRINT *,'SP_L = ',((SP_L(L,M),L=1,3),M=1,LSPL),
!     &    ' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 2*LSPS
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SP_S,NVAR,MPI_REAL8,
     &    STATUS,IERR)
!        PRINT *,'SP_S = ',((SP_S(L,M),L=1,2),M=1,LSPS),
!     &    ' ID = ',ID
!
!---    Read parameter for reactive species  ---
!
        DO I = 1,3
          DO J = 1,LSPS
            NVAR = NFCGC(ID+1)
            OFFSET = IOFFSET + NC*NBYTR + NBYTB
            IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
            CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &        STATUS,IERR)
            DO N = 1,NFCGC(ID+1)
              RS_S(I,J,N) = VARX(N)
            ENDDO
          ENDDO
        ENDDO
!        PRINT *,'RS_S = ',((RS_S(I,J,1),I=1,3),J=1,LSPS),' ID = ',ID
!
!---    Read parameter for reactive species  ---
!
        DO M = 1,LEQC+LEQK
          NVAR = NFCGC(ID+1)
          OFFSET = IOFFSET + NC*NBYTR + NBYTB
          IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
          CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &      STATUS,IERR)
          DO N = 1,NFCGC(ID+1)
            YSPG(N,M) = VARX(N)
          ENDDO
        ENDDO
!        PRINT *,'YSPG = ',(YSPG(1,M),M=1,LEQC+LEQK),' ID = ',ID
!
!---    Read parameter for reactive species  ---
!
        DO M = 1,LEQC+LEQK
          NVAR = NFCGC(ID+1)
          OFFSET = IOFFSET + NC*NBYTR + NBYTB
          IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
          CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &      STATUS,IERR)
          DO N = 1,NFCGC(ID+1)
            YSPL(N,M) = VARX(N)
          ENDDO
        ENDDO
!        PRINT *,'YSPL = ',(YSPL(1,M),M=1,LEQC+LEQK),' ID = ',ID
!
!---    Read parameter for reactive species  ---
!
        DO M = 1,LSPR
          NVAR = NFCGC(ID+1)
          OFFSET = IOFFSET + NC*NBYTR + NBYTB
          IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
          CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &      STATUS,IERR)
          DO N = 1,NFCGC(ID+1)
            SP_C(N,M) = VARX(N)
          ENDDO
        ENDDO
!        PRINT *,'SP_C = ',(SP_C(1,M),M=1,LSPR),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LRCE
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,RCNME,NVAR,MPI_CHAR,
     &    STATUS,IERR)
!        PRINT *,'RCNME = ',(RCNME(L),L=1,LRCE),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LRCK
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,RCNMK,NVAR,MPI_CHAR,
     &    STATUS,IERR)
!        PRINT *,'RCNMK = ',(RCNMK(L),L=1,LRCK),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LEQC
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNMC,NVAR,MPI_CHAR,
     &    STATUS,IERR)
!        PRINT *,'SPNMC = ',(SPNMC(L),L=1,LEQC),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LEQK
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNMK,NVAR,MPI_CHAR,
     &    STATUS,IERR)
!        PRINT *,'SPNMK = ',(SPNMK(L),L=1,LEQK),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LSPG
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNMG,NVAR,MPI_CHAR,
     &    STATUS,IERR)
!        PRINT *,'SPNMG = ',(SPNMG(L),L=1,LSPG),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LSPL
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNML,NVAR,MPI_CHAR,
     &    STATUS,IERR)
!        PRINT *,'SPNML = ',(SPNML(L),L=1,LSPL),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LSPS
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNMS,NVAR,MPI_CHAR,
     &    STATUS,IERR)
!        PRINT *,'SPNMS = ',(SPNMS(L),L=1,LSPS),' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LSPE
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNME,NVAR,MPI_CHAR,
     &    STATUS,IERR)
!        PRINT *,'SPNME = ',(SPNME(L),L=1,LSPE),' ID = ',ID
      ENDIF
!
!---  Deallocate local temporary real and integer nodal arrays 
!     (including ghost cells)  ---
!
      DEALLOCATE( VARX,STAT=ISTAT )
      CHMSG = 'VARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
      DEALLOCATE( IVARX,STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_REACT_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_SOLU_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Read binary input.bin file for solution control and output data
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE OUTPU
      USE GRID
      USE GLB_PAR
      USE FDVP
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
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
      CHARACTER(64), DIMENSION(2) :: FNSFX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_SOLU_W'
!
!---  Open solu.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'solu.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Initialize rror messaging variables  ---
!
      DO M = 1,4
        I_ERR(M) = 0
        M_ERR(M) = ''
      ENDDO
      I_ERR(4) = NP + 1
      R_ERR = 0.D+0
!
!---  Initialize cummulative offset  ---
!
      IOFFSET = 0
!
!---  Read parameter variables (duplicated across processors)  ---
!
      NVAR = 341
      ALLOCATE( IVARX(1:NVAR),STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IVARX,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      L_BH = IVARX(1)
      L_CW = IVARX(2)
      L_DP = IVARX(3)
      L_EC = IVARX(4)
      L_FRC = IVARX(5)
      L_LV = IVARX(6)
      L_LW = IVARX(7)
      L_SFC = IVARX(8)
      LAD = IVARX(9)
      LALC = IVARX(10)
      LAN = IVARX(11)
      LANI = IVARX(12)
      LANW = IVARX(13)
      LATM = IVARX(14)
      LBAL = IVARX(15)
      LBC = IVARX(16)
      LBC_BH = IVARX(17)
      LBC_EC = IVARX(18)
      LBC_FRC = IVARX(19)
      LBC_GM = IVARX(20)
      LBCA = IVARX(21)
      LBCC = IVARX(22)
      LBCG = IVARX(23)
      LBCGC = IVARX(24)
      LBCH = IVARX(25)
      LBCI = IVARX(26)
      LBCIN = IVARX(27)
      LBCIN_GM = IVARX(28)
      LBCL = IVARX(29)
      LBCN = IVARX(30)
      LBCN2 = IVARX(31)
      LBCS = IVARX(32)
      LBCT = IVARX(33)
      LBCU = IVARX(34)
      LBCV = IVARX(35)
      LBCV_GM = IVARX(36)
      LBD = IVARX(37)
      LBN_BH = IVARX(38)
      LBN_BHC = IVARX(39)
      LBR = IVARX(40)
      LBTM = IVARX(41)
      LBTM_GM = IVARX(42)
      LC = IVARX(43)
      LCAT = IVARX(44)
      LCDC = IVARX(45)
      LCDP = IVARX(46)
      LCDS = IVARX(47)
      LCH_HT = IVARX(48)
      LCHEM = IVARX(49)
      LCKN = IVARX(50)
      LCMP = IVARX(51)
      LCN = IVARX(52)
      LCN_HT = IVARX(53)
      LCOAX_BH = IVARX(54)
      LCP_HT = IVARX(55)
      LD = IVARX(56)
      LDCO2 = IVARX(57)
      LEPD = IVARX(58)
      LEQC = IVARX(59)
      LEQE = IVARX(60)
      LEQK = IVARX(61)
      LESITE = IVARX(62)
      LF_FRC = IVARX(63)
      LF_FRCC = IVARX(64)
      LFC_BH = IVARX(65)
      LFC_FRC = IVARX(66)
      LFD = IVARX(67)
      LFD_DP = IVARX(68)
      LFD_EC = IVARX(69)
      LFDA = IVARX(70)
      LFDC = IVARX(71)
      LFDCR = IVARX(72)
      LFDD = IVARX(73)
      LFDG = IVARX(74)
      LFDGC = IVARX(75)
      LFDH = IVARX(76)
      LFDI = IVARX(77)
      LFDL = IVARX(78)
      LFDM = IVARX(79)
      LFDN = IVARX(80)
      LFDN2 = IVARX(81)
      LFDNH = IVARX(82)
      LFDR = IVARX(83)
      LFDRG = IVARX(84)
      LFDRL = IVARX(85)
      LFDRN = IVARX(86)
      LFDS = IVARX(87)
      LFDT = IVARX(88)
      LFEN = IVARX(89)
      LFILES = IVARX(90)
      LFW = IVARX(91)
      LFX = IVARX(92)
      LFX_MPI = IVARX(93)
      LFXY = IVARX(94)
      LFY = IVARX(95)
      LFY_MPI = IVARX(96)
      LFYZ = IVARX(97)
      LFZ = IVARX(98)
      LFZ_MPI = IVARX(99)
      LFZX = IVARX(100)
      LG = IVARX(101)
      LGC = IVARX(102)
      LHBW = IVARX(103)
      LHE_HT = IVARX(104)
      LHF_HT = IVARX(105)
      LHYD = IVARX(106)
      LI_BH = IVARX(107)
      LINC = IVARX(108)
      LINH = IVARX(109)
      LIS = IVARX(110)
      LJA = IVARX(111)
      LJB = IVARX(112)
      LJC = IVARX(113)
      LJC_BH = IVARX(114)
      LJC_GM = IVARX(115)
      LJD = IVARX(116)
      LJE = IVARX(117)
      LJF = IVARX(118)
      LJG = IVARX(119)
      LJG_BCF = IVARX(120)
      LJG_BCM = IVARX(121)
      LJG_BH = IVARX(122)
      LJG_FCB = IVARX(123)
      LJG_FCM = IVARX(124)
      LJG_FRC = IVARX(125)
      LJG_GM = IVARX(126)
      LJG_MCB = IVARX(127)
      LJG_MCF = IVARX(128)
      LJH = IVARX(129)
      LJH_BCF = IVARX(130)
      LJH_BCM = IVARX(131)
      LJH_BH = IVARX(132)
      LJH_FCB = IVARX(133)
      LJH_FCM = IVARX(134)
      LJH_FRC = IVARX(135)
      LJH_GM = IVARX(136)
      LJH_MCB = IVARX(137)
      LJH_MCF = IVARX(138)
      LJI = IVARX(139)
      LJJ = IVARX(140)
      LJK = IVARX(141)
      LJK_BCF = IVARX(142)
      LJK_BCM = IVARX(143)
      LJK_BH = IVARX(144)
      LJK_FCB = IVARX(145)
      LJK_FCM = IVARX(146)
      LJK_FRC = IVARX(147)
      LJK_MCB = IVARX(148)
      LJK_MCF = IVARX(149)
      LJL = IVARX(150)
      LJL_BCF = IVARX(151)
      LJL_BCM = IVARX(152)
      LJL_BH = IVARX(153)
      LJL_FCB = IVARX(154)
      LJL_FCM = IVARX(155)
      LJL_FRC = IVARX(156)
      LJL_MCB = IVARX(157)
      LJL_MCF = IVARX(158)
      LJM = IVARX(159)
      LJN = IVARX(160)
      LJN_BH = IVARX(161)
      LJO = IVARX(162)
      LJO_GM = IVARX(163)
      LL = IVARX(164)
      LM = IVARX(165)
      LMC = IVARX(166)
      LMCG = IVARX(167)
      LMNP = IVARX(168)
      LMPH = IVARX(169)
      LN = IVARX(170)
      LN_BH = IVARX(171)
      LN_BHC = IVARX(172)
      LN_CW = IVARX(173)
      LN_LW = IVARX(174)
      LN2 = IVARX(175)
      LNAF = IVARX(176)
      LNC_FRC = IVARX(177)
      LNCF = IVARX(178)
      LNEU = IVARX(179)
      LNGC = IVARX(180)
      LNHC = IVARX(181)
      LNNF = IVARX(182)
      LNNGC = IVARX(183)
      LNOTES = IVARX(184)
      LNW = IVARX(185)
      LNWN = IVARX(186)
      LNWS = IVARX(187)
      LNWT = IVARX(188)
      LNWV = IVARX(189)
      LO_PH = IVARX(190)
      LO_TH = IVARX(191)
      LOBDS = IVARX(192)
      LOBDT = IVARX(193)
      LOUPV = IVARX(194)
      LP_MPI = IVARX(195)
      LP_TA = IVARX(196)
      LPC = IVARX(197)
      LPCF = IVARX(198)
      LPE_HT = IVARX(199)
      LPF_EOR = IVARX(200)
      LPH = IVARX(201)
      LPLANT = IVARX(202)
      LPOLYC = IVARX(203)
      LPOLYN = IVARX(204)
      LPP_HT = IVARX(205)
      LPT = IVARX(206)
      LPTA = IVARX(207)
      LPTM = IVARX(208)
      LPX_MPI = IVARX(209)
      LPY_MPI = IVARX(210)
      LPZ_MPI = IVARX(211)
      LR = IVARX(212)
      LRC = IVARX(213)
      LRCE = IVARX(214)
      LRCG = IVARX(215)
      LRCK = IVARX(216)
      LRCL = IVARX(217)
      LRCN = IVARX(218)
      LRCS = IVARX(219)
      LRCT = IVARX(220)
      LREF = IVARX(221)
      LREK = IVARX(222)
      LREL = IVARX(223)
      LREM = IVARX(224)
      LRFN = IVARX(225)
      LRK = IVARX(226)
      LRPGC = IVARX(227)
      LRPL = IVARX(228)
      LRPLC = IVARX(229)
      LRPNC = IVARX(230)
      LS = IVARX(231)
      LSALC = IVARX(232)
      LSCHR = IVARX(233)
      LSEC = IVARX(234)
      LSEE = IVARX(235)
      LSEK = IVARX(236)
      LSF = IVARX(237)
      LSFCA = IVARX(238)
      LSFCC = IVARX(239)
      LSFCN = IVARX(240)
      LSFCP = IVARX(241)
      LSFCT = IVARX(242)
      LSFDOM = IVARX(243)
      LSFV = IVARX(244)
      LSFVGC = IVARX(245)
      LSOLSR = IVARX(246)
      LSOLU = IVARX(247)
      LSOLU_BH = IVARX(248)
      LSOLU_CW = IVARX(249)
      LSP = IVARX(250)
      LSPBC = IVARX(251)
      LSPC_CW = IVARX(252)
      LSPE = IVARX(253)
      LSPG = IVARX(254)
      LSPILL = IVARX(255)
      LSPK = IVARX(256)
      LSPL = IVARX(257)
      LSPLK = IVARX(258)
      LSPN = IVARX(259)
      LSPR = IVARX(260)
      LSPS = IVARX(261)
      LSPT = IVARX(262)
      LSR = IVARX(263)
      LSR_BH = IVARX(264)
      LSR_FRC = IVARX(265)
      LSRX = IVARX(266)
      LSRY = IVARX(267)
      LSRZ = IVARX(268)
      LSTC = IVARX(269)
      LSTM = IVARX(270)
      LSTM_BH = IVARX(271)
      LSTM_FRC = IVARX(272)
      LSU = IVARX(273)
      LSV = IVARX(274)
      LSW = IVARX(275)
      LSX = IVARX(276)
      LSXC = IVARX(277)
      LSXG = IVARX(278)
      LSXGC = IVARX(279)
      LSXL = IVARX(280)
      LSXLC = IVARX(281)
      LSXN = IVARX(282)
      LSXN2 = IVARX(283)
      LSXNC = IVARX(284)
      LSXS = IVARX(285)
      LSXT = IVARX(286)
      LSY = IVARX(287)
      LSYC = IVARX(288)
      LSYG = IVARX(289)
      LSYGC = IVARX(290)
      LSYL = IVARX(291)
      LSYLC = IVARX(292)
      LSYN = IVARX(293)
      LSYN2 = IVARX(294)
      LSYNC = IVARX(295)
      LSYS = IVARX(296)
      LSYT = IVARX(297)
      LSZ = IVARX(298)
      LSZC = IVARX(299)
      LSZG = IVARX(300)
      LSZGC = IVARX(301)
      LSZL = IVARX(302)
      LSZLC = IVARX(303)
      LSZN = IVARX(304)
      LSZN2 = IVARX(305)
      LSZNC = IVARX(306)
      LSZS = IVARX(307)
      LSZT = IVARX(308)
      LSZW = IVARX(309)
      LT = IVARX(310)
      LT_BH = IVARX(311)
      LT_FRC = IVARX(312)
      LT_FRCC = IVARX(313)
      LT_PH = IVARX(314)
      LT_TA = IVARX(315)
      LT_TH = IVARX(316)
      LTBL = IVARX(317)
      LTC_FRC = IVARX(318)
      LTP_HT = IVARX(319)
      LUGR = IVARX(320)
      LUK = IVARX(321)
      LUK_BH = IVARX(322)
      LUK_CW = IVARX(323)
      LUK_SFC = IVARX(324)
      LUKW = IVARX(325)
      LVIC_FRC = IVARX(326)
      LVPLOT = IVARX(327)
      LVREF = IVARX(328)
      LWELL = IVARX(329)
      LWF_CW = IVARX(330)
      LWF_LW = IVARX(331)
      LWI_CW = IVARX(332)
      LWI_LW = IVARX(333)
      LWN_CW = IVARX(334)
      LWN_LW = IVARX(335)
      LWSI = IVARX(336)
      LWT_CW = IVARX(337)
      LWTI = IVARX(338)
      LWTP_CW = IVARX(339)
      LXP_FRC = IVARX(340)
      LXYZG = IVARX(341)
      DEALLOCATE( IVARX,STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read time variables (duplicated across processors)  ---
!
      NVAR = 18
      ALLOCATE( VARX(1:NVAR),STAT=ISTAT )
      CHMSG = 'VARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      TM = VARX(1)
      TMMX = VARX(2)
      TMPR = VARX(3)
      DT = VARX(4)
      DTI = VARX(5)
      DTMX = VARX(6)
      DTMN = VARX(7)
      DTAF = VARX(8)
      DTCF = VARX(9)
      DTO = VARX(10)
      DTSO = VARX(11)
      RSDMX = VARX(12)
      RLXF = VARX(13)
      CRNTMXC = VARX(14)
      RTOL_PETSC = VARX(15)
      ATOL_PETSC = VARX(16)
      DTOL_PETSC = VARX(17)
      MAXITS_PETSC = INT(VARX(18))
      DEALLOCATE( VARX,STAT=ISTAT )
      CHMSG = 'VARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for solution control variables  ---
!
      CALL ALLOC_SOLTN
!
!---  Initialize solution control variables  ---
!
      CALL INTLZ_SOLTN
!
!---  Read time stepping variables (duplicated across processors)  ---
!
      NVAR = LEPD
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TMPS,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = LEPD
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TMPE,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = LEPD
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TMPD,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = LEPD
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TMPX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = LEPD
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TMPN,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = LEPD
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TMPA,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = LEPD
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TMPC,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = LEPD
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,RSDM,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = LUK*(1+LWELL+LSPILL)
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,RSD,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = 20
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,WFMN,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read solution control integers (duplicated across processors)  ---
!
      NVAR = 39
      ALLOCATE( IVARX(1:NVAR),STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IVARX,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      IVRSN = IVARX(1)
      ISIC = IVARX(2)
      ICNV = IVARX(3)
      IEO = IVARX(4)
      ILES = IVARX(5)
      IOM = IVARX(6)
      ICODE = IVARX(7)
      IEQT = IVARX(8)
      IEQW = IVARX(9)
      IEQA = IVARX(10)
      IEQN = IVARX(11)
      IEQO = IVARX(12)
      IEQC = IVARX(13)
      IEQS = IVARX(14)
      IEQD = IVARX(15)
      IEQDO = IVARX(16)
      IEQHA = IVARX(17)
      IEQHN = IVARX(18)
      IEQHO = IVARX(19)
      IEQDA = IVARX(20)
      IAQU = IVARX(21)
      IGAS = IVARX(22)
      INAPL = IVARX(23)
      NEPD = IVARX(24)
      MEPD = IVARX(25)
      IEPD = IVARX(26)
      NRIMX = IVARX(27)
      NSTEP = IVARX(28)
      NRST = IVARX(29)
      NITER = IVARX(30)
      NTSR = IVARX(31)
      NGC = IVARX(32)
      MXSTEP = IVARX(33)
      IUNM = IVARX(34)
      IUNKG = IVARX(35)
      IUNS = IVARX(36)
      IUNK = IVARX(37)
      IUNMOL = IVARX(38)
      ISVC = IVARX(39)
      ISVF = 2*ISVC + 1
      DO M = 1,4
        M_ERR(M) = ''
      ENDDO
      DEALLOCATE( IVARX,STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read solution control integers (duplicated across processors)  ---
!
      NVAR = 100
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISLC,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      IEDLS = 1
      IF( ISLC(4).EQ.1 ) IEDLS = 3
!
!---  Read interfacial averaging indices
!     (duplicated across processors)  ---
!
      NVAR = 20
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IDMN,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Allocate memory for output control variables  ---
!
      CALL ALLOC_OUTPU
!
!---  Read output unit conversions
!     (duplicated across processors)  ---
!
      NVAR = 1
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,CNVTM,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = 1
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,CNVLN,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = LVPLOT
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,CNVPLOT,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = LVREF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,CNVREF,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      NVAR = 2*LSF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,CNVSF,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read output units
!     (duplicated across processors)  ---
!
      NVAR = 64
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,UNTM,NVAR,MPI_CHAR,
     &  STATUS,IERR )
      NVAR = 64
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,UNLN,NVAR,MPI_CHAR,
     &  STATUS,IERR )
      NVAR = 64*LVPLOT
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,UNPLOT,NVAR,MPI_CHAR,
     &  STATUS,IERR )
      NVAR = 64*LVREF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,UNREF,NVAR,MPI_CHAR,
     &  STATUS,IERR )
      NVAR = 64*2*LSF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,UNSF,NVAR,MPI_CHAR,
     &  STATUS,IERR )
!
!---  Read output control variables
!     (duplicated across processors)  ---
!
      NVAR = LPTM
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,PRTM,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read output control integer arrays
!     (duplicated across processors)  ---
!
      NVAR = LVPLOT
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IPLOT,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      NVAR = LVREF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IREF,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      NVAR = LREF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NDREF,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      NVAR = LSF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISFT,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      NVAR = LSF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISFF,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      NVAR = LSF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISFD,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      NVAR = LSF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISFGP,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read output control integers (duplicated across processors)  ---
!
      NVAR = 14
      ALLOCATE( IVARX(1:NVAR),STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IVARX,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      NPRTM = IVARX(1)
      NVPLOT = IVARX(2)
      NREF = IVARX(3)
      NVREF = IVARX(4)
      ICNO = IVARX(5)
      ICNS = IVARX(6)
      NSF = IVARX(7)
      NSFGP = IVARX(8)
      IHSF = IVARX(9)
      IHSF = 0
      IFQS = IVARX(10)
      IFQO = IVARX(11)
      ISGNS = IVARX(12)
      ISGNO = IVARX(13)
      ISGNP = IVARX(14)
!      PRINT *,'NSF = ',NSF,' ID = ',ID
!
!---  Read output character strings, duplicated over all processors  ---
!
      NVAR = 64*LVREF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,CHREF,NVAR,MPI_CHAR,
     &  STATUS,IERR )
!
!---  Read surface flux file names, duplicated over all processors  ---
!
      NVAR = 64*LSF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,FNSF,NVAR,MPI_CHAR,
     &  STATUS,IERR )
!      PRINT *,'FNSF = ',(FNSF(M),M=1,NSF),' ID = ',ID
!
!---  Read surface flux header character strings,
!     duplicated over all processors  ---
!
      NVAR = 64*2*LSF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,CHSF,NVAR,MPI_CHAR,
     &  STATUS,IERR )
!      PRINT *,'CHSF = ',((CHSF(L,M),L=1,2),M=1,NSF),' ID = ',ID
!
!---  Allocate memory for NSFN  ---
!
      ALLOCATE( NSFN(1:NP),STAT=ISTAT )
      CHMSG = 'NSFN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read local number of surface flux nodes on each processor,
!     duplicated over all processors  ---
!
      NVAR = NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NSFN,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      DO M = 1,NP
!        PRINT *,'NSFN(',M,') = ',NSFN(M),' ID = ',ID
!      ENDDO
!
!---  Determine offsets for reading ISFN and ISFS  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NSFN(I)
      ENDDO
      NCP = 0
      DO N = 1,NP
        NCP = NCP + NSFN(N)
      ENDDO
!      PRINT *,'NC = ',NC,' NCP = ',NCP,' ID = ',ID
!
!---  Allocate memory for ISFN  ---
!
      NSFNX = MAX( NSFN(ID+1),1 )
      ALLOCATE( ISFN(1:NSFNX),STAT=ISTAT )
      CHMSG = 'ISFN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for ISFS  ---
!
      NSFNX = MAX( NSFN(ID+1),1 )
      ALLOCATE( ISFS(1:NSFNX),STAT=ISTAT )
      CHMSG = 'ISFS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for ISFB  ---
!
      NSFNX = MAX( NSFN(ID+1),1 )
      ALLOCATE( ISFB(1:NSFNX),STAT=ISTAT )
      CHMSG = 'ISFB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read local surface-flux node number  ---
!
      NVAR = NSFN(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISFN,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'ISFN = ',(ISFN(M),M=1,NSFN(ID+1)),' ID = ',ID
!
!---  Read local surface-flux number  ---
!
      NVAR = NSFN(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISFS,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'ISFS = ',(ISFS(M),M=1,NSFN(ID+1)),' ID = ',ID
!
!---  Close the solu.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_SOLU_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_SORC_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Read binary input.bin file for source data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE SOURC
      USE SOLTN
      USE PROP
      USE OUTPU
      USE MPI
      USE GRID
      USE GLB_PAR
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
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_SORC_W'
!
!---  Open sorc.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'sorc.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Allocate memory for NSR  ---
!
      ALLOCATE( NSR(1:NP),STAT=ISTAT )
      CHMSG = 'NSR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read number of source nodes on each processor  ---
!
      NVAR = NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NSR,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Allocate local source arrays  ---
!
      CALL ALLOC_SOURC
!
!---  Initialize source variables  ---
!
      CALL INTLZ_SOURC
!
!---  Read source variables (duplicated across processors)  ---
!
      LX = 8+LSOLU+LSPT+LNGC
      NVAR = LX*LSTM*LSR
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,SRC,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!
!---  Read array of source nodes  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NSR(I)
      ENDDO
      NCP = 0
      DO N = 1,NP
        NCP = NCP + NSR(N)
      ENDDO
!
!---  Index array of source field nodes  ---
!
      NVAR = NSR(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISRN,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Index array of source number of time points  ---
!
      NVAR = NSR(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISRM,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Index array of source inputs  ---
!
      NVAR = NSR(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISRIN,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!
!---  Index array of source types  ---
!
      NVAR = NSR(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISRT,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!      DO N = 1,NSR(ID+1)
!        PRINT *,'ISRT(',N,') = ',ISRT(N),' ID = ',ID
!      ENDDO
!
!---  Index array of source types  ---
!
      NVAR = 13*NSR(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISRDM,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!      DO N = 1,NSR(ID+1)
!        DO M = 1,13
!          PRINT *,'ISRDM(',M,',',N,') = ',ISRDM(M,N),' ID = ',ID
!        ENDDO
!      ENDDO
!
!---  Close the sorc.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_SORC_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_STATE_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Read binary input.bin file for state condition data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE SOLTN
      USE PROP
      USE MPI
      USE HYST
      USE GRID
      USE GLB_PAR
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
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_STATE_W'
!
!---  Open state.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'state.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Allocate local memory state condition arrays 
!     (including ghost cells)  ---
!
      CALL ALLOC_FDVP
      CALL ALLOC_FDVH
!
!---  Initialize global array memory for general field variables  ---
!
      CALL INTLZ_FDVP
      CALL INTLZ_FDVH
!
!---  Allocate array memory for hysteretic k-s-P function 
!     variables  ---
!
      CALL ALLOC_HYST
!
!---  Initialize global array memory for Hysteretic k-s-P function 
!     variables  ---
!
      CALL INTLZ_HYST
!
!---  Allocate local temporary state condition arrays 
!     (including ghost cells)  ---
!
      ALLOCATE( VARX(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'VARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IVARX(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Set local starting point for local copies of nodal variables  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NFCGC(I)
      ENDDO
!
!---  Read local copies of T array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NC*NBYTR + NBYTB
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        T(2,N) = VARX(N)
      ENDDO
!      PRINT *,'T(2,',ND(11),') = ',T(2,11),' ID = ',ID
!
!---  Read local copies of PL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        PL(2,N) = VARX(N)
      ENDDO
!      PRINT *,'PL(2,',ND(11),') = ',PL(2,11),' ID = ',ID
!
!---  Read local copies of PG array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        PG(2,N) = VARX(N)
      ENDDO
!      PRINT *,'PG(2,',ND(11),') = ',PG(2,11),' ID = ',ID
!
!---  Read local copies of SG array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        SG(2,N) = VARX(N)
      ENDDO
!      PRINT *,'SG(2,',ND(11),') = ',SG(2,11),' ID = ',ID
!
!---  Read local copies of SL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        SL(2,N) = VARX(N)
      ENDDO
!      PRINT *,'SL(2,',ND(11),') = ',SL(2,11),' ID = ',ID
!
!---  Read local copies of SGT array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        SGT(2,N) = VARX(N)
      ENDDO
!      PRINT *,'SGT(2,',ND(11),') = ',SGT(2,11),' ID = ',ID
!
!---  Read local copies of ASLMIN array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        ASLMIN(2,N) = VARX(N)
      ENDDO
!      PRINT *,'ASLMIN(2,',ND(11),') = ',ASLMIN(2,11),' ID = ',ID
!
!---  Read local copies of PORD array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        PORD(2,N) = VARX(N)
      ENDDO
!      PRINT *,'PORD(2,',ND(11),') = ',PORD(2,11),' ID = ',ID
!
!---  Read local copies of PORT array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        PORT(2,N) = VARX(N)
      ENDDO
!      PRINT *,'PORT(2,',ND(11),') = ',PORT(2,11),' ID = ',ID
!
!---  Read local copies of NPHAZ array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IVARX,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        NPHAZ(2,N) = IVARX(N)
      ENDDO
!      PRINT *,'NPHAZ(2,',ND(11),') = ',NPHAZ(2,11),' ID = ',ID
!
!---  Deallocate local temporary state condition arrays 
!     (including ghost cells)  ---
!
      DEALLOCATE( VARX,STAT=ISTAT )
      CHMSG = 'VARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
      DEALLOCATE( IVARX,STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Close the state.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_STATE_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_TPORT_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Read binary input.bin file for solute transport data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE HYST
      USE GRID
      USE GLB_PAR
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
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_TPORT_W'
!
!---  Open tpor.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'tpor.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Allocate and initialize memory, global and local, for
!     solute transport arrays  ---
!
      CALL ALLOC_TRNSPT
      CALL INTLZ_TRNSPT
!
!---  Set local starting point for local copies of nodal variables  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NFCGC(I)
      ENDDO
!
!---  Allocate local temporary real and integer nodal arrays
!     (including ghost cells)  ---
!
      ALLOCATE( VARX(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'VARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IVARX(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read number of solutes (duplicated across processors)  ---
!
      NVAR = 1
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NSOLU,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'NSOLU = ',NSOLU,' ID = ',ID
!
!---  Read solute parameters, indices, and variables
!     if solutes are modeled  ---
!
      IF( NSOLU.GT.0 ) THEN
!
!---    Read solute names (duplicated across processors)  ---
!
        NVAR = 64*NSOLU
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SOLUT,NVAR,MPI_CHAR,
     &    STATUS,IERR )
!        PRINT *,'SOLUT = ',(SOLUT(NSL),NSL=1,NSOLU),' ID = ',ID
!
!---    Loop over the number of solutes  ---
!
        DO NSL = 1,NSOLU
!
!---      Read local copies of C array (including ghost cells)  ---
!
          NVAR = NFCGC(ID+1)
          OFFSET = IOFFSET + NC*NBYTR + NBYTB
          IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
          CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &      STATUS,IERR )
          DO N = 1,NFCGC(ID+1)
            C(N,NSL) = VARX(N)
          ENDDO
!          PRINT *,'C(2,',NSL,') = ',C(2,NSL),'C(3,',NSL,') = ',C(3,NSL),
!     &      'C(4,',NSL,') = ',C(4,NSL),' ID = ',ID
        ENDDO
!
!---    Read solute aqueous diffusion coefficient
!       (duplicated across processors)  ---
!
        NVAR = NSOLU
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SMDL,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        PRINT *,'SMDL = ',(SMDL(NSL),NSL=1,NSOLU),' ID = ',ID
!
!---    Read solute half-life, s
!       (duplicated across processors)  ---
!
        NVAR = NSOLU
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,HLF,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        PRINT *,'HLF = ',(HLF(NSL),NSL=1,NSOLU),' ID = ',ID
!
!---    Read solute chain-decay fraction
!       (duplicated across processors)  ---
!
        NVAR = NSOLU*NSOLU
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,CHDF,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        PRINT *,'CHDF = ',((CHDF(MSL,NSL),MSL=1,NSOLU),NSL=1,NSOLU),
!     &    ' ID = ',ID
!
!---    Read number of Bateman chain decay series
!       (duplicated across processors)  ---
!
        NVAR = 1
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,NBCDS,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'NBCDS = ',NBCDS,' ID = ',ID
!
!---    Read Bateman chain decay series index
!       (duplicated across processors)  ---
!
        NVAR = LCDC+LSOLU
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCDS,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'IBCDS = ',(IBCDS(M),M=1,LCDC+LSOLU),' ID = ',ID
!
!---    Read Bateman number of solutes in chain decay path
!       (duplicated across processors)  ---
!
        NVAR = LCDC
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,NBCDP,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'NBCDP = ',(NBCDP(M),M=1,LCDC),' ID = ',ID
!
!---    Read Bateman chain decay path indices
!       (duplicated across processors)  ---
!
        NVAR = LCDS*LCDP*LCDC
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCDP,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'IBCDP = ',
!     &    (((IBCDP(K,L,M),K=1,LCDS),L=1,LCDP),M=1,LCDC),' ID = ',ID
!
!---    Read Courant number calculation index
!       (duplicated across processors)  ---
!
        NVAR = 1
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ICRNT,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'ICRNT = ',ICRNT,' ID = ',ID
!        PRINT *,'IBCDS = ',(IBCDS(M),M=1,NSOLU+NBCDS),' ID = ',ID
!
!---    Read maximum Courant number
!       (duplicated across processors)  ---
!
        NVAR = 1
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,CRNTMXT,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        PRINT *,'CRNTMXT = ',CRNTMXT,' ID = ',ID
!
!---    Loop over the number of solutes  ---
!
        DO NSL = 1,NSOLU
!
!---      Loop over the number of solid-aqueous partition
!         coefficient indices  ---
!
          DO M = 1,5
!
!---        Read local copies of solute solid-aqueous partition
!           coefficient (including ghost cells)  ---
!
            NVAR = NFCGC(ID+1)
            OFFSET = IOFFSET + NC*NBYTR + NBYTB
            IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
            CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &        STATUS,IERR )
            DO N = 1,NFCGC(ID+1)
              PCSL(M,N,NSL) = VARX(N)
            ENDDO
!            PRINT *,'PCSL(',M,',',ND(11),',',NSL,') = ',
!     &        PCSL(M,11,NSL),' ID = ',ID
          ENDDO
        ENDDO
!
!---    Read local copies of longitudinal dispersivity
!       (including ghost cells)  ---
!
        NVAR = NFCGC(ID+1)
        OFFSET = IOFFSET + NC*NBYTR + NBYTB
        IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,DISPL,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        PRINT *,'DISPL(1) = ',DISPL(1),' ID = ',ID
!
!---    Read local copies of transverse dispersivity
!       (including ghost cells)  ---
!
        NVAR = NFCGC(ID+1)
        OFFSET = IOFFSET + NC*NBYTR + NBYTB
        IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,DISPT,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        PRINT *,'DISPT(1) = ',DISPT(1),' ID = ',ID
      ENDIF
!
!---  Loop over the number of solutes and reactive species  ---
!
      DO M = 1,LSOLU+LSPT
!
!---    Read local copies of solute or reactive species initial
!       condition indices (including ghost cells)  ---
!
        NVAR = NFCGC(ID+1)
        OFFSET = IOFFSET + NC*NBYTI + NBYTB
        IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IVARX,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
        DO N = 1,NFCGC(ID+1)
          ICT(N,M) = IVARX(N)
        ENDDO
!        PRINT *,'ICT(',ND(11),',',M,') = ',ICT(11,M),' ID = ',ID
      ENDDO
!
!---  Deallocate local temporary real and integer nodal arrays
!     (including ghost cells)  ---
!
      DEALLOCATE( VARX,STAT=ISTAT )
      CHMSG = 'VARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
      DEALLOCATE( IVARX,STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Close the tpor.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_TPORT_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKLS_W( HDGL,RKLX,SLP,SLPF,SLPM,SLX,N,IPHX,M )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Compute the aqueous relative permeability from the
!     aqueous saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
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
      SUB_LOG(ISUB_LOG) = '/RKLS_W'
      ESLZ = MAX( ESLX,0.D+0 )
!
!---  Constant relative permeability function  ---
!
      IF( MOD( IRPL(N),100 ).EQ.0 ) THEN
!
!---    Single-pressure dual-porosity saturation function  ---
!
        IF( ISCHR(N).EQ.3 .OR. ISCHR(N).EQ.4 ) THEN
          RKLM = RPLC(2,N)
          RKLF = RPLC(1,N)
          RKLX = ( PERM(4,N)*RKLM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKLF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
!
!---    Triple-curve saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.301 .OR. ISCHR(N).EQ.302 ) THEN
          IF( IPHX.EQ.2 ) THEN
            RKLX = 0.D+0
            IF( SLX.GE.SCHR(4,N) ) RKLX = RPLC(1,N)
          ELSE
            RKLX = RPLC(2,N)
          ENDIF
        ENDIF
!
!---  Mualem-irreducible porosity distribution function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.21 ) THEN
        SLPX = (SLP-SLP*SCHR(4,N)+SCHR(4,N)-RPLC(1,N))/
     &    (1.D+0-RPLC(1,N))
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(N).EQ.1  ) THEN
          RKLX = SQRT(SLPX)*((1.D+0-(1.D+0-SLPX**(1.D+0/RPLC(2,N)))
     &      **RPLC(2,N))**2)
!
!---    Brooks-Corey saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.2 ) THEN
          RKLX = SLPX**(2.5D+0 + 2.0D+0/RPLC(2,N))
        ENDIF
!
!---  Mualem-Anisotropy reference porosity distribution function  ---
!
      ELSEIF( IRPL(N).EQ.301 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(N).EQ.1 .OR. ISCHR(N).EQ.101 ) THEN
          RKLX = (SLP**RPLC(3,N))*((1.D+0-(1.D+0-SLP**
     &      (1.D+0/RPLC(2,N)))**RPLC(2,N))**2)
!
!---    Brooks-Corey saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 ) THEN
          RKLX = (SLP**RPLC(3,N))*SLP**(2.D+0 + 2.0D+0/RPLC(2,N))
        ENDIF
!
!---  Mualem porosity distribution function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.1 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(N).EQ.1 .OR. ISCHR(N).EQ.101 ) THEN
          RKLX = SQRT(SLP)*((1.D+0-(1.D+0-SLP**(1.D+0/RPLC(2,N)))
     &      **RPLC(2,N))**2)
!
!---    Brooks-Corey saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 ) THEN
          RKLX = SLP**(2.5D+0 + 2.0D+0/RPLC(2,N))
!
!---    van Genuchten triple-curve saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.301) THEN
          IF( IPHX.EQ.2 ) THEN
            RKLX = 0.D+0
            IF( SLX.GT.SCHR(4,N) ) RKLX = SQRT(SLP)*
     &        ((1.D+0-(1.D+0-SLP**(1.D+0/RPLC(1,N)))**RPLC(1,N))**2)
          ELSE
            RKLX = SQRT(SLP)*((1.D+0-(1.D+0-SLP**(1.D+0/RPLC(2,N)))
     &        **RPLC(2,N))**2)
          ENDIF
!
!---    Brooks-Corey triple-curve saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.302 ) THEN
          IF( IPHX.EQ.2 ) THEN
            RKLX = 0.D+0
            IF( SLX.GT.SCHR(4,N) )
     &        RKLX = SLP**(2.5D+0 + 2.0D+0/RPLC(1,N))
          ELSE
            RKLX = SLP**(2.5D+0 + 2.0D+0/RPLC(2,N))
          ENDIF
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
        ELSEIF( ISCHR(N).EQ.3  ) THEN
          RKLM = 1.D+0 - (SLPM**(1.D+0/RPLC(2,N)))
          RKLM = SQRT(SLPM)*((1.D+0-(RKLM**RPLC(2,N)))**2)
          RKLF = 1.D+0 - (SLPF**(1.D+0/RPLC(1,N)))
          RKLF = SQRT(SLPF)*((1.D+0-(RKLF**RPLC(1,N)))**2)
          RKLX = ( PERM(4,N)*RKLM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKLF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
!
!---    Single-pressure dual-porosity Brooks and Corey  ---
!
        ELSEIF( ISCHR(N).EQ.4  ) THEN
          RKLM = SLPM**(2.5D+0 + 2.0D+0/RPLC(2,N))
          RKLF = SLPF**(2.5D+0 + 2.0D+0/RPLC(1,N))
          RKLX = ( PERM(4,N)*RKLM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKLF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
        ENDIF
!
!---  Burdine porosity distribution function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.2 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(N).EQ.1 .OR. ISCHR(N).EQ.101 ) THEN
          RKLX = (SLP**2)*(1.D+0-(1.D+0-SLP**(1.D+0/RPLC(2,N)))
     &        **RPLC(2,N))
!
!---    Brooks-Corey saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 ) THEN
          RKLX = SLP**(3.0D+0 + 2.0D+0/RPLC(2,N))
!
!---    Triple-curve van Genuchten saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.301 ) THEN
          IF( IPHX.EQ.2 ) THEN
            RKLX = 0.D+0
            IF( SLX.GE.SCHR(4,N) ) RKLX = (SLP**2)*(1.D+0-
     &        (1.D+0-SLP**(1.D+0/RPLC(1,N)))**RPLC(1,N))
          ELSE
            RKLX = (SLP**2)*(1.D+0-(1.D+0-SLP**(1.D+0/RPLC(2,N)))
     &        **RPLC(2,N))
          ENDIF
!
!---    Triple-curve Brooks-Corey saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.302 ) THEN
          IF( IPHX.EQ.2 ) THEN
            RKLX = 0.D+0
            IF( SLX.GE.SCHR(4,N) )
     &        RKLX = SLP**(3.0D+0 + 2.0D+0/RPLC(1,N))
          ELSE
            RKLX = SLP**(3.0D+0 + 2.0D+0/RPLC(2,N))
          ENDIF
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
        ELSEIF( ISCHR(N).EQ.3 ) THEN
          RKLM = (SLPM**2)*(1.D+0-(1.D+0-SLPM**(1.D+0/RPLC(2,N)))
     &      **RPLC(2,N))
          RKLF = (SLPF**2)*(1.D+0-(1.D+0-SLPF**(1.D+0/RPLC(1,N)))
     &      **RPLC(1,N))
          RKLX = ( PERM(4,N)*RKLM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKLF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
!
!---    Single-pressure dual-porosity Brooks and Corey  ---
!
        ELSEIF( ISCHR(N).EQ.4 ) THEN
          RKLM = SLPM**(3.0D+0 + 2.0D+0/RPLC(2,N))
          RKLF = SLPF**(3.0D+0 + 2.0D+0/RPLC(1,N))
          RKLX = ( PERM(4,N)*RKLM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKLF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
        ENDIF
!
!---  Corey relative permeability function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.3 ) THEN
        RKLX = SLP**4
!
!---  Fatt and Klikoff relative permeability function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.4 ) THEN
        RKLX = SLP**3
!
!---  Haverkamp relative permeability function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.5 ) THEN
        IF( HDGL.LE.RPLC(3,N) ) THEN
          RKLX = 1.D+0
        ELSE
          ALPHAX = RPLC(1,N)/RPLC(4,N)
          HDGLX = HDGL/RPLC(4,N)
          RKLX = ALPHAX/(ALPHAX +
     &      (HDGLX**RPLC(2,N)))
        ENDIF
!
!---  Touma and Vauclin relative permeability function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.6 ) THEN
        RKLX = RPLC(1,N)*(SLP**RPLC(2,N))
!
!---  Free Corey relative permeability function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.7 ) THEN
        SLRX = RPLC(3,N)
        SGRX = RPLC(4,N)
        SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SGRX),0.D+0 ),1.D+0 )
        RKLX = RPLC(1,N)*(SLPX**(RPLC(2,N)))
!
!---  Rijtema-Gardner modified exponential function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.9 ) THEN
        RKLX = EXP( RPLC(1,N)*HDGL + RPLC(2,N) )
!
!---  Tabular function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.10 ) THEN
        ITBX = 0
        IF( M.NE.2 ) ITBX = 1
        RKLX = FNTBLY( SLX,IRLTBL(1,N),IRLTBL(2,N),ITBX )
!
!---  Cubic-spline tabular function versus saturation  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.11 ) THEN
        IF( ID.EQ.0 ) PRINT *,'Cubic-spline tabular function ' // 
     &    'versus saturation not currently available.'
        CALL MPI_FINALIZE(IERR)
        STOP
!        ITX = 1
!        ITBX = 0
!        RKLX = FSPLNY( SLX,IRLTBLT(1,N,ITX),IRLTBLT(2,N,ITX) )
!
!---  Linear tabular function versus capillary head  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.12 ) THEN
        IF( ID.EQ.0 ) PRINT *,'Linear tabular function versus ' // 
     &    'capillary head not currently available.'
        CALL MPI_FINALIZE(IERR)
        STOP
!        ITX = 1
!        ITBX = 0
!        RKLX = FNTBLY( HDGL,IRLTBLT(1,N,ITX),IRLTBLT(2,N,ITX),
!     &     ITBX )
!
!---  Cubic-spline tabular function versus capillary head  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.13 ) THEN
        IF( ID.EQ.0 ) PRINT *,'Cubic-spline tabular function versus ' // 
     &    'capillary head not currently available.'
        CALL MPI_FINALIZE(IERR)
        STOP
!        ITX = 1
!        ITBX = 0
!        RKLX = FSPLNY( HDGL,IRLTBLT(1,N,ITX),IRLTBLT(2,N,ITX) )
!
!---  Linear tabular function versus log capillary head  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.14 ) THEN
        IF( ID.EQ.0 ) PRINT *,'Linear tabular function versus ' // 
     &    'log capillary head not currently available.'
        CALL MPI_FINALIZE(IERR)
        STOP
!        ITX = 1
!        ITBX = 0
!        HDGLX = LOG(HDGL)
!        RKLX = FNTBLY( HDGLX,IRLTBLT(1,N,ITX),IRLTBLT(2,N,ITX),
!     &    ITBX )
!
!---  Cubic-spline tabular function versus log capillary head  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.15 ) THEN
        IF( ID.EQ.0 ) PRINT *,'Cubic-spline tabular function versus ' // 
     &    'log capillary head not currently available.'
        CALL MPI_FINALIZE(IERR)
        STOP
!        ITX = 1
!        ITBX = 0
!        HDGLX = LOG(HDGL)
!        RKLX = FSPLNY( HDGLX,IRLTBLT(1,N,ITX),IRLTBLT(2,N,ITX) )
!
!---  Polynomial function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.19 ) THEN
!
!---    Convert head units for polynomial function basis  ---
!
        HDGLU = HDGL/RPLC(1,N)
        IF( HDGLU.LT.CPLY_RL(1,1,N) ) THEN
          RKLX = 1.D+0
        ELSEIF( HDGLU.GE.CPLY_RL(2,NPLY_RL(N),N) ) THEN
          RKLX = 0.D+0
        ELSE
          DO NP = 1,NPLY_RL(N)
            IF( HDGLU.GE.CPLY_RL(1,NP,N) .AND.
     &        HDGLU.LT.CPLY_RL(2,NP,N) ) THEN
                RKLX = 0.D+0
                NPOLYC = LPOLYC
                DO NC = 5,NPOLYC
                  RKLX = RKLX + CPLY_RL(NC,NP,N)*
     &              (LOG10(HDGLU)**(NC-5))
                ENDDO
                EXIT
            ENDIF
          ENDDO
!
!---  Normalize absolute conductivity with saturated
!     conductivity  ---
!
          RKLX = (1.D+1**RKLX)/RPLC(2,N)
        ENDIF
!
!---  Modified-Mualem porosity distribution function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.22 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(N).EQ.1 .OR. ISCHR(N).EQ.101 ) THEN
          RKLX = (SLP**RPLC(1,N))*((1.D+0-(1.D+0-SLP**
     &      (1.D+0/RPLC(2,N)))**RPLC(2,N))**2)
!
!---    Brooks-Corey saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 ) THEN
          RKLX = (SLP**RPLC(1,N))*SLP**(2.D+0 + 2.0D+0/RPLC(2,N))
        ENDIF
      ENDIF
!
!---  Polmann anisotropy permeability function  ---
!
      IF( IRPL(N).GE.100 .AND. IRPL(N).LT.200 ) THEN
        PSI = HDGL*1.D+2
        SKLX = RPLC(5,N) - RPLC(10,N)*PSI -
     &    RPLC(6,N)*RPLC(9,N)*( RPLC(7,N) -
     &    (RPLC(7,N)**2)*PSI - (RPLC(8,N)**2)*PSI)/
     &    (1.D+0 + RPLC(10,N)*RPLC(9,N))
        SIGMA = RPLC(6,N)*(((1.D+0 - RPLC(7,N)*PSI)**2) +
     &    (RPLC(8,N)**2)*(PSI**2))/
     &    (1.D+0 + RPLC(10,N)*RPLC(9,N))
        SKHX = EXP( SKLX + 5.D-1*SIGMA )
        SKVX = EXP( SKLX - 5.D-1*SIGMA )
        ANISOX = MIN( MAX( SKHX/SKVX,0.D+0 ),RPLC(11,N) )
        ANISOX = MAX( ANISOX,RPLC(12,N) )
        RKLX = RKLX*ANISOX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RKLS_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Compute the maximum relative residuals
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE OUTPU
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(LUK*(1+LWELL+LSPILL)) :: RSDLX	
      REAL*8, DIMENSION(3) :: VARX	
      INTEGER, DIMENSION(LUK*(1+LWELL+LSPILL)) :: NSDLX,NPHLX,NPHX
      INTEGER, DIMENSION(LUK*(1+LWELL+LSPILL)) :: IDLX,IDX
      INTEGER, DIMENSION(NFCGC(ID+1)) :: IRSDX
      INTEGER, DIMENSION(3) :: IVARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER*64 PH_CND(2)
!
!----------------------Data Statements---------------------------------!
!
      DATA PH_CND /'Saturated','Unsaturated'/
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL_W'
!
!---  Zero local and global maximum residuals  ---
!
      DO M = 1,ISVC
        RSD(M) = 0.D+0
        RSDLX(M) = 0.D+0
        NSD(M) = 0
        NSDLX(M) = 0
        NPHLX(M) = 0
      ENDDO
!
!---  Loop over local nodes  ---
!
      NMD = 0
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
        MPW = NMD + IEQW
        N_DB = ND(N)
!
!---    Increment equation counter for next active node  ---
!
        NMD = NMD + ISVC
!
!---    Skip selected nodes in the residual calculation 
!       (not implemented)  ---
!
        IF( ISKP(N).EQ.1 ) CYCLE
!
!---    Water saturated system prior to iteration  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          ACP = PORD(2,N)*RHOL(2,N)*SL(2,N)*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPW))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/ACP) )
          IF( RSDX.GT.RSDLX(IEQW) ) THEN
            RSDLX(IEQW) = RSDX
            NSDLX(IEQW) = N
            NPHLX(IEQW) = NPHAZ(2,N)
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---    Water-gas system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          ACP = PORD(2,N)*RHOL(2,N)*SL(2,N)*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPW))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/ACP) )
          IF( RSDX.GT.RSDLX(IEQW) ) THEN
            RSDLX(IEQW) = RSDX
            NSDLX(IEQW) = N
            NPHLX(IEQW) = NPHAZ(2,N)
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
        ENDIF
      ENDDO
!
!---  Maximum global residuals  ---
!
      CALL MPI_ALLREDUCE( RSDLX,RSD,ISVC,MPI_REAL8,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
!
!---  Identify processor with maximum residual  ---
!
      DO M = 1,ISVC
        IDLX(M) = -1
        IF( ABS((RSDLX(M)-RSD(M))/EPSL).LT.EPSL ) IDLX(M) = ID
        IF( ID.EQ.IDLX(M) ) THEN
          NSD(M) = ND( NSDLX(M) )
          NPHX(M) = NPHLX(M)
        ELSE
          NSD(M) = 0
          NPHX(M) = 0
        ENDIF
      ENDDO
      CALL MPI_ALLREDUCE( IDLX,IDX,ISVC,MPI_INTEGER,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
      DO M = 1,ISVC
        CALL MPI_BCAST( NSD(M),1,MPI_INTEGER,IDX(M),MPI_COMM_WORLD,
     &    IERR )
        CALL MPI_BCAST( NPHX(M),1,MPI_INTEGER,IDX(M),MPI_COMM_WORLD,
     &    IERR )
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
!---  Unconverged solution and Newton-Raphson iteration 
!     limit exceeded  ---
!
      IF( ICNV.EQ.1 ) THEN
        OFFSET = IOFFSET_REF
        IF( ID.EQ.0 ) THEN
          IF( RSDX.GE.1.D+2 .AND. RSD_CW.LT.1.D+2 ) THEN
            PRINT *,'           ---  Excessive Residual  ---'
            IVARX(1) = -5
          ELSEIF( RSD_CW.LT.RSDMX ) THEN
            PRINT *,'           ---  Convergence Failure  ---'
            IVARX(1) = -1
          ENDIF
        ENDIF
!
!---    Debug printing to the screen  ---
!
        IF( ID.EQ.0 ) THEN
          NX = NSD(IEQW)
          IF( NX.GT.0 ) THEN
            NPX = NPHX(IEQW)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            PRINT *,'  Water Mass Equation Maximum Residual = ',
     &        RSD(IEQW),': Node = ',NX,': Phase Condition = ',
     &        PH_CND(NPX)(1:NCHX)
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
          NCH = INDEX( UNTM(1:),'  ' ) - 1
          IF( ID.EQ.0 ) PRINT *,'  Time Step Reduced From ',
     &      DTX*CNVTM,UNTM(1:NCH),' to ',DT*CNVTM,UNTM(1:NCH)
          DO N = 1,NFCGC(ID+1)
            PL(2,N) = PL(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
          ENDDO





          IVARX(1) = -1
          VARX(2) = DTX*CNVTM
          VARX(3) = DT*CNVTM
!
!---    Number of time step reductions failure: stop simulation  ---
!
        ELSE
          IF( ID.EQ.0 ) PRINT *,'          ---  Time Step Reduction ' // 
     &        'Limit Exceeded  ---'
          ICNV = 4
!
!---      Write a time-step reduction limit exceeded index in the 
!         NSTEP location to output.bin.  ---
!
          IVARX(1) = -2
          NVAR = 1
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &      MPI_INTEGER,STATUS,IERR)
          OFFSET = OFFSET + NVAR*NBYTI
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
          RETURN
        ENDIF
!
!---    Write a convergence failure index in the NSTEP location
!       plus write the global node numbers and phase condition indices
!       for the location of maximum residuals for the water, air, 
!       and salt equations to output.bin  ---
!
        NVAR = 3
        IVARX(2) = NSD(IEQW)
        IVARX(3) = NPHX(IEQW)
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &    MPI_INTEGER,STATUS,IERR)
        OFFSET = OFFSET + NVAR*NBYTI
        IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
!
!---    Write maximum residuals for the water, air, and salt 
!       equations and time step reductions to output.bin  ---
!
        NVAR = 3
        VARX(1) = RSD(IEQW)
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,VARX,NVAR,
     &    MPI_REAL8,STATUS,IERR)
        OFFSET = OFFSET + NVAR*NBYTR
        IOFFSET_REF = IOFFSET_REF + NVAR*NBYTR
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDL_W group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORC_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!

!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Compute source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PROP
      USE HYST
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
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORC_W'
!
!---  Zero source terms  ---
!
      DO NS = 1,NSR(ID+1)
        N = ISRN(NS)
        DO M = 2,ISVC+2
          SRCW(M,N) = 0.D+0
        ENDDO
      ENDDO
!
!---  Loop over sources  ---
!
      DO NS = 1,NSR(ID+1)
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = ISRIN(NS)
        IF( TMZ.LE.SRC(1,1,MB) ) CYCLE
        IF( ISRM(NS).EQ.1 ) THEN
          DO N = 1,8+NSOLU
            SRX(N) = SRC(N,1,MB)
          ENDDO

        ELSE
          IFIND = 0
          DO M = 2,ISRM(NS)
            IF( TMZ.LE.SRC(1,M,MB) ) THEN
              DTSR = MIN( SRC(1,M,MB)-TMZ,DT )
              TFSR = (TMZ-0.5D+0*DTSR-SRC(1,M-1,MB))/
     &          (SRC(1,M,MB)-SRC(1,M-1,MB))
              DO N = 1,8+NSOLU
                SRX(N) = SRC(N,M-1,MB) + 
     &            TFSR*(SRC(N,M,MB)-SRC(N,M-1,MB))
              ENDDO
!
!---          Stair step the slug and pulse well sources  ---
!
              IF( ISRT(NS).GE.24 .AND. ISRT(NS).LE.27 ) THEN
                DO N = 2,8+NSOLU
                  SRX(N) = SRC(N,M-1,NS)
                ENDDO
              ENDIF
              IFIND = 1
              EXIT
            ENDIF

          ENDDO
          IF( IFIND.EQ.0 ) CYCLE
        ENDIF



!
!---    Local node associated with processor-dependent source  ---
!
        N = ISRN(NS)
        N_DB = ND(N)
!
!---    Loop over increment indices  ---
!
        DO M = 2,ISVC+2
!
!---      Aqueous volumetric source  ---
!
          IF( ISRT(NS).EQ.2 ) THEN
            SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOL(M,N)
!
!---      Aqueous volumetric density source  ---
!
          ELSEIF( ISRT(NS).EQ.3 ) THEN
            SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOL(M,N)*VOL(N)
!
!---      Aqueous mass source  ---
!
          ELSEIF( ISRT(NS).EQ.4 ) THEN
            SRCW(M,N) = SRCW(M,N) + SRX(4)
!
!---      Aqueous mass density source  ---
!
          ELSEIF( ISRT(NS).EQ.5 ) THEN
            SRCW(M,N) = SRCW(M,N) + SRX(4)*VOL(N)
!
!---      Z-direction vertical injection well  ---
!
          ELSEIF( ISRT(NS).EQ.31 ) THEN
!
!---        Geometric factors  ---
!
            RDW = SRX(3)
            RDE = SQRT( AFZ(1,N)/GPI/SRX(4) )
            ACWX = 2.D+0*GPI*RDW*DZGF(N)
            DRD2 = (RDE**2-RDW**2)
!
!---        Well pressure  ---
!
            IF( M.EQ.2 ) THEN
              IF( K.EQ.ISRDM(5,NS) ) THEN
                PLWX = SRX(2)
                PX = PLWX+PATM
                CALL WATLQD(T(2,N),PX,RHOLX)
              ELSE
                PX = PLWX+PATM
                CALL WATLQD(T(2,N),PX,RHOLX)
                NX = ICM(1,N)
                IF( NX.NE.0 ) THEN
                  GB = (ZP(N)-ZP(NX))*GRAVZ
                  PLWX = PLWX - RHOLX*GB
                ENDIF
              ENDIF
            ENDIF
            IF( (PERM(1,N)/EPSL).GT.EPSL ) THEN
              IF( (PERM(2,N)/EPSL).GT.EPSL ) THEN
                PERMX = SQRT( PERM(1,N)*PERM(2,N) )
              ELSE
                PERMX = PERM(1,N)
              ENDIF
            ELSE
              PERMX = PERM(2,N)
            ENDIF
            HCLX = 2.D+0*GPI*PERMX*DRD2*DZGF(N)/
     &        (VISL(M,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
!
!---        Injection  ---
!
            IF( PLWX-PL(M,N).GT.EPSL ) THEN
              QLX = SRX(4)*(PLWX-PL(M,N))*HCLX
              SRCW(M,N) = SRCW(M,N) + QLX*RHOLX
!
!---          Withdrawl  ---
!
            ELSE
              QLX = SRX(4)*(PLWX-PL(M,N))*RKL(M,N)*HCLX
              SRCW(M,N) = SRCW(M,N) + QLX*RHOL(M,N)
            ENDIF
!
!---      X-direction horizontal injection well  ---
!
          ELSEIF( ISRT(NS).EQ.32 ) THEN
!
!---        Geometric factors  ---
!
            RDW = SRX(3)
            RDE = SQRT( AFX(1,N)/GPI/SRX(4) )
            ACWX = 2.D+0*GPI*RDW*DXGF(N)
            DRD2 = (RDE**2-RDW**2)
!
!---        Well pressure  ---
!
            IF( M.EQ.2 ) THEN
              IF( I.EQ.ISRDM(1,NS) ) THEN
                PLWX = SRX(2)
                PX = PLWX+PATM
                CALL WATLQD(T(2,N),PX,RHOLX)
              ELSE
                PX = PLWX+PATM
                CALL WATLQD(T(2,N),PX,RHOLX)
                NX = ICM(3,N)
                IF( NX.NE.0 ) THEN
                  GB = (XP(N)-XP(NX))*GRAVX
                  PLWX = PLWX - RHOLX*GB
                ENDIF
              ENDIF
            ENDIF
            IF( (PERM(2,N)/EPSL).GT.EPSL ) THEN
              IF( (PERM(3,N)/EPSL).GT.EPSL ) THEN
                PERMX = SQRT( PERM(2,N)*PERM(3,N) )
              ELSE
                PERMX = PERM(2,N)
              ENDIF
            ELSE
              PERMX = PERM(3,N)
            ENDIF
            HCLX = 2.D+0*GPI*PERMX*DRD2*DXGF(N)/
     &        (VISL(M,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
!
!---        Injection  ---
!
            IF( PLWX-PL(M,N).GT.EPSL ) THEN
              QLX = SRX(4)*(PLWX-PL(M,N))*HCLX
              SRCW(M,N) = SRCW(M,N) + QLX*RHOLX
!
!---          Withdrawl  ---
!
            ELSE
              QLX = SRX(4)*(PLWX-PL(M,N))*RKL(M,N)*HCLX
              SRCW(M,N) = SRCW(M,N) + QLX*RHOL(M,N)
            ENDIF
!
!---      Y-direction horizontal injection well  ---
!
          ELSEIF( ISRT(NS).EQ.33 ) THEN
!
!---        Geometric factors  ---
!
            RDW = SRX(3)
            RDE = SQRT( AFY(1,N)/GPI/SRX(4) )
            ACWX = 2.D+0*GPI*RDW*DYGF(N)
            DRD2 = (RDE**2-RDW**2)
!
!---        Well pressure  ---
!
            IF( M.EQ.2 ) THEN
              IF( J.EQ.ISRDM(3,NS) ) THEN
                PLWX = SRX(2)
                PX = PLWX+PATM
                CALL WATLQD(T(2,N),PX,RHOLX)
              ELSE
                PX = PLWX+PATM
                CALL WATLQD(T(2,N),PX,RHOLX)
                NX = ICM(2,N)
                IF( NX.NE.0 ) THEN
                  GB = (YP(N)-YP(NX))*RP(I)*GRAVY
                  PLWX = PLWX - RHOLX*GB
                ENDIF
              ENDIF
            ENDIF
            IF( (PERM(1,N)/EPSL).GT.EPSL ) THEN
              IF( (PERM(3,N)/EPSL).GT.EPSL ) THEN
                PERMX = SQRT( PERM(1,N)*PERM(3,N) )
              ELSE
                PERMX = PERM(1,N)
              ENDIF
            ELSE
              PERMX = PERM(3,N)
            ENDIF
            HCLX = 2.D+0*GPI*PERMX*DRD2*DYGF(N)/
     &        (VISL(M,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
!
!---        Injection  ---
!
            IF( PLWX-PL(M,N).GT.EPSL ) THEN
              QLX = SRX(4)*(PLWX-PL(M,N))*HCLX
              SRCW(M,N) = SRCW(M,N) + QLX*RHOLX
!
!---        Withdrawl  ---
!
            ELSE
              QLX = SRX(4)*(PLWX-PL(M,N))*RKL(M,N)*HCLX
              SRCW(M,N) = SRCW(M,N) + QLX*RHOL(M,N)
            ENDIF
          ENDIF

        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORC_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SP_W( N,M,PGX,PLX,SLX,RKLX,ASLX,ASLMINX,
     &  ASGTX,SGRMX,HDGL,SLP,SLPF,SLPM,INDX,IPHX,PGOX,PLOX,SLOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Compute the aqueous saturation from the gas/aqueous capillary
!     pressure, and compute the aqueous relative permeability from the
!     aqueous saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
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
      SUB_LOG(ISUB_LOG) = '/SP_W'
!
!---  van Genuchten saturation function  ---
!
      IF( ISCHR(N).EQ.1 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        CN = MAX( SCHR(3,N),SMALL )
        IF( SCHR(14,N).LE.ZERO ) THEN
          IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(14,N)
        ENDIF
        SLP = (1.D+0/(1.D+0 + (SCHR(1,N)*HDGL)**CN))**CM
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SMP = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        SLX = SLP*(1.D+0-SMP) + SMP
        ASLX = SLP
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(N).EQ.2 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        CL = MAX( SCHR(3,N),SMALL )
        IF( HDGL.LE.SCHR(1,N) ) THEN
          SLP = 1.D+0
        ELSE
          SLP = (SCHR(1,N)/HDGL)**CL
        ENDIF
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SMP = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        SLX = SLP*(1.D+0-SMP) + SMP
        ASLX = SLP
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Single-pressure dual-porosity
!     van Genuchten saturation function  ---
!
      ELSEIF( ISCHR(N).EQ.3 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        CN = MAX( SCHR(3,N),SMALL )
        IF( SCHR(14,N).LE.ZERO ) THEN
          IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(14,N)
        ENDIF
        SLPM = (1.D+0/(1.D+0 + (SCHR(1,N)*HDGL)**CN))**CM
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SMPM = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        SDPM(N) = SLPM*(1.D+0-SMPM) + SMPM
        CN = MAX( SCHR(6,N),SMALL )
        IF( SCHR(15,N).LE.ZERO ) THEN
          IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(15,N)
        ENDIF
        SLPF = (1.D+0/(1.D+0 + (SCHR(5,N)*HDGL)**CN))**CM
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SMPF = MAX( (1.D+0-HSCL)*SCHR(7,N),ZERO )
        SDPF(N) = SLPF*(1.D+0-SMPF) + SMPF
        PORDM = (1.D+0-POR(4,N))*POR(2,N)/
     &    ( POR(4,N) + (1.D+0-POR(4,N))*POR(2,N) + SMALL )
        PORDF = POR(4,N)/
     &    ( POR(4,N) + (1.D+0-POR(4,N))*POR(2,N) + SMALL )
        SLX = SDPF(N)*PORDF + SDPM(N)*PORDM
        ASLX = SLPF*PORDF + SLPM*PORDM
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Single-pressure dual-porosity
!     Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(N).EQ.4 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        CL = MAX( SCHR(3,N),SMALL )
        IF( HDGL.LE.SCHR(1,N) ) THEN
          SLPM = 1.D+0
        ELSE
          SLPM = (SCHR(1,N)/HDGL)**CL
        ENDIF
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SMPM = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        SDPM(N) = SLPM*(1.D+0-SMPM) + SMPM
        CL = MAX( SCHR(6,N),SMALL )
        IF( HDGL.LE.SCHR(5,N) ) THEN
          SLPF = 1.D+0
        ELSE
          SLPF = (SCHR(5,N)/HDGL)**CL
        ENDIF
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SMPF = MAX( (1.D+0-HSCL)*SCHR(7,N),ZERO )
        SDPF(N) = SLPF*(1.D+0-SMPF) + SMPF
        PORDM = (1.D+0-POR(4,N))*POR(2,N)/
     &    ( POR(4,N) + (1.D+0-POR(4,N))*POR(2,N) + SMALL )
        PORDF = POR(4,N)/
     &    ( POR(4,N) + (1.D+0-POR(4,N))*POR(2,N) + SMALL )
        SLX = SDPF(N)*PORDF + SDPM(N)*PORDM
        ASLX = SLPF*PORDF + SLPM*PORDM
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Haverkamp saturation function  ---
!
      ELSEIF( ABS(ISCHR(N)).EQ.5 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( HDGL.LE.SCHR(1,N) ) THEN
          SLP = 1.D+0
        ELSE
          ALPHAX = SCHR(2,N)/SCHR(5,N)
          IF( ISCHR(N).EQ.-5 ) THEN
            HDGLX = LOG(HDGL/SCHR(5,N))
          ELSE
            HDGLX = HDGL/SCHR(5,N)
          ENDIF
          SLP = ALPHAX/(ALPHAX
     &      + (HDGLX**SCHR(3,N)))
        ENDIF
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SMP = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        SLX = SLP*(1.D+0-SMP) + SMP
        ASLX = SLP
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Russo saturation function  ---
!
      ELSEIF( ISCHR(N).EQ.9 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        SLP = (EXP(-5.D-1*SCHR(1,N)*HDGL)*
     &    (1.D+0 + 5.D-1*SCHR(1,N)*HDGL))**(2.D+0/(SCHR(3,N)+2.D+0))
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SMP = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        SLX = SLP*(1.D+0-SMP) + SMP
        ASLX = SLP
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Linear or linear-log interpolation function  ---
!
      ELSEIF( ISCHR(N).EQ.10 .OR. ISCHR(N).EQ.12 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( ISCHR(N).EQ.12 ) HDGL = LOG(HDGL)
        ITBX = 0
        IF( M.NE.2 ) ITBX = 1
        SLP = FNTBLY( HDGL,ISLTBL(1,N),ISLTBL(2,N),ITBX )
        SLX = SLP
        SGX = MAX( 1.D+0-SLX,ZERO )
        ASLX = SLP
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Cubic-spline or cubic-spline-log interpolation function  ---
!
      ELSEIF( ISCHR(N).EQ.11 .OR. ISCHR(N).EQ.13 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( ISCHR(N).EQ.13 ) HDGL = LOG(HDGL)
        SLP = FSPLNY( HDGL,ISLTBL(1,N),ISLTBL(2,N) )
        SLX = SLP
        SGX = MAX( 1.D+0-SLX,ZERO )
        ASLX = SLP
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Polynomial function  ---
!
      ELSEIF( ISCHR(N).EQ.19 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
!
!---    Convert head units for polynomial function basis  ---
!
        HDGLU = HDGL/SCHR(1,N)
        IF( HDGLU.LT.CPLY_SL(1,1,N) ) THEN
          SLX = 1.D+0
        ELSEIF( HDGLU.GE.CPLY_SL(2,NPLY_SL(N),N) ) THEN
          SLX = 0.D+0
        ELSE
          L1: DO NP = 1,NPLY_SL(N)
            IF( HDGLU.GE.CPLY_SL(1,NP,N) .AND.
     &        HDGLU.LT.CPLY_SL(2,NP,N) ) THEN
                SLX = 0.D+0
                NPOLYC = LPOLYC
                L2: DO NC = 5,NPOLYC
                  SLX = SLX + CPLY_SL(NC,NP,N)*(LOG10(HDGLU)**(NC-5))
                ENDDO L2
                EXIT L1
            ENDIF
          ENDDO L1
        ENDIF
        SLP = SLX
        SGX = MAX( 1.D+0-SLX,0.D+0 )
        ASLX = SLP
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Cambridge function
!
      ELSEIF( ISCHR(N).EQ.41 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        CN = MAX( SCHR(3,N),SMALL )
        SLP = ((SCHR(1,N)-HDGL)/(SCHR(1,N)+SMALL))**(1.D+0/CN)
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(HDOD),ZERO )*REALX
        SMP = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        SLX = SLP*(1.D+0-SMP) + SMP
        ASLX = SLP
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  van Genuchten saturation function w/ gas entrapment  ---
!
      ELSEIF( ISCHR(N).EQ.101 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        CN = MAX( SCHR(3,N),SMALL )
        IF( SCHR(14,N).LE.ZERO ) THEN
          IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(14,N)
        ENDIF
        ASLX = (1.D+0/(1.D+0 + (SCHR(1,N)*HDGL)**CN))**CM
        IF( INDX.EQ.2 ) RETURN
        ASLM = MIN( ASLX,ASLMINX )
        IF( SGRMX.GT.EPSL ) THEN
          R = 1.D+0/SGRMX - 1.D+0
          ASGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &      (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
        ELSE
          ASGTX = 0.D+0
        ENDIF
        SLP = ASLX - ASGTX
        SMP = MAX( SCHR(4,N),ZERO )
        SLX = SLP*(1.D+0-SMP) + SMP
!
!---  Brooks and Corey saturation function w/ gas entrapment  ---
!
      ELSEIF( ISCHR(N).EQ.102 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        CL = MAX( SCHR(3,N),SMALL )
        IF( HDGL.LE.SCHR(1,N) ) THEN
          ASLX = 1.D+0
        ELSE
          ASLX = (SCHR(1,N)/HDGL)**CL
        ENDIF
        IF( INDX.EQ.2 ) RETURN
        ASLM = MIN( ASLX,ASLMINX )
        IF( SGRMX.GT.EPSL ) THEN
          R = 1.D+0/SGRMX - 1.D+0
          ASGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &      (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
        ELSE
          ASGTX = 0.D+0
        ENDIF
        SLP = ASLX - ASGTX
        SMP = MAX( SCHR(4,N),ZERO )
        SLX = SLP*(1.D+0-SMP) + SMP
!
!---  van Genuchten triple curve saturation function  ---
!
      ELSEIF( ISCHR(N).EQ.301 ) THEN
!
!---  Drainage scanning (including main drainage)  ---
!
        IF( IPHX.EQ.-1 ) THEN
          HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
          CN = MAX( SCHR(3,N),SMALL )
          IF( SCHR(14,N).LE.ZERO ) THEN
            IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
              CM = 1.D+0 - 2.D+0/CN
            ELSE
              CM = 1.D+0 - 1.D+0/CN
            ENDIF
          ELSE
            CM = SCHR(14,N)
          ENDIF
          SLP = (1.D+0/(1.D+0 + (SCHR(1,N)*HDGL)**CN))**CM
          SMP = SCHR(4,N)
          IF( SLOX.GT.EPSL ) THEN
            SLPO = (SLOX-SMP)/(1.D+0-SMP)
            HDGL = MAX( (PGOX-PLOX)/RHORL/GRAV,1.D-14 )
            CN = MAX( SCHR(3,N),SMALL )
            IF( SCHR(14,N).LE.ZERO ) THEN
              IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
                CM = 1.D+0 - 2.D+0/CN
              ELSE
                CM = 1.D+0 - 1.D+0/CN
              ENDIF
            ELSE
              CM = SCHR(14,N)
            ENDIF
            SLPHO = (1.D+0/(1.D+0 + (SCHR(1,N)*HDGL)**CN))**CM
            SLP = SLP*SLPO/(SLPHO+SMALL)
          ENDIF
          SLX = SLP*(1.D+0-SMP) + SMP
!
!---  Main wetting  ---
!
        ELSEIF( IPHX.EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
          CN = MAX( SCHR(5,N),SMALL )
          IF( SCHR(7,N).LE.ZERO ) THEN
            IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
              CM = 1.D+0 - 2.D+0/CN
            ELSE
              CM = 1.D+0 - 1.D+0/CN
            ENDIF
          ELSE
            CM = SCHR(7,N)
          ENDIF
          SLP = (1.D+0/(1.D+0 + (SCHR(2,N)*HDGL)**CN))**CM
          SMP = SCHR(6,N)
          SLX = SLP*(1.D+0-SMP) + SMP
!
!---  Wetting scanning (including boundary wetting scanning)  ---
!
        ELSE
          HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
          CN = MAX( SCHR(3,N),SMALL )
          IF( SCHR(14,N).LE.ZERO ) THEN
            IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
              CM = 1.D+0 - 2.D+0/CN
            ELSE
              CM = 1.D+0 - 1.D+0/CN
            ENDIF
          ELSE
            CM = SCHR(14,N)
          ENDIF
          SLP = (1.D+0/(1.D+0 + (SCHR(12,N)*HDGL)**CN))**CM
          SMP = SCHR(4,N)
          IF( SLOX.GT.EPSL ) THEN
            SLPO = (SLOX-SMP)/(1.D+0-SMP)
            SLPM = (SCHR(10,N)-SMP)/(1.D+0-SMP)
            HDGL = MAX( (PGOX-PLOX)/RHORL/GRAV,1.D-14 )
            CN = MAX( SCHR(3,N),SMALL )
            IF( SCHR(14,N).LE.ZERO ) THEN
              IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
                CM = 1.D+0 - 2.D+0/CN
              ELSE
                CM = 1.D+0 - 1.D+0/CN
              ENDIF
            ELSE
              CM = SCHR(14,N)
            ENDIF
            SLPHO = (1.D+0/(1.D+0 + (SCHR(12,N)*HDGL)**CN))**CM
            SLP = (SLP-SLPM)*(SLPO-SLPM)/(SLPHO-SLPM+SMALL) + SLPM
            SLP = MIN( SLP,SLPM )
          ENDIF
          SLX = SLP*(1.D+0-SMP) + SMP
        ENDIF
        ASLX = SLP
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
!
!---  Brooks and Corey triple curve saturation function  ---
!
      ELSEIF( ISCHR(N).EQ.302 ) THEN
!
!---  Drainage scanning (including main drainage)  ---
!
        IF( IPHX.EQ.-1 ) THEN
          HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
          CL = MAX( SCHR(3,N),SMALL )
          IF( HDGL.LE.SCHR(1,N) ) THEN
            SLP = 1.D+0
          ELSE
            SLP = (SCHR(1,N)/HDGL)**CL
          ENDIF
          SMP = SCHR(4,N)
          IF( SLOX.GT.EPSL ) THEN
            SLPO = (SLOX-SMP)/(1.D+0-SMP)
            HDGL = MAX( (PGOX-PLOX)/RHORL/GRAV,1.D-14 )
            CL = MAX( SCHR(3,N),SMALL )
            IF( HDGL.LE.SCHR(1,N) ) THEN
              SLPHO = 1.D+0
            ELSE
              SLPHO = (SCHR(1,N)/HDGL)**CL
            ENDIF
            SLP = SLP*SLPO/(SLPHO+SMALL)
          ENDIF
          SLX = SLP*(1.D+0-SMP) + SMP
!
!---  Main wetting  ---
!
        ELSEIF( IPHX.EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
          CL = MAX( SCHR(5,N),SMALL )
          IF( HDGL.LE.SCHR(2,N) ) THEN
            SLP = 1.D+0
          ELSE
            SLP = (SCHR(2,N)/HDGL)**CL
          ENDIF
          SMP = SCHR(6,N)
          SLX = SLP*(1.D+0-SMP) + SMP
!
!---  Wetting scanning (including boundary wetting scanning)  ---
!
        ELSE
          HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
          CL = MAX( SCHR(3,N),SMALL )
          IF( HDGL.LE.SCHR(12,N) ) THEN
            SLP = 1.D+0
          ELSE
            SLP = (SCHR(12,N)/HDGL)**CL
          ENDIF
          SMP = SCHR(4,N)
          IF( SLOX.GT.EPSL ) THEN
            SLPM = (SCHR(10,N)-SMP)/(1.D+0-SMP)
            SLPO = (SLOX-SMP)/(1.D+0-SMP)
            HDGL = MAX( (PGOX-PLOX)/RHORL/GRAV,1.D-14 )
            CL = MAX( SCHR(3,N),SMALL )
            IF( HDGL.LE.SCHR(12,N) ) THEN
              SLPHO = 1.D+0
            ELSE
              SLPHO = (SCHR(12,N)/HDGL)**CL
            ENDIF
            SLP = (SLP-SLPM)*(SLPO-SLPM)/(SLPHO-SLPM+SMALL) + SLPM
            SLP = MIN( SLP,SLPM )
          ENDIF
          SLX = SLP*(1.D+0-SMP) + SMP
        ENDIF
        ASLX = SLP
        ASGTX = 0.D+0
        ASLM = MIN( ASLX,ASLMINX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SP_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRPGL_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Compute the total trapping number for gas entrapment in the
!     aqueous phase.
!
!     Pennell, K.D., G.A. Pope, L.M. Abriola.  1996.
!     "Influence of Viscous and Buoyancy Forces on the Mobilization
!     of Residual Tetrachloroethylene during Surfactant Flushing."
!     Environ. Sci. Technol.  30(4):1328-1335.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE JACOB
      USE HYST
      USE GRID
      USE FLUX
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
      IF( ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TRPGL_W'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = ND(N)
        DFMLX = SQRT( ABS(UL(1,1,N))*ABS(UL(1,2,N)) +
     &    ABS(VL(1,1,N))*ABS(VL(1,2,N)) +
     &    ABS(WL(1,1,N))*ABS(WL(1,2,N)))
        ULX = 5.D-1*(UL(1,1,N)+UL(1,2,N))
        VLX = 5.D-1*(VL(1,1,N)+VL(1,2,N))
        WLX = 5.D-1*(WL(1,1,N)+WL(1,2,N))
        ULGX = 5.D-1*(UL(1,1,N)*GRVX(1,N)+UL(1,2,N)*GRVX(2,N))
        VLGX = 5.D-1*(VL(1,1,N)*GRVY(1,N)+VL(1,2,N)*GRVY(2,N))
        WLGX = 5.D-1*(WL(1,1,N)*GRVZ(1,N)+WL(1,2,N)*GRVZ(2,N))
        DFALX = ((ULGX + VLGX + WLGX)/GRAV)/
     &    ( SQRT( ULX**2 + VLX**2 + WLX**2 ) + SMALL )
        SKL = SQRT((PERM(1,N)*(5.D-1*(GRVX(1,N)+GRVX(2,N))))**2 +
     &    (PERM(2,N)*5.D-1*(GRVY(1,N)+GRVY(2,N)))**2 +
     &    (PERM(3,N)*5.D-1*(GRVZ(1,N)+GRVZ(2,N)))**2)/GRAV
        BNDX = (RHOL(2,N)-RHORG)*GRAV*SKL*RKL(2,N)
        CAPX = DFMLX*VISL(2,N)
        TRPGL(2,N) = SQRT(CAPX**2 + 2.D+0*CAPX*BNDX*DFALX + BNDX**2)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TRPGL_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Update the primary variables on field cells w/o ghost cells.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE OUTPU
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(4) :: VARZ	
      REAL*8, DIMENSION(3) :: VARX	
      INTEGER, DIMENSION(3) :: IVARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_W'
      OFFSET = IOFFSET_REF
      IF( ICNV.EQ.1 ) THEN
        IF( NTSR.LT.4 .OR. (DTCF*DT).GT.DTMN ) THEN
          NTSR = NTSR + 1
          DTX = DT
          TM = TM - (1.D+0-DTCF)*DT
          DT = DTCF*DT
          DTO = DT
          DTI = 1.D+0/DT
          NCH = INDEX( UNTM(1:),'  ') - 1
          IF( ID.EQ.0 ) PRINT *,'Linear System Solver Failure: ' //
     &      'Time Step  Reduced from ',DTX*CNVTM,UNTM(1:NCH),' to ',
     &      DT*CNVTM,UNTM(1:NCH)
!
!---      Write a linear system solver failure index of -3 in the 
!         NSTEP location plus time step reductions to output.bin  ---
!
          NVAR = 1
          IVARX(1) = -3
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &      MPI_INTEGER,STATUS,IERR)
          OFFSET = OFFSET + NVAR*NBYTI
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
          NVAR = 2
          VARX(1) = DTX*CNVTM
          VARX(2) = DT*CNVTM
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,VARX,NVAR,
     &      MPI_REAL8,STATUS,IERR)
          OFFSET = OFFSET + NVAR*NBYTR
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTR
!
!---      Reset principal variables to old time step values
!         on field and ghost cells  ---
!
          DO N = 1,NFCGC(ID+1)
            PL(2,N) = PL(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
          ENDDO





!
!---    Number of time step reductions failure: stop simulation  ---
!
        ELSE
!
!---      Reset principal variables to old time step values
!         on field and ghost cells  ---
!
          DO N = 1,NFCGC(ID+1)
            PL(2,N) = PL(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
          ENDDO





          IF( ID.EQ.0 )  PRINT *,'          ---  Time Step ' // 
     &      'Reduction Limit Exceeded  ---'
          ICNV = 4
          IVARX(1) = -2
!
!---      Write a time-step reduction limit exceeded index in the 
!         NSTEP location to output.bin.  ---
!
          NVAR = 1
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &      MPI_INTEGER,STATUS,IERR)
          OFFSET = OFFSET + NVAR*NBYTI
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
        ENDIF
!
!---    Reset subroutine string sequence  ---
!
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Update primary variables on field cells  ---
!
      IERR = 0
      IERRL = NFLD + 1
      IERRG = NFLD + 1
      NMD = 0
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
        IF( IERRL.NE.(NFLD+1) ) CYCLE
        MPW = NMD + IEQW
        N_DB = ND(N)
        PAE = 0.D+0
        IF( MOD(ISCHR(N),10).EQ.2 ) PAE = SCHR(1,N)*RHORL*GRAV
!
!---    Limit aqueous pressure updates to changes in aqueous
!       saturation of 0.125 for unsaturated conditions   ---
!
        SLX = SL(2,N)-SIGN(0.125D+0,(SL(2,N)-0.5D+0))
        IF( (ISCHR(N).NE.301) .AND. (ISCHR(N).NE.302) .AND.
     &    (ISCHR(N).NE.3) .AND. (ISCHR(N).NE.4) .AND.
     &    (PG(2,N)-PL(2,N)-PAE.GT.0.D+0) ) THEN
          CPGLO = PG(1,N)-PL(1,N)
          CALL CAP_W( N,SLX,SGT(2,N),CPGL,SL(1,N),CPGLO,IPH(2,N) )
          DPW = ABS( CPGL-PG(2,N)+PL(2,N) )
          DPW = MIN( DPW,ABS(BLU(MPW)) )
          DPW = SIGN( DPW,BLU(MPW) )*RLXF
        ELSE
          DPW = BLU(MPW)*RLXF
        ENDIF
!
!---    Relax aqueous pressure updates when transitioning from
!       unsaturated to saturated states   ---
!
        IF( PG(2,N)-PL(2,N)-PAE.GT.ZERO .AND.
     &      PG(2,N)-PL(2,N)-PAE-DPW.LT.ZERO ) DPW = 6.D-1*DPW
        PL(2,N) = PL(2,N) + DPW
!
!---    Increment equation counter for next active node  ---
!
        NMD = NMD + ISVC
!
!---    Check for excessive aqueous pressure   ---
!
        IF( PL(2,N).GT.PMX-PATM ) IERR = 1
!
!---    Excess changes in primary variables, skip updates  ---
!
        IF( IERR.EQ.1 .AND. IERRL.EQ.(NFLD+1) ) THEN
          IERRL = ND(N)
          VARZ(1) = REAL( N )
          VARZ(2) = REAL( NPHAZ(2,N) )
          VARZ(3) = DPW
          VARZ(4) = PL(2,N)+PATM
        ENDIF
      ENDDO
!
!---  Identify minimum global node with errors  ---
!
      CALL MPI_ALLREDUCE( IERRL,IERRG,1,MPI_INTEGER,MPI_MIN,
     &  MPI_COMM_WORLD,IERR )
!
!---  Node identified with an excessive primary variable change  ---
!
      IF( IERRG.NE.(NFLD+1) ) THEN
        ICNV = 1
        ID2X = 0
!
!---    Identify processor with error  ---
!
        DO MP = 1,NP
          IF( IERRL.EQ.IERRG ) ID2X = ID
        ENDDO
        CALL MPI_ALLREDUCE( ID2X,IDX,1,MPI_INTEGER,MPI_MAX,
     &    MPI_COMM_WORLD,IERR )
!
!---    Send node state to P0  ---
!
        IF( ID.EQ.IDX ) THEN
          CALL MPI_SEND( VARZ,29,MPI_REAL8,0,18,MPI_COMM_WORLD,IERR )
        ENDIF
        IF( ID.EQ.0 ) THEN
          CALL MPI_RECV( VARZ,4,MPI_REAL8,IDX,18,MPI_COMM_WORLD,
     &      STATUS,IERR )
          N = INT(VARZ(1))
          NPHZX = MOD(INT(VARZ(2)),100)
          PRINT *,'---  Primary Variable(s) Error  ---'
          PRINT *,'  Node = ',ND(N)   
          PRINT *,'  DPW = ',VARZ(3)
          PRINT *,'  Aqueous Pressure, Pa = ',VARZ(4)
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
          NCH = INDEX( UNTM(1:),'  ' ) - 1
          IF( ID.EQ.0 ) PRINT *,'  Time Step Reduced From ',
     &      DTX*CNVTM,UNTM(1:NCH),' to ',DT*CNVTM,UNTM(1:NCH)
          DO N = 1,NFCGC(ID+1)
            PL(2,N) = PL(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
          ENDDO





!
!---      Write a excessive primary variable index in the NSTEP location
!         plus write the global node numbers and phase condition indices
!         for the location of maximum residuals for the water, air, 
!         and salt equations to output.bin  ---
!
          NVAR = 3
          IVARX(1) = -4
          IVARX(2) = ND(IERRG)
          IVARX(3) = NPHZX
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &      MPI_INTEGER,STATUS,IERR)
          OFFSET = OFFSET + NVAR*NBYTI
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
!
!---      Write maximum residual for the water equation
!         and time step reductions to output.bin  ---
!
          NVAR = 3
          VARX(1) = VARZ(4)
          VARX(2) = DTX*CNVTM
          VARX(3) = DT*CNVTM
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,VARX,NVAR,
     &      MPI_REAL8,STATUS,IERR)
          OFFSET = OFFSET + NVAR*NBYTR
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTR
!
!---    Number of time step reductions failure: stop simulation  ---
!
        ELSE
          DO N = 1,NFCGC(ID+1)
            PL(2,N) = PL(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
          ENDDO





          IF( ID.EQ.0 ) PRINT *,'          ---  Time Step Reduction ' // 
     &        'Limit Exceeded  ---'
          ICNV = 4
!
!---      Write a time-step reduction limit exceeded index in the 
!         NSTEP location to output.bin.  ---
!
          IVARX(1) = -2
          NVAR = 1
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &      MPI_INTEGER,STATUS,IERR)
          OFFSET = OFFSET + NVAR*NBYTI
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_W group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_GC_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Update the primary variables on ghost cells
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June, 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE OUTPU
      USE JACOB
      USE HYST
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
      INTEGER	STATUS(MPI_STATUS_SIZE)	
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_GC_W'
      NPVX = 1
!
!---  Load sending buffer for bottom ghost cells for processors
!     with bottom ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(1,ID+1).GT.0 ) THEN
        MCS = 0
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(1,ID+1) = ',NCGC(1,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(1,ID+1)
          SBFB(NCS+1) = PL(2,NLSGC(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(1,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(1,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Bottom Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFB,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post Bottom Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending bottom ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(6,ID+1).GT.0 ) THEN
!        NRCVX = NCGC(6,ID+1)*NPVX
        IDSNDX = NPGC(6,ID+1) - 1
        NRCVX = NCGC(1,IDSNDX+1)*NPVX
        IDRCVX = NPGC(1,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Bottom Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFB,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post Bottom Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,5
          MCR = MCR + NCGC(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGC(6,ID+1) = ',NCGC(6,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(6,ID+1)
          PL(2,NLRGC(M+MCR)) = RBFB(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for south ghost cells for processors
!     with south ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(2,ID+1).GT.0 ) THEN
        MCS = NCGC(1,ID+1)
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(2,ID+1) = ',NCGC(2,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(2,ID+1)
          SBFS(NCS+1) = PL(2,NLSGC(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(2,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(2,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'South Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFS,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post South Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending south ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(5,ID+1).GT.0 ) THEN
!        NRCVX = NCGC(5,ID+1)*NPVX
        IDSNDX = NPGC(5,ID+1) - 1
        NRCVX = NCGC(2,IDSNDX+1)*NPVX
        IDRCVX = NPGC(2,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'South Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFS,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post South Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,4
          MCR = MCR + NCGC(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGC(5,ID+1) = ',NCGC(5,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(5,ID+1)
          PL(2,NLRGC(M+MCR)) = RBFS(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for west ghost cells for processors
!     with west ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(3,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,2
          MCS = MCS + NCGC(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(3,ID+1) = ',NCGC(3,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(3,ID+1)
!          IF( ID.EQ.2 ) 
!     &      PRINT *,'S2: ND(',NLSGC(M+MCS),') = ',ND(NLSGC(M+MCS))
          SBFW(NCS+1) = PL(2,NLSGC(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(3,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(3,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'West Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFW,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post West Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending west ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(4,ID+1).GT.0 ) THEN
!        NRCVX = NCGC(4,ID+1)*NPVX
        IDSNDX = NPGC(4,ID+1) - 1
        NRCVX = NCGC(3,IDSNDX+1)*NPVX
        IDRCVX = NPGC(3,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'West Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFW,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post West Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,3
          MCR = MCR + NCGC(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGC(4,ID+1) = ',NCGC(4,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(4,ID+1)
!          IF( ID.EQ.1 ) 
!     &      PRINT *,'R1: ND(',NLSGC(M+MCR),') = ',ND(NLSGC(M+MCR))
          PL(2,NLRGC(M+MCR)) = RBFW(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for east ghost cells for processors
!     with east ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(4,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,3
          MCS = MCS + NCGC(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(4,ID+1) = ',NCGC(4,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(4,ID+1)
          SBFE(NCS+1) = PL(2,NLSGC(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(4,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(4,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'East Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFE,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post East Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending east ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(3,ID+1).GT.0 ) THEN
!        NRCVX = NCGC(3,ID+1)*NPVX
        IDSNDX = NPGC(3,ID+1) - 1
        NRCVX = NCGC(4,IDSNDX+1)*NPVX
        IDRCVX = NPGC(4,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'East Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFE,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post East Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,2
          MCR = MCR + NCGC(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGC(3,ID+1) = ',NCGC(3,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(3,ID+1)
          PL(2,NLRGC(M+MCR)) = RBFE(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for north ghost cells for processors
!     with north ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(5,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,4
          MCS = MCS + NCGC(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(5,ID+1) = ',NCGC(5,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(5,ID+1)
          SBFN(NCS+1) = PL(2,NLSGC(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(5,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(5,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'North Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFN,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post North Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending north ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(2,ID+1).GT.0 ) THEN
!        NRCVX = NCGC(2,ID+1)*NPVX
        IDSNDX = NPGC(2,ID+1) - 1
        NRCVX = NCGC(5,IDSNDX+1)*NPVX
        IDRCVX = NPGC(5,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'North Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFN,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post North Receive: IERR = ',IERR,' ID = ',ID
        MCR = NCGC(1,ID+1)
!        PRINT *,' MCR = ',MCR,' NCGC(2,ID+1) = ',NCGC(2,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(2,ID+1)
          PL(2,NLRGC(M+MCR)) = RBFN(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for top ghost cells for processors
!     with top ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(6,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,5
          MCS = MCS + NCGC(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(6,ID+1) = ',NCGC(6,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(6,ID+1)
          SBFT(NCS+1) = PL(2,NLSGC(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(6,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(6,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Top Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFT,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post Top Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending top ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(1,ID+1).GT.0 ) THEN
!        NRCVX = NCGC(1,ID+1)*NPVX
        IDSNDX = NPGC(1,ID+1) - 1
        NRCVX = NCGC(6,IDSNDX+1)*NPVX
        IDRCVX = NPGC(6,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Top Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFT,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post Top Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
!        PRINT *,' MCR = ',MCR,' NCGC(1,ID+1) = ',NCGC(1,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(1,ID+1)
          PL(2,NLRGC(M+MCR)) = RBFT(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_GC_W group
!
      RETURN
      END
      



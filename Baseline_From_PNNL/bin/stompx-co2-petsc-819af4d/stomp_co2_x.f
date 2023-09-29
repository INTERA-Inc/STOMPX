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
!     STOMP-CO2
!
!     This engineering program numerically simulates the transport
!     of H2O, NaCl and CO2 through multifluid subsurface environments
!     under isothermal conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 September 2021
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
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSCKSP
      USE PETSC_STOMP

!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)





!
!
!  Include file for Fortran use of the Mat package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!




!
!
!  Include file for Fortran use of the Vec package in PETSc
!




!
!
!  Include file for Fortran use of the AO (application ordering) package in PETSc
!




!
!
!  Include file for Fortran use of the IS (index set) package in PETSc
!




!
!
!  Part of the base include file for Fortran use of PETSc.
!  Note: This file should contain only define statements and
!  not the declaration of variables.

! No spaces for #defines as some compilers (PGI) also adds
! those additional spaces during preprocessing - bad for fixed format
!
























































!
!  Include file for Fortran use of the PetscViewer package in PETSc
!













!
!  No includes needed for logging

!
!  Include file for Fortran use of the Bag package in PETSc
!







!
! The real*8,complex*16 notatiton is used so that the
! PETSc double/complex variables are not affected by
! compiler options like -r4,-r8, sometimes invoked
! by the user. NAG compiler does not like integer*4,real*8
























!
! Fortran does not support unsigned, though ISO_C_BINDING
! supports INTEGER(KIND=C_SIZE_T). We don't use that here
! only to avoid importing the module.





!

!






!

!


!



!
!     Macro for templating between real and complex
!










!
!    Allows the matrix Fortran Kernels to work with single precision
!    matrix data structures
!

!
!     PetscLogDouble variables are used to contain double precision numbers
!     that are not used in the numerical computations, but rather in logging,
!     timing etc.
!


!
!     Macros for error checking
!




















!
!  Include file for Fortran use of the type(tPetscViewer) package in PETSc
!




































































!
!  Matrix types
!


!
! MatMFFDType values
!



!
! MatSolverTypes
!


!
! GPU Storage Formats for CUSPARSE
!



!
! GPU Storage Formats for HIPSPARSE
!



!
! sparsity reducing ordering for STRUMPACK
!



!
!
!  Include file for Fortran use of the type(tVec) package in PETSc
!

!
!
!  Include file for Fortran use of the KSP package in PETSc
!




!
!
!  Include file for Fortran use of the PC (preconditioner) package in PETSc
!




!
!
!  Include file for Fortran use of the type(tMat) package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!


!
!  Include file for Fortran use of the DM package in PETSc
!




!
!
!  Include file for Fortran use of the type(tIS) (index set) package in PETSc
!

!
!
!  Include file for Fortran use of the type(tVec) package in PETSc
!

!
!
!  Include file for Fortran use of the type(tMat) package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!



















!
! GAMG types
!



!
! GAMG classical types
!



!
! Various preconditioners
!









!
!  Various Krylov subspace methods
!

!
!  Various Initial guesses for Krylov subspace methods
!


!
!----------------------Type Declarations-------------------------------!
!
      LOGICAL HALT,PLOT,RESTART




      integer(kind=selected_int_kind(5)) :: IERR

      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      INTEGER, DIMENSION(:), ALLOCATABLE :: ID_NNZ,IO_NNZ
      CHARACTER*132 CHMSGX(2)
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = 1
      SUB_LOG(1) = 'STOMPX-CO2'
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
!---  Initialize PETSc ---
!
      CALL PetscInitialize(PETSC_NULL_CHARACTER,IERR )

!
!---  Print banner to screen ---
!
      IF( ID.EQ.0 ) THEN
      PRINT *,' Welcome to ...'
      PRINT *,' '
      PRINT *,'                     STOMPX-CO2'
      PRINT *,'        Subsurface Transport Over Multiple Phases'
      PRINT *,'                 CO2 Operational Mode'
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
      CALL READ_BIN_CO2
!      ISKIP = 1
!      IF( ISKIP.EQ.0 ) THEN
!
!---  Initialize global arrays  ---
!
      CALL INTLZ_CO2
!      PRINT *,'Post INTLZ_CO2: ID = ',ID
!
!---  For geomechanics set k iterate value of pore pressure  ---
!
      IF( ISLC(50).NE.0 ) THEN
        INDX = 2
        CALL PRESS_GM( INDX )
!        PRINT *,'Post PRESS_GM: ID = ',ID
      ENDIF
!
!---  Configure the arrays for compressed sparse row matrix storage  ---
!
      CALL JCBP_CO2
!      PRINT *,'Post JCBP_CO2: ID = ',ID
!
!---  For geomechanics simulations compute Jacobian matrix pointers  --
!
      IF( ISLC(50).NE.0 ) CALL JCBP_GM
!      PRINT *,'Post JCBP_GM: ID = ',ID
!
!---  Compute primary variable increments  ---
!
      CALL INCRM_CO2
!      PRINT *,'Post INCRM_CO2: ID = ',ID
!
!---  Hydrologic and thermodynamic properties  ---
!
      CALL PROP_CO2
!      PRINT *,'Post PROP_CO2: ID = ',ID
      CALL BCP_CO2
!      PRINT *,'Post BCP_CO2: ID = ',ID
!
!---  Compute initial solute concentrations  ---
!
      CALL CISC_CO2
!      PRINT *,'Post CISC_CO2: ID = ',ID
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
!      ISKIP = 1
!      IF( ISKIP.EQ.0 ) THEN
!
!---  Compute initial fluxes on non-boundary and boundary surfaces  ---
!
      ISVF = 1
      CALL FLUX_CO2
!      PRINT *,'Post FLUX_CO2: ID = ',ID
      ISVF = 2*ISVC + 1
      CALL BCF_CO2
!      PRINT *,'Post BCF_CO2: ID = ',ID
!
!---  Surface flux integrator for zero time step  ---
!
      DTX = DT
      DT = 0.D+0
      CALL SFIN_CO2
!      PRINT *,'Post SFIN_CO2: ID = ',ID
      DT = DTX

!
!---  Create PETSc matrix, solver, and solution and problem vectors
!     for coupled flow  ---
!
!      PRINT *,'NUKFL(ID+1) = ',NUKFL(ID+1),' ID = ',ID
!      PRINT *,'NUKFO(ID+1) = ',NUKFO(ID+1),' ID = ',ID
!      PRINT *,'NUKFG = ',NUKFG,' ID = ',ID
      INDX = 0
      CALL STOMP_PETSC_CREATE( F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC,
     &  F_SOL_VEC_S,F_SCATTER,INDX )
!
!---    Create Lis matrix, solver, and solution and problem vectors
!       for solute transport  ---
!
      IF( IEQC.NE.0 .OR. ISLC(40).NE.0 ) THEN
!        PRINT *,'NUKTL(ID+1) = ',NUKTL(ID+1),' ID = ',ID
!        PRINT *,'NUKTO(ID+1) = ',NUKTO(ID+1),' ID = ',ID
!        PRINT *,'NUKTG = ',NUKTG,' ID = ',ID
        INDX = 1
        CALL STOMP_PETSC_CREATE( T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC,
     &    T_SOL_VEC_S,T_SCATTER,INDX )
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
        CALL STOMP_PETSC_CREATE( G_KSP,G_MAT,G_RHS_VEC,G_SOL_VEC,
     &    G_SOL_VEC_S,G_SCATTER,INDX )
      ENDIF

!      PRINT *,'Post STOMP_PETSC_CREATE: ID = ',ID
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
!      ISKIP = 1
!      IF( ISKIP.EQ.0 ) THEN
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
          CALL REFNOD_CO2
          IF( ID.EQ.0 ) THEN
            IF( ICNO.EQ.10 ) THEN
              ICNO = 0
              NCH = INDEX(UNTM(1:),'  ') - 1
              IF( ISLC(50).NE.0 ) THEN
                PRINT *,'       Step         Itr  MLp           Time' // 
     &            ' [',UNTM(1:NCH),']             Timestep [',
     &            UNTM(1:NCH),']'
              ELSE
              PRINT *,'       Step         Itr           Time' //
     &          ' [',UNTM(1:NCH),']             Timestep [',
     &          UNTM(1:NCH),']'
            ENDIF
            ENDIF
            ICNO = ICNO + 1
            IF( ISLC(50).NE.0 ) THEN
              PRINT *,NSTEP,K_GM(2),K_GM(1),TM*CNVTM,DT*CNVTM
            ELSE
              PRINT *,NSTEP,NITER,TM*CNVTM,DT*CNVTM
            ENDIF
          ENDIF
        ENDIF
!
!---    Normalize mineral species concentrations after initial
!       output for normal simulations  ---
!
        IF( ISLC(40).EQ.1 .AND. (NSTEP-NRST).EQ.0 ) CALL NMNSP
!
!---    Update porosity and permeability in response to geomechanical
!       stress  ---
!
        IF( ISLC(50).NE.0 .AND. NSTEP.EQ.0 ) THEN
          CALL PORSTY_GM
          CALL PERMRF_GM
        ENDIF
!
!---    Load old time step arrays  ---
!
        CALL LDO_CO2
!        PRINT *,'Post LDO_CO2: ID = ',ID
!
!---    Load old time step arrays for the coupled-well model  ---
!
        IF( L_CW.EQ.1 ) THEN
          CALL LDO_COUP_WELL
!          PRINT *,'Post LDO_COUP_WELL: ID = ',ID
        ENDIF
!
!---    Load old time step arrays for the volumetric stress
!       and pore pressure  ---
!
        IF( ISLC(50).NE.0 ) THEN
          INDX = 1
          CALL LD_GM( INDX )
        ENDIF
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
          CALL WRPLOT_CO2
          IF( ISLC(18).LT.1 ) CALL WRRST_CO2
        ENDIF
!
!---    Compute the next time step and increment time step counter  ---
!
        DTSO = DT
!        PRINT *,'Pre TMSTEP: ID = ',ID
        CALL TMSTEP
!        PRINT *,'Post TMSTEP: ID = ',ID
        IF( NSTEP.EQ.0 ) DTSO = DT
        NSTEP = NSTEP + 1
        IF( NSTEP-NRST.GT.MXSTEP ) THEN
          IF( ID.EQ.0 ) PRINT *,'Simulation Stopped:  Time Step Limit'
          NSTEP = NSTEP - 1
          EXIT TSLOOP
        ENDIF
!
!---    Reset the time step reduction counter  ---
!
        NTSR = 0
!
!---    Top of sequential flow and transport and geomechanics  ---
!
        K_GM(1) = 0
        K_GM(2) = 0
        GMLOOP: DO
          K_GM(1) = K_GM(1) + 1
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
          CALL BCP_CO2
!          PRINT *,'Post BCP_CO2: ID = ',ID
!
!---      Compute coupled-well fluxes  ---
!
          IF( L_CW.EQ.1 ) CALL FLUX_COUP_WELL
!          PRINT *,'Post FLUX_COUP_WELL: ID = ',ID
!
!---      Compute source contributions  ---
!
          CALL SORC_CO2
!          PRINT *,'Post SORC_CO2: ID = ',ID
!
!---      Compute fluxes on non-boundary and boundary surfaces  ---
!
          CALL FLUX_CO2
!          PRINT *,'Post FLUX_CO2: ID = ',ID
          CALL BCF_CO2
!          PRINT *,'Post BCF_CO2: ID = ',ID
!
!---      Zero Jacobian matrix  ---
!
          CALL JCBZ_CO2
!          PRINT *,'Post JCBZ_CO2:ID = ',ID
!
!---      Load Jacobian matrix for the water, CO2
!         and salt mass equations (zero flux boundary)  ---
!
          CALL JCB_CO2
!          PRINT *,'Post JCB_CO2:ID = ',ID
!
!---      Modify the Jacobian matrix for boundary conditions  ---
!
          CALL BCJ_CO2
!          PRINT *,'Post BCJ_CO2:ID = ',ID
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
!          PRINT *,'Post JCB_SV:ID = ',ID
!
!---      Check for fatal execution errors and stop simulation
!         if detected  ---
!
          CALL CHK_ERROR
!
!---      Solve the linear system A x = b for coupled flow  ---
!








          INDX = 0
          CALL STOMP_PETSC_SOLVE( F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC,
     &      F_SOL_VEC_S,F_SCATTER,INDX )
!          PRINT *,'Post STOMP_PETSC_SOLVE:ID = ',ID
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
!            PRINT *,'Post UPDT_COUP_WELL:ID = ',ID
          ENDIF
!
!---      Update primary variables on field nodes w/o ghost cells  ---
!
          CALL UPDT_CO2
!          PRINT *,'Post UPDT_CO2:ID = ',ID







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
          CALL RSDL_CO2
!          PRINT *,'Post RSDL_CO2:ID = ',ID
!
!---      Compute primary variable increments, saturation,
!         relative permeability, porosity, tortuosity,
!         thermodynamic properties for interior nodes,
!         except immediately after a new time step  ---
!
          CALL INCRM_CO2
!          PRINT *,'Post INCRM_CO2:ID = ',ID
          CALL PROP_CO2
!          PRINT *,'Post PROP_CO2:ID = ',ID
!
!---      For geomechanics simulations alter permeability with
!         porosity  --
!
          IF( ISLC(50).NE.0 ) CALL PERMRF_GM
!          PRINT *,'Post PERMRF_GM:ID = ',ID
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
!---    Solve geomechanics  ---
!
        IF( ISLC(50).NE.0 ) THEN
!
!---      Set k+1 iterate value of pore pressure  ---
!
          INDX = 3
          CALL PRESS_GM( INDX )
!
!---      Static porothermoelastic geomechanics  ---
!
          CALL STATIC_GM
!
!---      Convergence check for sequential coupled flow and transport
!         and geomechanics  ---
!
          CALL RSDL_GM
          IF( RSD_GM.GT.RSDM_GM(IEPD) ) THEN
!
!---        Load k level arrays for the volumetric stress
!           and pore pressure  ---
!
            INDX = 2
            CALL LD_GM( INDX )
!
!---        Update porosity and permeability for geomechical stress  ---
!
            CALL PORSTY_GM
            CALL PERMRF_GM
            CYCLE GMLOOP
          ENDIF
!
!---      Update porosity and permeability for geomechical stress  ---
!
          CALL PORSTY_GM
          CALL PERMRF_GM
          EXIT GMLOOP
        ELSE
          EXIT GMLOOP
        ENDIF
        ENDDO GMLOOP
!
!---    Integrate coupled-equation sources ---
!
        CALL SORIC_CO2
!
!---    Compute current fluxes for transport solutions or flux
!       integrations  ---
!
        ISVF = 1
        CALL FLUX_CO2
        CALL BCF_CO2
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
              CALL SPRP_CO2( NSL )
!
!---          Solute transport ---
!
              CALL TPORT_CO2( NSL )
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
!---                Skip transport for linked aqueous CO2   ---
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
                      CALL TPORT_CO2( NSL )
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
!---                Skip transport for linked aqueous CO2   ---
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
                      CALL TPORT_CO2( NSL )
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
        CALL SFIN_CO2
!        PRINT *,'Post SFIN_CO2: ID = ',ID
!
!---  Proceed to new time step  ---
!
      ENDDO TSLOOP

!
!---  Destroy solution vector, problem vector, matrix, solver, and
!     solution-vector scatter scheme for the PETSc solver for coupled
!     flow and heat transport  ---
!
      INDX = 0
      CALL STOMP_PETSC_DESTROY( F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC,
     &  F_SOL_VEC_S,F_SCATTER,INDX )
!
!---  Destroy solution vector, problem vector, matrix, solver, and
!     solution-vector scatter scheme for the PETSc solver for 
!     solute transport  ---
!
      IF( IEQC.NE.0 .OR. ISLC(40).NE.0 ) THEN
        INDX = 1
        CALL STOMP_PETSC_DESTROY( T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC,
     &   T_SOL_VEC_S,T_SCATTER,INDX )
      ENDIF
!
!---  Destroy solution vector, problem vector, matrix, solver, and
!     solution-vector scatter scheme for the PETSc solver for 
!     geomechanics  ---
!
      IF( ISLC(50).NE.0 ) THEN
        INDX = 2
        CALL STOMP_PETSC_DESTROY( G_KSP,G_MAT,G_RHS_VEC,G_SOL_VEC,
     &   G_SOL_VEC_S,G_SCATTER,INDX )
      ENDIF
!
!---  Finalize PETSc execution  ---
!
      CALL PetscFinalize( IERR )

!
!---  Write plot_xxx.bin file  ---
!
      CALL WRPLOT_CO2
!
!---  Write restart_xxx.bin file  ---
!
      CALL WRRST_CO2
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
      CALL MPI_FINALIZE( IERR )
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
!     STOMPX-CO2
!
!     Compute boundary surface fluxes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021
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
      SUB_LOG(ISUB_LOG) = '/BCF_CO2'
      DO NB = 1,NBC(ID+1)
        N = IBCN(NB)
        IF( IBCD(NB).EQ.-3 ) THEN
          DO M = 1,ISVF
            WL(M,1,N) = 0.D+0
            WG(M,1,N) = 0.D+0
            WLA(M,1,N) = 0.D+0
            WLW(M,1,N) = 0.D+0
            WGA(M,1,N) = 0.D+0
            WGW(M,1,N) = 0.D+0
            WDLA(M,1,N) = 0.D+0
            WDGW(M,1,N) = 0.D+0
            WS(M,1,N) = 0.D+0
            WDS(M,1,N) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          DO M = 1,ISVF
            VL(M,1,N) = 0.D+0
            VG(M,1,N) = 0.D+0
            VLA(M,1,N) = 0.D+0
            VLW(M,1,N) = 0.D+0
            VGA(M,1,N) = 0.D+0
            VGW(M,1,N) = 0.D+0
            VDLA(M,1,N) = 0.D+0
            VDGW(M,1,N) = 0.D+0
            VS(M,1,N) = 0.D+0
            VDS(M,1,N) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          DO M = 1,ISVF
            UL(M,1,N) = 0.D+0
            UG(M,1,N) = 0.D+0
            ULA(M,1,N) = 0.D+0
            ULW(M,1,N) = 0.D+0
            UGA(M,1,N) = 0.D+0
            UGW(M,1,N) = 0.D+0
            UDLA(M,1,N) = 0.D+0
            UDGW(M,1,N) = 0.D+0
            US(M,1,N) = 0.D+0
            UDS(M,1,N) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          DO M = 1,ISVF
            UL(M,2,N) = 0.D+0
            UG(M,2,N) = 0.D+0
            ULA(M,2,N) = 0.D+0
            ULW(M,2,N) = 0.D+0
            UGA(M,2,N) = 0.D+0
            UGW(M,2,N) = 0.D+0
            UDLA(M,2,N) = 0.D+0
            UDGW(M,2,N) = 0.D+0
            US(M,2,N) = 0.D+0
            UDS(M,2,N) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          DO M = 1,ISVF
            VL(M,2,N) = 0.D+0
            VG(M,2,N) = 0.D+0
            VLA(M,2,N) = 0.D+0
            VLW(M,2,N) = 0.D+0
            VGA(M,2,N) = 0.D+0
            VGW(M,2,N) = 0.D+0
            VDLA(M,2,N) = 0.D+0
            VDGW(M,2,N) = 0.D+0
            VS(M,2,N) = 0.D+0
            VDS(M,2,N) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          DO M = 1,ISVF
            WL(M,2,N) = 0.D+0
            WG(M,2,N) = 0.D+0
            WLA(M,2,N) = 0.D+0
            WLW(M,2,N) = 0.D+0
            WGA(M,2,N) = 0.D+0
            WGW(M,2,N) = 0.D+0
            WDLA(M,2,N) = 0.D+0
            WDGW(M,2,N) = 0.D+0
            WS(M,2,N) = 0.D+0
            WDS(M,2,N) = 0.D+0
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
              WL(M,1,N) = BCX(3)
            ENDDO
!
!---        Isobrine option and salt boundary conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSB( N,NB )
            ENDIF
            CALL DFFLAWB( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLB( N,NB )
!
!---        Isobrine option and salt boundary conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSB( N,NB )
            ENDIF
            CALL DFFLAWB( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              WG(M,1,N) = BCX(4)
            ENDDO
            CALL DFFGAWB( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGB( N,NB )
            CALL DFFGAWB( N,NB )
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
              VL(M,1,N) = BCX(3)
            ENDDO
!
!---        Isobrine option and Salt Boundary Conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSS( N,NB )
            ENDIF
            CALL DFFLAWS( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLS( N,NB )
!
!---        Isobrine option and Salt Boundary Conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSS( N,NB )
            ENDIF
            CALL DFFLAWS( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              VG(M,1,N) = BCX(4)
            ENDDO
            CALL DFFGAWS( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGS( N,NB )
            CALL DFFGAWS( N,NB )
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
              UL(M,1,N) = BCX(3)
            ENDDO
!
!---        Isobrine option and Salt Boundary Conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSW( N,NB )
            ENDIF
            CALL DFFLAWW( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLW( N,NB )
!
!---        Isobrine option and Salt Boundary Conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSW( N,NB )
            ENDIF
            CALL DFFLAWW( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              UG(M,1,N) = BCX(4)
            ENDDO
            CALL DFFGAWW( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGW( N,NB )
            CALL DFFGAWW( N,NB )
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              UL(M,2,N) = BCX(3)
            ENDDO
!
!---        Isobrine option and Salt Boundary Conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSE( N,NB )
            ENDIF
            CALL DFFLAWE( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLE( N,NB )
!
!---        Isobrine option and Salt Boundary Conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSE( N,NB )
            ENDIF
            CALL DFFLAWE( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              UG(M,2,N) = BCX(4)
            ENDDO
            CALL DFFGAWE( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGE( N,NB )
            CALL DFFGAWE( N,NB )
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              VL(M,2,N) = BCX(3)
            ENDDO
!
!---        Isobrine option and Salt Boundary Conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSN( N,NB )
            ENDIF
            CALL DFFLAWN( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLN( N,NB )
!
!---        Isobrine option and Salt Boundary Conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLSN( N,NB )
            ENDIF
            CALL DFFLAWN( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              VG(M,2,N) = BCX(4)
            ENDDO
            CALL DFFGAWN( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGN( N,NB )
            CALL DFFGAWN( N,NB )
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Aqueous Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQW,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              WL(M,2,N) = BCX(3)
            ENDDO
!
!---        Isobrine option and Salt Boundary Conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLST( N,NB )
            ENDIF
            CALL DFFLAWT( N,NB )
          ELSEIF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL DRCVLT( N,NB )
!
!---        Isobrine option and Salt Boundary Conditions  ---
!
            IF( ISLC(32).EQ.0 .AND. IBCT(IEQS,NB).NE.3 ) THEN
              CALL DFFLST( N,NB )
            ENDIF
            CALL DFFLAWT( N,NB )
          ENDIF
!
!---      Gas Neumann else Dirichlet, Saturated, Unit Gradient  ---
!
          IF( IBCT(IEQA,NB).EQ.2 ) THEN
            DO M = 1,ISVF
              WG(M,2,N) = BCX(4)
            ENDDO
            CALL DFFGAWT( N,NB )
          ELSEIF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL DRCVGT( N,NB )
            CALL DFFGAWT( N,NB )
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     STOMPX-CO2
!
!     Modify the Jacobian matrix for boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 17 August 2021
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
      SUB_LOG(ISUB_LOG) = '/BCJ_CO2'
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
            CALL JCBLWB_CO2( N,NB )
            CALL JCBLAB_CO2( N,NB )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWB_CO2( N,NB )
            CALL JCBGAB_CO2( N,NB )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSB_CO2( N,NB )
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWS_CO2( N,NB )
            CALL JCBLAS_CO2( N,NB )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWS_CO2( N,NB )
            CALL JCBGAS_CO2( N,NB )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSS_CO2( N,NB )
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWW_CO2( N,NB )
            CALL JCBLAW_CO2( N,NB )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWW_CO2( N,NB )
            CALL JCBGAW_CO2( N,NB )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSW_CO2( N,NB )
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWE_CO2( N,NB )
            CALL JCBLAE_CO2( N,NB )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWE_CO2( N,NB )
            CALL JCBGAE_CO2( N,NB )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSE_CO2( N,NB )
          ENDIF
!
!---    North boundary
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWN_CO2( N,NB )
            CALL JCBLAN_CO2( N,NB )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWN_CO2( N,NB )
            CALL JCBGAN_CO2( N,NB )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBSN_CO2( N,NB )
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(IEQW,NB).NE.3 ) THEN
            CALL JCBLWT_CO2( N,NB )
            CALL JCBLAT_CO2( N,NB )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(IEQA,NB).NE.3 ) THEN
            CALL JCBGWT_CO2( N,NB )
            CALL JCBGAT_CO2( N,NB )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(IEQS,NB).NE.3 ) CALL JCBST_CO2( N,NB )
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCJ_CO2 group
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
!     STOMPX-CO2
!
!     Compute saturation, relative permeability and thermodynamic
!     properties for boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 01 November 2021
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
        ENDDO
      ENDIF
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
             IFIND = 1
             EXIT
            ENDIF
          ENDDO
          IF( IFIND.EQ.0 ) CYCLE
        ENDIF
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
        IBD = ABS(IBCD(NB))
        N_DB = -NB
!
!---    Assign gas-entry pressure for non Brooks-Corey;
!       Brooks-Corey; Brooks-Corey, Dual Porosity; and
!       Brooks-Corey, Entrapment  ---
!
        IF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 .OR.
     &      ISCHR(N).EQ.202 ) THEN
          ENPR = SCHR(1,N)*RHORL*GRAV
        ELSEIF( ISCHR(N).EQ.4 ) THEN
          ENPR = MIN( SCHR(1,N),SCHR(5,N) )*RHORL*GRAV
        ELSE
          ENPR = 0.D+0
        ENDIF
!
!---    Initial trapped gas saturation for the van Genuchten or
!       Brooks/Corey entrapment model  ---
!
        IF( ISCHR(N).EQ.101 .OR. ISCHR(N).EQ.102 .OR.
     &      ISCHR(N).EQ.201 .OR. ISCHR(N).EQ.202 ) THEN
!
!---      w/ Webb extension  ---
!
          IF( ISM(N).EQ.2 ) THEN
            ESGTMX = SCHR(15,N)
!
!---      w/o Webb extension  ---
!
          ELSE
            ESGTMX = SCHR(15,N)/(1.D+0-SCHR(4,N))
          ENDIF
        ELSE
          ESGTMX = 0.D+0
        ENDIF
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
!---    Loop over secondary variable indices  ---
!
        DO M = 2,ISVC+2
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
            PGX = PGX + BCX(4)*DB*VISG(M,N)/PERM(IBD,N)
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
            PLX = PLX + BCX(3)*DB*VISL(M,N)/PERM(IBD,N)
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
            IF( IBCM(NB).EQ.1 .AND. (NSTEP-NRST).GT.1 ) EXIT
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
            CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),N )
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
              DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSX(2)))
              DO NC = 1,32
                IF( YLSX(2).GE.XLSMX ) THEN
                  YLSB(M,NB) = XLSMX
                  EXIT
                ENDIF
!
!---            Single-variable Newton-Raphson loop (YLS)  ---
!
                DO L = 2,3
                  YLSX(L) = YLSX(2)
                  IF( L.EQ.3 ) YLSX(L) = YLSX(2) + DYLSX
                  CALL SP_B( TX,YLSX(L),PSBX )
                  PX = MAX( PGX,PSBX )
                  CALL DENS_B( TX,PX,YLSX(L),RHOBX )
                  CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSX(L),N )
                  PVWB(M,NB) = BCX(5)*PVBX
                  XLSSX = YLSX(L)
                  CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAX,XGWX,
     &              XLASX,XLSSX,XLWSX,XMGAX,XMGWX,XMLAX,XMLSX,XMLWX )
                  IF( BCX(7).LT.0 ) THEN
                    XLAX = MIN(-BCX(7),XLASX) ! mass fraction
                  ELSE
                    XLAX = BCX(7)*XLASX       ! relative saturation
                  ENDIF
                  XLSX = YLSX(L) + (XLSSX-YLSX(L))*(XLAX/XLASX)
                  CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
                  GX(L) = XLSX - RHOLSX/RHOLX
                ENDDO
                RX(1) = (GX(3)-GX(2))/DYLSX
                RPX(1) = -GX(2)
                CYLSX = RPX(1)/RX(1)
                YLSX(2) = MAX( YLSX(2)+CYLSX,0.D+0 )
                IF( ABS(CYLSX).LE.(1.D-9*XLSMX) ) THEN
                  YLSB(M,NB) = YLSX(2)
                  EXIT
                ENDIF
              ENDDO
              IF( NC.GT.32 ) THEN
                M_ERR(1) = 'Unconverged Boundary Condition: '
     &            // 'Salt Aqu. Concentration = '
                M_ERR(2) = ' at Boundary Surface: '
                CALL PATH
                R_ERR = BCX(6)
                I_ERR(1) = NB
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
                YLSB(M,NB) = RHOLSX/(RHOLWX+RHOLSX)
              ENDIF
              CALL SP_B( TX,YLSB(M,NB),PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( TX,PX,YLSB(M,NB),RHOBX )
              CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),N )
              PVWB(M,NB) = BCX(5)*PVBX
              PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
              XLSSX = YLSB(M,NB)
              CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAB(M,NB),
     &          XGWB(M,NB),XLASX,XLSSX,XLWSX,XMGAB(M,NB),XMGWB(M,NB),
     &          XMLAX,XMLSX,XMLWX )
              IF( BCX(7).LT.0 ) THEN
                XLAB(M,NB) = MIN(-BCX(7),XLASX) ! mass fraction
              ELSE
                XLAB(M,NB) = BCX(7)*XLASX       ! relative saturation
              ENDIF
              XLSB(M,NB) = YLSB(M,NB) +
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
              CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),N )
              PVWB(M,NB) = BCX(5)*PVBX
              PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
              XLSSX = YLSB(M,NB)
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
              DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSX(2)))
              DO NC = 1,32
                IF( YLSX(2).GE.XLSMX ) THEN
                  YLSB(M,NB) = XLSMX
                  EXIT
                ENDIF
!
!---            Single-variable Newton-Raphson loop (YLS)  ---
!
                DO L = 2,3
                  YLSX(L) = YLSX(2)
                  IF( L.EQ.3 ) YLSX(L) = YLSX(2) + DYLSX
                  CALL SP_B( TX,YLSX(L),PSBX )
                  PX = MAX( PGX,PSBX )
                  CALL DENS_B( TX,PX,YLSX(L),RHOBX )
                  CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSX(L),N )
                  PVWB(M,NB) = BCX(5)*PVBX
                  XLSSX = YLSX(L)
                  CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAX,XGWX,
     &              XLASX,XLSSX,XLWSX,XMGAX,XMGWX,XMLAX,XMLSX,XMLWX )
                  IF( BCX(7).LT.0 ) THEN
                    XLAX = MIN(-BCX(7),XLASX) ! mass fraction
                  ELSE
                    XLAX = BCX(7)*XLASX       ! relative saturation
                  ENDIF
                  XLSX = YLSX(L) + (XLSSX-YLSX(L))*(XLAX/XLASX)
                  GX(L) = BCX(6) - XLSX
                ENDDO
                RX(1) = (GX(3)-GX(2))/DYLSX
                RPX(1) = -GX(2)
                CYLSX = RPX(1)/RX(1)
                YLSX(2) = MAX( YLSX(2)+CYLSX,0.D+0 )
                IF( ABS(CYLSX).LE.(1.D-9*XLSMX) ) THEN
                  YLSB(M,NB) = YLSX(2)
                  EXIT
                ENDIF
              ENDDO
              IF( NC.GT.32 ) THEN
                M_ERR(1) = 'Unconverged Boundary Condition: '
     &            // 'Salt Aqu. Mass Fraction = '
                M_ERR(2) = ' at Boundary Surface: '
                CALL PATH
                R_ERR = BCX(6)
                I_ERR(1) = NB
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
                YLSB(M,NB) = BCX(6)
              ENDIF
              CALL SP_B( TX,YLSB(M,NB),PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( TX,PX,YLSB(M,NB),RHOBX )
              CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),N )
              PVWB(M,NB) = BCX(5)*PVBX
              PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
              XLSSX = YLSB(M,NB)
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
              DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSX(2)))
              DO NC = 1,32
                IF( YLSX(2).GE.XLSMX ) THEN
                  YLSB(M,NB) = XLSMX
                  EXIT
                ENDIF
!
!---            Single-variable Newton-Raphson loop (YLS)  ---
!
                DO L = 2,3
                  YLSX(L) = YLSX(2)
                  IF( L.EQ.3 ) YLSX(L) = YLSX(2) + DYLSX
                  CALL SP_B( TX,YLSX(L),PSBX )
                  PX = MAX( PGX,PSBX )
                  CALL DENS_B( TX,PX,YLSX(L),RHOBX )
                  CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSX(L),N )
                  PVWB(M,NB) = BCX(5)*PVBX
                  XLSSX = YLSX(L)
                  CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVWB(M,NB),XGAX,XGWX,
     &              XLASX,XLSSX,XLWSX,XMGAX,XMGWX,XMLAX,XMLSX,XMLWX )
                  IF( BCX(7).LT.0 ) THEN
                    XLAX = MIN(-BCX(7),XLASX) ! mass fraction
                  ELSE
                    XLAX = BCX(7)*XLASX       ! relative saturation
                  ENDIF
                  XLSX = YLSX(L) + (XLSSX-YLSX(L))*(XLAX/XLASX)
                  GX(L) = (BCX(6)*WTMS/(1.D+0 + BCX(6)*WTMS)) - XLSX
                ENDDO
                RX(1) = (GX(3)-GX(2))/DYLSX
                RPX(1) = -GX(2)
                CYLSX = RPX(1)/RX(1)
                YLSX(2) = MAX( YLSX(2)+CYLSX,0.D+0 )
                IF( ABS(CYLSX).LE.(1.D-9*XLSMX) ) THEN
                  YLSB(M,NB) = YLSX(2)
                  EXIT
                ENDIF
              ENDDO
              IF( NC.GT.32 ) THEN
                M_ERR(1) = 'Unconverged Boundary Condition: '
     &            // 'Salt Aqu. Molality = '
                M_ERR(2) = ' at Boundary Surface: '
                CALL PATH
                R_ERR = BCX(6)
                I_ERR(1) = NB
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
                YLSB(M,NB) = BCX(6)*WTMS/(1.D+0 + BCX(6)*WTMS)
              ENDIF
              CALL SP_B( TX,YLSB(M,NB),PSBX )
              PX = MAX( PGX,PSBX )
              CALL DENS_B( TX,PX,YLSB(M,NB),RHOBX )
              CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),N )
              PVWB(M,NB) = BCX(5)*PVBX
              PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
              XLSSX = YLSB(M,NB)
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
              CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),N )
              PVWB(M,NB) = BCX(5)*PVBX
              PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
              XLSB(M,NB) = YLSB(M,NB)
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
            CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX,YLSB(M,NB),N )
            PVWB(M,NB) = BCX(5)*PVBX
            PVAB(M,NB) = MAX( PGX-PVBX,0.D+0 )
            XLSSX = YLSB(M,NB)
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
!
!---      Porous-media porosity  ---
!
          CALL PORSTY_CO2E( PX,PCMP(N),PORDB(M,NB),PORTB(M,NB),N )
          PORDB(M,NB) = MAX( PORDB(M,NB),EPSL )
          PORTB(M,NB) = MAX( PORTB(M,NB),PORDB(M,NB) )
!
!---      Surface tension and saturation  ---
!
          CALL SFT_L( TX,XLSB(M,NB),SFTLX )
          BTGLB(M,NB) = 1.D+0
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &      BTGLB(M,NB) = SCHR(16,N)/SFTLX
          INDX = 0
!
!---      No trapped gas on the boundary surface  ---
!
          ASLMINX = 1.D+0
          CALL KSP_CO2( N,PGX,PLX,SGB(M,NB),SGTX,SLB(M,NB),SLFX,SLMX,
     &      RKLB(M,NB),RKGB(M,NB),ASLX,ASLMINX,ESGTX,ESGTMX,
     &      SLRX,BTGLB(M,NB),INDX )
!
!---      Aqueous and gas tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( SLB(M,NB),SGB(M,NB),
     &      PORDB(M,NB),TORLB(M,NB),TORGB(M,NB),N )
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
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     STOMPX-CO2
!
!     Compute the gas/aqueous capillary pressure from the aqueous
!     saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 Oct 2021.
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
      SUB_LOG(ISUB_LOG) = '/CAP_CO2'
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
            GOTO 222
          ENDIF
          CL = MAX( SCHR(3,N),SMALL )
          IF( SLX.GT.SCHR(4,N) ) THEN
            HDGL = SCHR(1,N)*(1.D+0/ESLX)**(1.D+0/CL)
          ELSE
            HDGL = SCHR(12,N)
            SLX = SLRX + 1.D-9
          ENDIF
  222     CONTINUE
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
  300   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,N)),ZERO )*REALX
!
!---    Matrix saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        DSLRX = -SCHR(4,N)/(HDGL*LOG(SCHR(12,N)))*REALX
        ASLX = 1.D+0/((1.D+0 + (SCHR(1,N)*HDGL)**CNM)**CMM)
        DASLX = -CMM*SCHR(1,N)*CNM*((SCHR(1,N)*HDGL)**(CNM-1.D+0))
     &  /((1.D+0 + (SCHR(1,N)*HDGL)**CNM)**(CMM+1.D+0))
        SLMZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLMZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
!
!---    Fracture saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(7,N),ZERO )
        DSLRX = -SCHR(7,N)/(HDGL*LOG(SCHR(12,N)))*REALX
        ASLX = 1.D+0/((1.D+0 + (SCHR(5,N)*HDGL)**CNF)**CMF)
        DASLX = -CMF*SCHR(5,N)*CNF*((SCHR(5,N)*HDGL)**(CNF-1.D+0))
     &  /((1.D+0 + (SCHR(5,N)*HDGL)**CNF)**(CMF+1.D+0))
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
          M_ERR(1) = 'Dual Porosity van Genuchten: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &      'Aqueous Saturation = '
          M_ERR(2) = ' at Node: '
          CALL PATH
          R_ERR = SLX
          I_ERR(1) = ND(N)
          I_ERR(2) = 1
          I_ERR(3) = 2
          I_ERR(4) = ID
!
!---      Use matrix properties to generate a guess for
!         capillary head  ---
!
          IF( SLX.GT.SCHR(4,N) ) THEN
            ASLX = (SLX-SCHR(4,N))/(1.D+0-SCHR(4,N))
            HDGL = (((1.D+0/ASLX)**(1.D+0/CMM)-1.D+0)**(1.D+0/CNM))/
     &        SCHR(1,N)
          ELSE
            HDGL = SCHR(12,N)
            SLX = SCHR(4,N) + 1.D-6
          ENDIF
          DH = 0.D+0
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 300
        CPGL = HDGL*RHORL*GRAV/BTGLX
!
!---  Dual porosity Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(N).EQ.4 ) THEN
        IF( (1.D+0-SLX)/EPSL.LT.EPSL ) THEN
          HDGL = SCHR(1,N)
          GOTO 410
        ENDIF
        CLM = MAX( SCHR(3,N),SMALL )
        CLF = MAX( SCHR(6,N),SMALL )
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
          HDGL = SCHR(1,N)*(1.D+0/ASLX)**(1.D+0/CLM)
        ELSE
          HDGL = SCHR(12,N)
          SLX = SCHR(4,N) + 1.D-9
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  400   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,N)),ZERO )*REALX
!
!---    Matrix saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        DSLRX = -SCHR(4,N)*REALX/(HDGL*LOG(SCHR(12,N)))
        HDGLX = MAX( SCHR(1,N),HDGL )
        ASLX = (SCHR(1,N)/HDGLX)**CLM
        DASLX = -CLM*(SCHR(1,N)/(HDGLX**2))
     &    *(SCHR(1,N)/HDGLX)**(CLM-1.D+0)
        SLMZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLMZ = DASLX*(1.D+0-SLRX) + DSLRX*(1.D+0-ASLX)
!
!---    Fracture saturation and partial derivative  ---
!
        SLRX = MAX( (1.D+0-HSCL)*SCHR(7,N),ZERO )
        DSLRX = -SCHR(7,N)*REALX/(HDGL*LOG(SCHR(12,N)))
        HDGLX = MAX( SCHR(5,N),HDGL )
        ASLX = (SCHR(5,N)/HDGLX)**CLF
        DASLX = -CLF*(SCHR(5,N)/(HDGLX**2))
     &    *(SCHR(5,N)/HDGLX)**(CLF-1.D+0)
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
          M_ERR(1) = 'Dual Porosity Brooks and Corey: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &      'Aqueous Saturation = '
          M_ERR(2) = ' at Node: '
          CALL PATH
          R_ERR = SLX
          I_ERR(1) = ND(N)
          I_ERR(2) = 1
          I_ERR(3) = 2
          I_ERR(4) = ID
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
          DH = 0.D+0
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 400
  410   CONTINUE
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
  500   CONTINUE
        NC = NC + 1
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,N)),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        DSLRX = -SCHR(4,N)*REALX/(HDGL*LOG(SCHR(12,N)))
        ASLX = SCHR(2,N)/(SCHR(2,N)+((HDGL-SCHR(1,N))/
     &    SCHR(5,N))**SCHR(3,N))
        DASLX = -(SCHR(2,N)*SCHR(3,N)*
     &    (((HDGL-SCHR(1,N))/SCHR(5,N))**(SCHR(3,N)-1.D+0))
     &    /SCHR(5,N))/((SCHR(2,N)+((HDGL-SCHR(1,N))/SCHR(5,N))
     &    **SCHR(3,N))**2)
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
          M_ERR(1) = 'Haverkamp: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &      'Aqueous Saturation = '
          M_ERR(2) = ' at Node: '
          CALL PATH
          R_ERR = SLX
          I_ERR(1) = ND(N)
          I_ERR(2) = 1
          I_ERR(3) = 2
          I_ERR(4) = ID
!
!---      Use matrix properties to generate a guess for
!         capillary head  ---
!
          IF( SLX.GT.SCHR(4,N) ) THEN
            ASLX = (SLX-SCHR(4,N))/(1.D+0-SCHR(4,N))
            HDGL = SCHR(1,N) + SCHR(5,N)*
     &        ((SCHR(2,N)/ASLX)-SCHR(2,N))**(1.D+0/SCHR(3,N))
          ELSE
            HDGL = SCHR(12,N)
            SLX = SCHR(4,N) + 1.D-9
          ENDIF
          DH = 0.D+0
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 500
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
 1100   CONTINUE
        NC = NC + 1
!
!---    No convergence on van Genuchten drainage-imbibition
!       saturation function  ---
!
        IF( NC.GT.32 ) THEN
          M_ERR(1) = 'van Genuchten Drainage-Imbibition: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &      'Aqueous Saturation = '
          M_ERR(2) = ' at Node: '
          CALL PATH
          R_ERR = SLX
          I_ERR(1) = ND(N)
          I_ERR(2) = 1
          I_ERR(3) = 2
          I_ERR(4) = ID
          HDGL = (ASLMINX/ASLX)*HDGLD + (1.D+0-(ASLMINX/ASLX))*HDGLI
          GOTO 1120
        ENDIF
        DO M = 1,2
          HDGLZ = HDGL
          DHDGLZ = SIGN( MAX( 1.D-6*HDGL,1.D-6 ),
     &      (5.D-1*SCHR(12,N) - HDGL) )
          IF( M.EQ.2 ) HDGLZ = HDGL + DHDGLZ
          HMPDZ = SCHR(9,N)
          HMPIZ = SCHR(11,N)
          SLRZ = SCHR(4,N)
!
!---      van Genuchten drainage-imbibition saturation function,
!         w/ Webb extension ASL = SL + SGT  ---
!
          IF( ISM(N).EQ.2 ) THEN
!
!---        Capillary head above the drainage matching point head,
!           use Webb extension  ---
!
            IF( HDGLZ.GT.HMPDZ ) THEN
              SMPDZ = SCHR(8,N)
              HDGLZ = MIN( HDGLZ,SCHR(12,N) )
              DMPDZ = SMPDZ/(LOG10(SCHR(12,N))-LOG10(HMPDZ))
              SLDZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPDZ
              ASLDZ = SLDZ
              ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
!
!---          Capillary head above the imbibition matching point head,
!             use Webb extension  ---
!
              IF( HDGLZ.GT.HMPIZ ) THEN
                SMPIZ = SCHR(10,N)
                DMPIZ = SMPIZ/(LOG10(SCHR(12,N))-LOG10(HMPIZ))
                SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPIZ
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
!---          Capillary head above the imbibition matching point head,
!             use Webb extension  ---
!
              IF( HDGLZ.GT.HMPIZ ) THEN
                SMPIZ = SCHR(10,N)
                DMPIZ = SMPIZ/(LOG10(SCHR(12,N))-LOG10(HMPIZ))
                SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPIZ
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
        ENDDO
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
 1200   CONTINUE
        NC = NC + 1
!
!---    No convergence on Brooks and Corey drainage-imbibition
!       saturation function  ---
!
        IF( NC.GT.32 ) THEN
          M_ERR(1) = 'Brooks and Corey Drainage-Imbibition: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &      'Aqueous Saturation = '
          M_ERR(2) = ' at Node: '
          CALL PATH
          R_ERR = SLX
          I_ERR(1) = ND(N)
          I_ERR(2) = 1
          I_ERR(3) = 2
          I_ERR(4) = ID
          HDGL = (ASLMINX/ASLX)*HDGLD + (1.D+0-(ASLMINX/ASLX))*HDGLI
          GOTO 1220
        ENDIF
        DO M = 1,2
          HDGLZ = HDGL
          DHDGLZ = SIGN( MAX( 1.D-6*HDGL,1.D-6 ),
     &      (5.D-1*SCHR(12,N) - HDGL) )
          IF( M.EQ.2 ) HDGLZ = HDGL + DHDGLZ
          HMPDZ = SCHR(9,N)
          HMPIZ = SCHR(11,N)
          SLRZ = SCHR(4,N)
!
!---      Brooks and Corey drainage-imbibition saturation function,
!         w/ Webb extension ASL = SL + SGT  ---
!
          IF( ISM(N).EQ.2 ) THEN
!
!---        Capillary head above the drainage matching point head,
!           use Webb extension  ---
!
            IF( HDGLZ.GT.HMPDZ ) THEN
              SMPDZ = SCHR(8,N)
              HDGLZ = MIN( HDGLZ,SCHR(12,N) )
              DMPDZ = SMPDZ/(LOG10(SCHR(12,N))-LOG10(HMPDZ))
              SLDZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPDZ
              ASLDZ = SLDZ
              ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
!
!---          Capillary head above the imbibition matching point head,
!             use Webb extension  ---
!
              IF( HDGLZ.GT.HMPIZ ) THEN
                SMPIZ = SCHR(10,N)
                DMPIZ = SMPIZ/(LOG10(SCHR(12,N))-LOG10(HMPIZ))
                SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPIZ
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
                CLI = MAX( SCHR(6,N),SMALL )
                IF( HDGLZ.LE.SCHR(5,N) ) THEN
                  ASLIZ = 1.D+0
                ELSE
                  ASLIZ = (SCHR(5,N)/HDGLZ)**CLI
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
              CLD = MAX( SCHR(3,N),SMALL )
              IF( HDGLZ.LE.SCHR(1,N) ) THEN
                ASLDZ = 1.D+0
              ELSE
                ASLDZ = (SCHR(1,N)/HDGLZ)**CLD
              ENDIF
              ASLDZ = ASLDZ*(1.D+0-SLRZ) + SLRZ
              ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
!
!---          Capillary head above the imbibition matching point head,
!             use Webb extension  ---
!
              IF( HDGLZ.GT.HMPIZ ) THEN
                SMPIZ = SCHR(10,N)
                DMPIZ = SMPIZ/(LOG10(SCHR(12,N))-LOG10(HMPIZ))
                SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,N)))*DMPIZ
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
                CLI = MAX( SCHR(6,N),SMALL )
                IF( HDGLZ.LE.SCHR(5,N) ) THEN
                  ASLIZ = 1.D+0
                ELSE
                  ASLIZ = (SCHR(5,N)/HDGLZ)**CLI
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
        ENDDO
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
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLUX_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Compute fluxes on internal and boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021
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
      SUB_LOG(ISUB_LOG) = '/FLUX_CO2'
!
!---  Gas volumetric flux (non-boundary surfaces)  ---
!
      CALL DRCVG
!
!---  Aqueous volumetric flux (non-boundary surfaces)  ---
!
      CALL DRCVL
!
!---  Salt diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-brine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS
!
!---  CO2 and water diffusive flux through gas
!     (non-boundary surfaces)  ---
!
      CALL DFFGAW
!
!---  CO2 and water diffusive flux through aqueous
!     (non-boundary surfaces)  ---
!
      CALL DFFLAW
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLUX_CO2 group
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
!     STOMPX-CO2
!
!     Establish hydrostatic boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 1 November 2021.
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
      DO K = 1,KC
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
          M_ERR(1) = 'Hydrostatic Boundary Condition: Salt ' //
     &      'Mass Fraction < 0.0 : XLS = '
          M_ERR(2) = ' at Node: '
          CALL PATH
          R_ERR = YLSX
          I_ERR(1) = ABS(N_DB)
          I_ERR(2) = 1
          I_ERR(3) = 2
          I_ERR(4) = ID
          YLSX = 0.D+0
        ELSEIF( YLSX.GT.XLSMX ) THEN
          M_ERR(1) = 'Hydrostatic Boundary Condition: Salt ' //
     &      'Mass Fraction > Solubility Limit : XLS = '
          M_ERR(2) = ' at Node: '
          CALL PATH
          R_ERR = YLSX
          I_ERR(1) = ABS(N_DB)
          I_ERR(2) = 1
          I_ERR(3) = 2
          I_ERR(4) = ID
          YLSX = XLSMX
        ENDIF
        PX = MAX( PGX,PLX )
        CHMSGX(1) = 'Unconverged Hydrostatic Boundary Conditions: ' //
     &    'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
        CHMSGX(2) = 'Hydrostatic Boundary Condition Transition: ' //
     &    'Vapor Pressure > Gas Pressure: PVB + PVA = '
        CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
        CALL DENS_B( TX,PLX,XLSX,RHOBX )
        CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
        PLX = PLX + RHOLX*GRAV*DISTZ
        PGX = PLX
      ENDDO
      PLX = PLX - PATM
      PGX = PLX
!
!---  Reset subroutine string sequence  ---
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
!     STOMPX-CO2
!
!     Compute primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 October 2021
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
      SUB_LOG(ISUB_LOG) = '/INCRM_CO2'
      PETA = 1.D-1
      EPSLX = 1.D-4
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = ND(N)
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
     &    BTGLX = SCHR(16,N)/SFTLX
!
!---    Assign gas entry pressure and minimum gas saturation
!       for transition to unsaturated conditions---
!
        SGMNX = 1.D-3
        ENPR = 0.D+0
        IF( ISCHR(N).EQ.1 .OR. ISCHR(N).EQ.101 .OR.
     &    ISCHR(N).EQ.201 ) THEN
          SGMNX = 1.D+1**(-3.D+0+LOG10(1.D+0/SCHR(1,N)))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ELSEIF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 .OR.
     &    ISCHR(N).EQ.202 ) THEN
          ENPR = SCHR(1,N)*RHORL*GRAV
          SGMNX = 1.D+1**(-3.D+0+LOG10(SCHR(1,N)))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ELSEIF( ISCHR(N).EQ.3 ) THEN
          SCHRX = MAX( SCHR(1,N),SCHR(5,N) )
          SGMNX = 1.D+1**(-3.D+0+LOG10(1.D+0/SCHRX))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ELSEIF( ISCHR(N).EQ.4 ) THEN
          ENPR = MIN( SCHR(1,N),SCHR(5,N) )*RHORL*GRAV
          SCHRX = MIN( SCHR(1,N),SCHR(5,N) )
          SGMNX = 1.D+1**(-3.D+0+LOG10(SCHRX))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ENDIF
!
!---    Maximum effective trapped gas saturation  ---
!
        IF( ISCHR(N).EQ.101 .OR. ISCHR(N).EQ.102 .OR.
     &    ISCHR(N).EQ.201 .OR. ISCHR(N).EQ.202 ) THEN
!
!---      w/ Webb extension  ---
!
          IF( ISM(N).EQ.2 ) THEN
            ESGTMX = SCHR(15,N)
!
!---      w/o Webb extension  ---
!
          ELSE
            ESGTMX = SCHR(15,N)/(1.D+0-SCHR(4,N))
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
     &          BTGL(2,N) = SCHR(16,N)/SFTLX
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
          CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVW(2,N),XLSX,N )
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
                  M_ERR(1) = 'Transition from Subcritical Liquid to '
     &              // 'Gas: Gas Pressure = '
                  M_ERR(2) = ' at Node: '
                  CALL PATH
                  R_ERR = PG(2,N)+PATM
                  I_ERR(1) = ND(N)
                  I_ERR(2) = 1
                  I_ERR(3) = 2
                  I_ERR(4) = ID
                ENDIF
                NPHAZ(2,N) = 3
              ELSEIF( I_VX.EQ.1 ) THEN
                IF( NPHAZ(2,N).EQ.2 ) THEN
                  M_ERR(1) = 'Transition from Subcritical Gas to '
     &              // 'Liquid: Gas Pressure = '
                  M_ERR(2) = ' at Node: '
                  CALL PATH
                  R_ERR = PG(2,N)+PATM
                  I_ERR(1) = ND(N)
                  I_ERR(2) = 1
                  I_ERR(3) = 2
                  I_ERR(4) = ID
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
     &        BTGL(2,N) = SCHR(16,N)/SFTLX
            SL(2,N) = 0.D+0
            SGT(2,N) = 0.D+0
            CALL CAP_CO2( ASLMIN(2,N),BTGL(2,N),PCX,SL(2,N),SGT(2,N),N )
            PL(2,N) = PG(2,N) - PCX
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL DENS_B( T(2,N),PX,XLSX,RHOBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVW(2,N),XLSX,N )
            PVA(2,N) = PGX - PVW(2,N)
!
!---        Check for transitions across the saturation line  ---
!
            CALL DENS_A( T(2,N),PVA(2,N),RHOAX,I_VX )
            IF( I_VX.EQ.0 ) THEN
              IF( NPHAZ(2,N).EQ.4 ) THEN
                M_ERR(1) = 'Transition from Subcritical Liquid to '
     &            // 'Gas: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
              ENDIF
              NPHAZ(2,N) = 8
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.2 ) THEN
                M_ERR(1) = 'Transition from Subcritical Gas to '
     &            // 'Liquid: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
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
                M_ERR(1) = 'Transition from Subcritical Liquid to '
     &            // 'Gas: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
              ENDIF
              NPHAZ(2,N) = 2
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.2 ) THEN
                M_ERR(1) = 'Transition from Subcritical Gas to '
     &            // 'Liquid: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
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
     &      BTGL(2,N) = SCHR(16,N)/SFTLX
          SLRX = SCHR(4,N)
!
!---      w/ Webb extension  ---
!
          IF( ISM(N).EQ.2 ) THEN
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
!---      Aqueous-gas 1  ---
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
            IF( ISM(N).EQ.2 ) THEN
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
                M_ERR(1) = 'Transition from Subcritical Liquid to '
     &            // 'Gas: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
              ENDIF
              NPHAZ(2,N) = 2
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.3 ) THEN
                M_ERR(1) = 'Transition from Subcritical Gas to '
     &            // 'Liquid: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
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
            IF( ISCHR(N).EQ.202 ) THEN
              PG(2,N) = PL(2,N) + (SCHR(5,N)/BTGLX)-EPSLX
            ELSE
              PG(2,N) = PL(2,N) + (ENPR/BTGLX)-EPSLX
            ENDIF
            CALL DENS_A( T(2,N),PGAX,RHOAX,I_VX )
            IF( I_VX.EQ.0 ) THEN
              IF( NPHAZ(2,N).EQ.5 ) THEN
                M_ERR(1) = 'Transition from Subcritical Liquid to '
     &            // 'Gas: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
              ENDIF
              NPHAZ(2,N) = 3
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.3 ) THEN
                M_ERR(1) = 'Transition from Subcritical Gas to '
     &            // 'Liquid: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
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
     &      BTGL(2,N) = SCHR(16,N)/SFTLX
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
          CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,N )
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
     &        BTGL(2,N) = SCHR(16,N)/SFTLX
            SGTX = 0.D+0
            CALL CAP_CO2( ASLMIN(2,N),BTGL(2,N),PCX,SL(2,N),SGTX,N )
            PL(2,N) = PG(2,N) - PCX
!
!---        Check for transitions across the saturation line  ---
!
            CALL DENS_A( T(2,N),PGAX,RHOAX,I_VX )
            IF( I_VX.EQ.0 ) THEN
              IF( NPHAZ(2,N).EQ.9 ) THEN
                M_ERR(1) = 'Transition from Subcritical Liquid to '
     &            // 'Gas: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
              ENDIF
              NPHAZ(2,N) = 2
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.8 ) THEN
                M_ERR(1) = 'Transition from Subcritical Gas to '
     &            // 'Liquid: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
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
                M_ERR(1) = 'Transition from Subcritical Liquid to '
     &            // 'Gas: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
              ENDIF
              NPHAZ(2,N) = 8
            ELSEIF( I_VX.EQ.1 ) THEN
              IF( NPHAZ(2,N).EQ.8 ) THEN
                M_ERR(1) = 'Transition from Subcritical Gas to '
     &            // 'Liquid: Gas Pressure = '
                M_ERR(2) = ' at Node: '
                CALL PATH
                R_ERR = PG(2,N)+PATM
                I_ERR(1) = ND(N)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
              ENDIF
              NPHAZ(2,N) = 9
            ELSE
              NPHAZ(2,N) = 10
            ENDIF
          ENDIF
        ENDIF
!
!---    Compute increments  ---
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
     &    BTGLX = SCHR(16,N)/SFTLX
!
!---    Assign gas-entry pressure for non Brooks-Corey;
!       Brooks-Corey; Brooks-Corey, Dual Porosity; and
!       Brooks-Corey, Entrapment  ---
!
        IF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 .OR.
     &      ISCHR(N).EQ.202 ) THEN
          ENPR = SCHR(1,N)*RHORL*GRAV
        ELSEIF( ISCHR(N).EQ.4 ) THEN
          ENPR = MIN( SCHR(1,N),SCHR(5,N) )*RHORL*GRAV
        ELSE
          ENPR = 0.D+0
        ENDIF
!
!---    Initial trapped gas saturation for the van Genuchten or
!       Brooks/Corey entrapment model  ---
!
        IF( ISCHR(N).EQ.101 .OR. ISCHR(N).EQ.102 .OR.
     &      ISCHR(N).EQ.201 .OR. ISCHR(N).EQ.202 ) THEN
!
!---      w/ Webb extension  ---
!
          IF( ISM(N).EQ.2 ) THEN
            ESGTMX = SCHR(15,N)
!
!---      w/o Webb extension  ---
!
          ELSE
            ESGTMX = SCHR(15,N)/(1.D+0-SCHR(4,N))
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
          IF( N.EQ.2 .AND. ID.EQ.0 ) THEN
            XLAMX = -1.D+0
          ENDIF
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
          IF( ISM(N).EQ.2 ) THEN
            SGTMX = ESGTMX
!
!---      w/o Webb extension  ---
!
          ELSE
            SGTMX = ESGTMX*(1.D+0-SCHR(4,N))
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
        DO M = 3,ISVC+2
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
              IF( ISCHR(N).EQ.202 ) THEN
                PG(M,N) = PL(M,N) + (SCHR(5,N)/BTGLX)-EPSLX
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
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     STOMPX-CO2
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
!
!---  Reset subroutine string sequence  ---
!
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
!     STOMPX-CO2
!
!     Load the current time step values into the old time step
!     variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 November 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
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
      SUB_LOG(ISUB_LOG) = '/LDO_CO2'
!
!---  Assign old time step values  ---
!
      DO N = 1,NFCGC(ID+1)
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
        RKL(1,N) = RKL(2,N)
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
        DO NSL = 1,NSOLU+NEQC+NEQK
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     STOMPX-CO2
!
!     Compute hydrologic, thermodynamic and physical properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 May 2002
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE LEAK_WELL
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
      SUB_LOG(ISUB_LOG) = '/PROP_CO2'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = ND(N)
!
!---    Initial trapped gas saturation for the van Genuchten or
!       Brooks/Corey entrapment model  ---
!
        IF( ISCHR(N).EQ.101 .OR. ISCHR(N).EQ.102 .OR.
     &      ISCHR(N).EQ.201 .OR. ISCHR(N).EQ.202 ) THEN
!
!---      w/ Webb extension  ---
!
          IF( ISM(N).EQ.2 ) THEN
            ESGTMX = SCHR(15,N)
!
!---      w/o Webb extension  ---
!
          ELSE
            ESGTMX = SCHR(15,N)/(1.D+0-SCHR(4,N))
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
            CALL VPL_B( T(M,N),PSBX,PCX,RHOBX,PVBX,XLSX,N )
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
            CALL VPL_B( T(M,N),PSBX,PCX,RHOBX,PVBX,XLSX,N )
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
            CALL VPL_B( T(M,N),PSBX,PCX,RHOBX,PVBX,XLSX,N )
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
     &        BTGL(M,N) = SCHR(16,N)/SFTLX
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
          CALL PORSTY_CO2E( PX,PCMP(N),PORD(M,N),PORT(M,N),N )
          PORD(M,N) = MAX( PORD(M,N),EPSL )
          PORT(M,N) = MAX( PORT(M,N),PORD(M,N) )
!
!---      Surface tension, saturation, relative permeability  ---
!
          CALL SFT_L( T(M,N),XLSX,SFTLX )
          BTGL(M,N) = 1.D+0
          IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &      BTGL(M,N) = SCHR(16,N)/SFTLX
          INDX = 0
          IF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.5
     &      .OR. NPHAZ(2,N).EQ.7 ) THEN
            SGT(M,N) = SG(M,N)
            INDX = 1
          ENDIF
          CALL KSP_CO2( N,PG(M,N),PL(M,N),SG(M,N),SGT(M,N),SL(M,N),
     &      SDPF(N),SDPM(N),RKL(M,N),RKG(M,N),ASLX,ASLMINX,
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
!
!---      Permeability reduction factor with simplified Verma and
!         Preuss model  ---
!
          IF( IPRF(N)/10.EQ.4 ) THEN
!
!---        Porosity with salt precipitation  ---
!
            PORD(M,N) = MAX( PORD(M,N)*(1.D+0-SS(M,N)),
     &        PORD(M,N)*PERM(5,N),1.D-12 )
!
!---        Permeability reduction factor with salt precipitation
!           and mineral precipitation or dissolution  ---
!
            PERMRF(M,N) = ((PORD(M,N)-PERM(5,N))/
     &        (PERM(4,N)-PERM(5,N)))**PERM(6,N)
!
!---      Permeability reduction factor with Verma and Preuss model  ---
!
          ELSEIF( IPRF(N).EQ.1 ) THEN
            CALL PERM_R( SS(M,N),PERMRF(M,N),PORD(M,N),N )
          ENDIF
!
!---      Aqueous and gas tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( SL(M,N),SG(M,N),
     &      PORD(M,N),TORL(M,N),TORG(M,N),N )
!
!---      Gas pressure for Darcy velocity calculation  ---
!
          PSO(M,N) = PG(M,N)
        ENDDO
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_BIN_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Read binary files from processor.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 September 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE REACT
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
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_BIN_CO2'
!
!---  Read solu.bin for solution control and output data  ---
!
      CALL READ_SOLU_CO2
!      PRINT *,'Post READ_SOLU_CO2: ID = ',ID
!
!---  Read a series of binary grid files for grid data  ---
!
      CALL READ_GRID_CO2
!      PRINT *,'Post READ_GRID_CO2: ID = ',ID
!
!---  Read mech.bin, hydr.bin satu.bin, perm.bin,
!     comp.bin, and tabl.bin for property data  ---
!
      CALL READ_PROP_CO2
!      PRINT *,'Post READ_PROP_CO2: ID = ',ID
!
!---  Read state.bin for state condition data  ---
!
      CALL READ_STATE_CO2
!      PRINT *,'Post READ_STATE_CO2: ID = ',ID
!
!---  Read boco.bin for boundary condition data  ---
!
      CALL READ_BOCO_CO2
!      PRINT *,'Post READ_BOCO_CO2: ID = ',ID
!
!---  Read sorc.bin for source data  ---
!
      CALL READ_SORC_CO2
!      PRINT *,'Post READ_SORC_CO2: ID = ',ID
!
!---  Read well.bin for coupled-well data ---
!
      CALL READ_COUP_WELL_CO2
!      PRINT *,'Post READ_COUP_WELL_CO2: ID = ',ID
!
!---  Read tpor.bin for transport data  ---
!
      IF( IEQC.NE.0 ) CALL READ_TPOR_CO2
!      PRINT *,'Post READ_TPOR_CO2: ID = ',ID
!
!---  Read geomechanial binary files  ---
!
      IF( ISLC(50).NE.0 ) CALL READ_BIN_GM
!      PRINT *,'Post READ_BIN_GM: ID = ',ID
!
!---  Read ecke.bin for ECKEChem (i.e., reactive transport data)  ---
!
      IF( ISLC(40).EQ.1 ) THEN
        CALL READ_REACT_CO2
      ELSE
        ALLOCATE( ISPLK(1:2*LNGC+LSPR+14),STAT=ISTAT )
        CHMSG = 'ISPLK'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        DO L = 1,2*LNGC+LSPR+14
          ISPLK(L) = 0
        ENDDO
      ENDIF
!      PRINT *,'Post READ_REACT_CO2: ID = ',ID
!
!---  Read ptcf.bin  ---
!
      IF( IACTV.EQ.2 ) CALL READ_PTZRCOEF
!      PRINT *,'Post READ_PTZR_CO2: IACTV = ',IACTV,'ID = ',ID
!
!---  Read co2_prop.bin  ---
!
      CALL READ_PE_CO2
!      PRINT *,'Post READ_PE_CO2: ID = ',ID
!
!---  Read restart.bin for restart simulations  ---
!
      IF( IEO.EQ.2 ) CALL RDRST_CO2
!      PRINT *,'Post RDRST_CO2: IEO = ',IEO,'ID = ',ID
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
!---  End of READ_BIN_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_BOCO_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Read binary boco_sorc.bin file for boundary condition data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 September 2021
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
      SUB_LOG(ISUB_LOG) = '/READ_BOCO_CO2'
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
!---  Array of initial condition boundary condition gas 
!     pressure  ---
!
      NVAR = NBC(ID+1)
      OFFSET = IOFFSET + NC*NBYTR + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NBC(ID+1)
        PGB(1,N) = VARX(N)
        VARX(N) = 0.D+0
      ENDDO
!
!---  Array of initial condition boundary condition aqueous salt 
!     mass fraction  ---
!
      NVAR = NBC(ID+1)
      OFFSET = IOFFSET + NC*NBYTR + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NBC(ID+1)
        YLSB(1,N) = VARX(N)
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
!---  End of READ_BOCO_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_COUP_WELL_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Read binary well.bin file for coupled-well data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 13 December 2021
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
      SUB_LOG(ISUB_LOG) = '/READ_COUP_WELL_CO2'
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
!     &    ' ID = ',ID
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
!        PRINT *,'LSOLU_CW = ',LSOLU_CW,'LN_CW = ',LN_CW,'ID = ',ID
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ISOLU_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!
!---    Read active species in coupled well
!       (duplicated across processors)  ---
!
        NVAR = LSPC_CW*LN_CW
!        PRINT *,'LSPC_CW = ',LSPC_CW,'LN_CW = ',LN_CW,'ID = ',ID
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ISPC_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!
!---    Read active component species in coupled well
!       (duplicated across processors)  ---
!
        NVAR = LEQC*LN_CW
!        PRINT *,'LEQC = ',LEQC,'LN_CW = ',LN_CW,'ID = ',ID
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,ISOLC_CW,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!
!---    Read active kinetic species in coupled well
!       (duplicated across processors)  ---
!
        NVAR = LEQK*LN_CW
!        PRINT *,'LEQK = ',LEQK,'LN_CW = ',LN_CW,'ID = ',ID
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
!---  End of READ_COUP_WELL_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_GRID_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Read a series of binary grid files for grid data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 September 2021
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
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER IVARX(10)
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_GRID_CO2'
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
!---  Allocate memory for NUKFP ---
!
      ALLOCATE( NUKFP(1:NP),STAT=ISTAT )
      CHMSG = 'NUKFP'
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
!---  Allocate memory for NUKTP ---
!
      ALLOCATE( NUKTP(1:NP),STAT=ISTAT )
      CHMSG = 'NUKTP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Open grid1.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'grid1.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
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
!---  End of READ_GRID_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_PE_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Read binary co2_prop.bin file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 September 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE NCG_PT
      USE GRID
      USE GLB_PAR
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(11) :: VARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_PE_CO2'
!
!---  Allocate memory for CO2 property variables  ---
!
      CALL ALLOC_NCG_PT
!
!---  Open co2_prop.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'co2_prop.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IFILE,IERR )
      INCG = 1
!
!---  Initialize cummulative offset  ---
!
      IOFFSET = 0
!
!---  Read number number of pressure points  ---
!
      NVAR = 1
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,IP_TA(INCG),NVAR,
     &  MPI_INTEGER,STATUS,IERR )
!
!---  Read pressure points  ---
!
      NVAR = IP_TA(INCG)
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,P_TA(1,INCG),NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read number of temperature points for each pressure point  ---
!
      NVAR = IP_TA(INCG)
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,IT_TA(1,INCG),NVAR,
     &  MPI_INTEGER,STATUS,IERR )
!
!---  Loop over the number of pressure points, reading temperatures  ---
!
      DO IPX = 1,IP_TA(INCG)
        NVAR = IT_TA(IPX,INCG)
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IFILE,OFFSET,T_TA(1,IPX,INCG),NVAR,
     &    MPI_REAL8,STATUS,IERR )
      ENDDO
!
!---  Loop over the number of pressure points, reading density  ---
!
      DO IPX = 1,IP_TA(INCG)
        NVAR = IT_TA(IPX,INCG)
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IFILE,OFFSET,RHO_TA(1,IPX,INCG),NVAR,
     &    MPI_REAL8,STATUS,IERR )
      ENDDO
!
!---  Loop over the number of pressure points, reading enthalpy  ---
!
      DO IPX = 1,IP_TA(INCG)
        NVAR = IT_TA(IPX,INCG)
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IFILE,OFFSET,H_TA(1,IPX,INCG),NVAR,
     &    MPI_REAL8,STATUS,IERR )
      ENDDO
!
!---  Loop over the number of pressure points, reading internal
!     energy  ---
!
      DO IPX = 1,IP_TA(INCG)
        NVAR = IT_TA(IPX,INCG)
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IFILE,OFFSET,U_TA(1,IPX,INCG),NVAR,
     &    MPI_REAL8,STATUS,IERR )
      ENDDO
!
!---  Loop over the number of pressure points, reading fugacity  ---
!
      DO IPX = 1,IP_TA(INCG)
        NVAR = IT_TA(IPX,INCG)
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IFILE,OFFSET,FUG_TA(1,IPX,INCG),NVAR,
     &    MPI_REAL8,STATUS,IERR )
      ENDDO
!
!---  Loop over the number of pressure points, reading entropy  ---
!
      DO IPX = 1,IP_TA(INCG)
        NVAR = IT_TA(IPX,INCG)
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IFILE,OFFSET,S_TA(1,IPX,INCG),NVAR,
     &    MPI_REAL8,STATUS,IERR )
      ENDDO
!
!---  Read index of vapor-liquid point  ---
!
      NVAR = IP_TA(INCG)
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,IV_TA(1,INCG),NVAR,
     &  MPI_INTEGER,STATUS,IERR )
!
!---  Read number of vapor-liquid data points  ---
!
      NVAR = 1
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,I_LV(INCG),NVAR,
     &  MPI_INTEGER,STATUS,IERR )
!
!---  Loop over the number of vapor-liquid data points, reading
!     temperature, pressure, liquid density, liquid enthalpy,
!     liquid internal energy, liquid entropy, vapor density,
!     vapor enthalpy, vapor internal energy, internal entropy,
!     vapor fugacity  ---
!
      DO IPX = 1,IP_TA(INCG)
        NVAR = 11
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IFILE,OFFSET,VARX,NVAR,
     &    MPI_REAL8,STATUS,IERR )
        T_LV(IPX,INCG) = VARX(1)
        P_LV(IPX,INCG) = VARX(2)
        RHOL_LV(IPX,INCG) = VARX(3)
        HL_LV(IPX,INCG) = VARX(4)
        UL_LV(IPX,INCG) = VARX(5)
        SL_LV(IPX,INCG) = VARX(6)
        RHOV_LV(IPX,INCG) = VARX(7)
        HV_LV(IPX,INCG) = VARX(8)
        UV_LV(IPX,INCG) = VARX(9)
        SV_LV(IPX,INCG) = VARX(10)
        FUG_LV(IPX,INCG) = VARX(11)
      ENDDO
!
!---  Close the co2_prop.bin file  ---
!
      CALL MPI_FILE_CLOSE( IFILE,IERR )
!
!---  Load pressure triple-point and critical-point indices  ---
!
      DO I=1,IP_TA(INCG)
        IF( ABS(P_TA(I,INCG)-(1.D-6*PTPA)).LT.EPSL ) IPTP(INCG) = I
        IF( ABS(P_TA(I,INCG)-(1.D-6*PCRA)).LT.EPSL ) IPCR(INCG) = I
      ENDDO
!
!---  Vapor temperature splines  ---
!
      DO IPX = 1,IP_TA(INCG)
        N = IT_TA(IPX,INCG)-IV_TA(IPX,INCG)+1
        CALL SPLINE( T_TA(IV_TA(IPX,INCG),IPX,INCG),
     &    RHO_TA(IV_TA(IPX,INCG),IPX,INCG),N,
     &    RHO_ST(IV_TA(IPX,INCG),IPX,INCG) )
        CALL SPLINE( T_TA(IV_TA(IPX,INCG),IPX,INCG),
     &    H_TA(IV_TA(IPX,INCG),IPX,INCG),N,
     &    H_ST(IV_TA(IPX,INCG),IPX,INCG) )
        CALL SPLINE( T_TA(IV_TA(IPX,INCG),IPX,INCG),
     &    U_TA(IV_TA(IPX,INCG),IPX,INCG),N,
     &    U_ST(IV_TA(IPX,INCG),IPX,INCG) )
        CALL SPLINE( T_TA(IV_TA(IPX,INCG),IPX,INCG),
     &    FUG_TA(IV_TA(IPX,INCG),IPX,INCG),N,
     &    FUG_ST(IV_TA(IPX,INCG),IPX,INCG) )
        CALL SPLINE( T_TA(IV_TA(IPX,INCG),IPX,INCG),
     &    S_TA(IV_TA(IPX,INCG),IPX,INCG),N,
     &    S_ST(IV_TA(IPX,INCG),IPX,INCG) )
       ENDDO
!
!---  Liquid temperature splines  ---
!
      DO IPX = 1,IP_TA(INCG)
        N = IV_TA(IPX,INCG)-1
        IF( N.GT.0 ) THEN
          CALL SPLINE( T_TA(1,IPX,INCG),
     &      RHO_TA(1,IPX,INCG),N,
     &      RHO_ST(1,IPX,INCG) )
          CALL SPLINE( T_TA(1,IPX,INCG),
     &      H_TA(1,IPX,INCG),N,
     &      H_ST(1,IPX,INCG) )
          CALL SPLINE( T_TA(1,IPX,INCG),
     &      U_TA(1,IPX,INCG),N,
     &      U_ST(1,IPX,INCG) )
          CALL SPLINE( T_TA(1,IPX,INCG),
     &      FUG_TA(1,IPX,INCG),N,
     &      FUG_ST(1,IPX,INCG) )
          CALL SPLINE( T_TA(1,IPX,INCG),
     &      S_TA(1,IPX,INCG),N,
     &      S_ST(1,IPX,INCG) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_PE_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_PROP_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Read binary mech.bin, hydr.bin, satu1.bin,  satu2.bin, perm.bin,
!     comp.bin, and tabl.bin files for property data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 September 2021
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
      SUB_LOG(ISUB_LOG) = '/READ_PROP_CO2'
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
     &  STATUS,IERR )
!
!---  Read local copies of RHOS array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NC*NBYTR + NBYTB
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,RHOS,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of CPS array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NC*NBYTR + NBYTB
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,CPS,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of CMP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*4
      OFFSET = IOFFSET + NBYTB + NC*4*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*4*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,CMP,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of POR array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*6
      OFFSET = IOFFSET + NBYTB + NC*6*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*6*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,POR,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of TOR array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*6
      OFFSET = IOFFSET + NBYTB + NC*6*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*6*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TOR,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of ITOR array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ITOR,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
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
     &  STATUS,IERR )
!
!---  Read local copies of IPRF array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IPRF,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
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
     &  STATUS,IERR )
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
     &  STATUS,IERR )
!
!---  Read local copies of ISM array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read local copies of ISLTBL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*4
      OFFSET = IOFFSET + NBYTB + NC*4*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*4*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISLTBL,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read local copies of IRGTBL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*4
      OFFSET = IOFFSET + NBYTB + NC*4*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*4*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IRGTBL,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read local copies of IRLTBL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*4
      OFFSET = IOFFSET + NBYTB + NC*4*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*4*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IRLTBL,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
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
!---  Read local copies of RPGC array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*LRPGC
      OFFSET = IOFFSET + NBYTB + NC*LRPGC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*LRPGC*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,RPGC,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of IRPG array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IRPG,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read local copies of RPLC array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*LRPLC
      OFFSET = IOFFSET + NBYTB + NC*LRPLC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*LRPLC*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,RPLC,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of IRPL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IRPL,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
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
     &  STATUS,IERR )
!
!---  Read local copies of TCMP array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TCMP,NVAR,MPI_REAL8,
     &  STATUS,IERR )
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
     &  STATUS,IERR )
!
!---  Read local copies of TBLX array  ---
!
      NVAR = LTBL
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TBLX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of TBLY array  ---
!
      NVAR = LTBL
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TBLY,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of TBLDDX array  ---
!
      NVAR = LTBL
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TBLDDX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of TBLDDY array  ---
!
      NVAR = LTBL
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,TBLDDY,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Close the table.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_PROP_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_REACT_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Read binary ecke.bin file for ECKEChem (reactive transport) 
!     data
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 January 2022
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
      SUB_LOG(ISUB_LOG) = '/READ_REACT_CO2'
!
!---  Open ecke.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'ecke.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
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
     &  STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &      STATUS,IERR )
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
     &      STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &    STATUS,IERR )
!        PRINT *,'RC_K = ',(RC_K(K,1,1),K=1,LSPK+11),' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 3*LSPL
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SP_L,NVAR,MPI_REAL8,
     &    STATUS,IERR )
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
     &    STATUS,IERR )
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
     &        STATUS,IERR )
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
     &      STATUS,IERR )
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
     &      STATUS,IERR )
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
     &      STATUS,IERR )
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
     &    STATUS,IERR )
!        PRINT *,'RCNME = ',(RCNME(L),L=1,LRCE),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LRCK
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,RCNMK,NVAR,MPI_CHAR,
     &    STATUS,IERR )
!        PRINT *,'RCNMK = ',(RCNMK(L),L=1,LRCK),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LEQC
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNMC,NVAR,MPI_CHAR,
     &    STATUS,IERR )
!        PRINT *,'SPNMC = ',(SPNMC(L),L=1,LEQC),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LEQK
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNMK,NVAR,MPI_CHAR,
     &    STATUS,IERR )
!        PRINT *,'SPNMK = ',(SPNMK(L),L=1,LEQK),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LSPG
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNMG,NVAR,MPI_CHAR,
     &    STATUS,IERR )
!        PRINT *,'SPNMG = ',(SPNMG(L),L=1,LSPG),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LSPL
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNML,NVAR,MPI_CHAR,
     &    STATUS,IERR )
!        PRINT *,'SPNML = ',(SPNML(L),L=1,LSPL),' ID = ',ID
!
!---    Read character string for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LSPS
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNMS,NVAR,MPI_CHAR,
     &    STATUS,IERR )
!        PRINT *,'SPNMS = ',(SPNMS(L),L=1,LSPS),' ID = ',ID
!
!---    Read parameter for reactive species
!       (duplicated across processors)  ---
!
        NVAR = 64*LSPE
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SPNME,NVAR,MPI_CHAR,
     &    STATUS,IERR )
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
!---  Close the ecke.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_REACT_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_SOLU_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Read binary grid.bin file for solution control and output 
!     data
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 September 2021
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
      SUB_LOG(ISUB_LOG) = '/READ_SOLU_CO2'
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
!---  End of READ_SOLU_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_SORC_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Read binary boco_sorc.bin file for source data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 September 2021
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
      SUB_LOG(ISUB_LOG) = '/READ_SORC_CO2'
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
     &  STATUS,IERR )
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
     &  STATUS,IERR )
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
     &  STATUS,IERR )
!
!---  Index array of source number of time points  ---
!
      NVAR = NSR(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISRM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Index array of source inputs  ---
!
      NVAR = NSR(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISRIN,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Index array of source types  ---
!
      NVAR = NSR(ID+1)
      OFFSET = IOFFSET + NC*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ISRT,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      DO N = 1,NSR(ID+1)
!        PRINT *,'ISRT(',N,') = ',ISRT(N),' ID = ',ID
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
!---  End of READ_SORC_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_STATE_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Read binary state.bin file for state condition data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 September 2021
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
      SUB_LOG(ISUB_LOG) = '/READ_STATE_CO2'
!
!---  Open state1.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'state1.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Allocate local memory state condition arrays
!     (including ghost cells)  ---
!
      CALL ALLOC_FDVP
!
!---  Initialize global array memory for general field variables  ---
!
      CALL INTLZ_FDVP
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
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        T(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of PL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        PL(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of PG array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        PG(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of SG array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        SG(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of SL array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        SL(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of SGT array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        SGT(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of ASLMIN array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        ASLMIN(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of TMS array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        TMS(2,N) = VARX(N)
      ENDDO
!
!---  Close the state1.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Open state2.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'state2.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read local copies of YLS array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        YLS(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of PVA array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        PVA(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of PVW array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        PVW(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of XLA array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        XLA(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of PORD array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        PORD(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of PORT array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        PORT(2,N) = VARX(N)
      ENDDO
!
!---  Read local copies of NPHAZ array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IVARX,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        NPHAZ(2,N) = IVARX(N)
      ENDDO
!
!---  Read local copies of SI (Eclipse gas saturation)
!     array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        SI(1,N) = VARX(N)
      ENDDO
!
!---  Read local copies of SI array (Eclipse pressure)
!     (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,VARX,NVAR,MPI_REAL8,
     &  STATUS,IERR )
      DO N = 1,NFCGC(ID+1)
        SI(2,N) = VARX(N)
      ENDDO
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
!---  Close the state2.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_STATE_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_TPOR_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Read binary well_tport.bin file for solute transport data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 12 January 2022
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
      SUB_LOG(ISUB_LOG) = '/READ_TPOR_CO2'
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
!---    Read solute gas diffusion coefficient
!       (duplicated across processors)  ---
!
        NVAR = NSOLU
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,SMDG,NVAR,MPI_REAL8,
     &    STATUS,IERR )
!        PRINT *,'SMDG = ',(SMDG(NSL),NSL=1,NSOLU),' ID = ',ID
!
!---    Read index for solute gas-aqueous partition coefficient
!       (duplicated across processors)  ---
!
        NVAR = NSOLU
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,IPCGL,NVAR,MPI_INTEGER,
     &    STATUS,IERR )
!        PRINT *,'IPCGL = ',(IPCGL(NSL),NSL=1,NSOLU),' ID = ',ID
!
!---    Read solute gas-aqueous partition coefficient parameters
!       (duplicated across processors)  ---
!
        NVAR = 5*NSOLU
        OFFSET = IOFFSET + NBYTB
        IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
        CALL MPI_FILE_READ_AT( IRD,OFFSET,PCGL,NVAR,MPI_REAL8,
     &    STATUS,IERR )
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
!---  End of READ_TPOR_CO2 group  ---
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
!     STOMPX-CO2
!
!     Gas relative permeability.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 October 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
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
      ASLZ = MAX( ASLX,0.D+0 )
!
!---  Gas relative permeability function  ---
!
      IF( IRPG(N).EQ.0 ) THEN
        RKGX = RPGC(3,N)
!
!---  Mualem porosity distribution function  ---
!
      ELSEIF( IRPG(N).EQ.1 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(N).EQ.1 .OR. ISCHR(N).EQ.101 ) THEN
          ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
          SGP = 1.D+0 - ASLZ
          RKGX = SQRT(SGP)*((1.D+0-ASLZ**(1.D+0/RPGC(3,N)))
     &      **RPGC(3,N))**2
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 .OR.
     &    ISCHR(N).EQ.202 ) THEN
          ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
          SGP = 1.D+0 - ASLZ
          RKGX = SQRT(SGP)*(1.D+0-ASLZ**(1.D+0+1.D+0/RPGC(3,N)))**2
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
        ELSEIF( ISCHR(N).EQ.3 ) THEN
          SGPM = 1.D+0-ASL_M
          RKGM = SQRT(SGPM)*((1.D+0-ASL_M**(1.D+0/RPGC(3,N)))
     &      **RPGC(3,N))**2
          SGPF = 1.D+0-ASL_F
          RKGF = SQRT(SGPF)*((1.D+0-ASL_F**(1.D+0/RPGC(1,N)))
     &      **RPGC(1,N))**2
          RKGX = ( PERM(4,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKGF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
          RKGY = ( PERM(5,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(8,N)*RKGF*POR(4,N) )/
     &      ( PERM(5,N)*(1.D+0-POR(4,N)) + PERM(8,N)*POR(4,N)
     &      + SMALL )
          RKGZ = ( PERM(6,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(9,N)*RKGF*POR(4,N) )/
     &      ( PERM(6,N)*(1.D+0-POR(4,N)) + PERM(9,N)*POR(4,N)
     &      + SMALL )
          RKGX = MAX( RKGX,RKGY,RKGZ )
!
!---    Single-pressure dual-porosity Brooks and Corey  ---
!
        ELSEIF( ISCHR(N).EQ.4 ) THEN
          SGPM = 1.D+0-ASL_M
          RKGM = SGPM**(2.5D+0 + 2.0D+0/RPGC(3,N))
          SGPF = 1.D+0-ASL_F
          RKGF = SGPF**(2.5D+0 + 2.0D+0/RPGC(1,N))
          RKGX = ( PERM(4,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKGF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
          RKGY = ( PERM(5,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(8,N)*RKGF*POR(4,N) )/
     &      ( PERM(5,N)*(1.D+0-POR(4,N)) + PERM(8,N)*POR(4,N)
     &      + SMALL )
          RKGZ = ( PERM(6,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(9,N)*RKGF*POR(4,N) )/
     &      ( PERM(6,N)*(1.D+0-POR(4,N)) + PERM(9,N)*POR(4,N)
     &      + SMALL )
          RKGX = MAX( RKGX,RKGY,RKGZ )
        ENDIF
!
!---  Burdine porosity distribution function  ---
!
      ELSEIF( IRPG(N).EQ.2 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(N).EQ.1 .OR. ISCHR(N).EQ.101 ) THEN
          ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
          SGP = 1.D+0 - ASLZ
          RKGX = (SGP**2)*((1.D+0-ASLZ**(1.D+0/RPGC(3,N)))
     &      **RPGC(3,N))
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 .OR.
     &    ISCHR(N).EQ.202 ) THEN
          ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
          SGP = 1.D+0 - ASLZ
          RKGX = (SGP**2)*(1.D+0-ASLZ**(1.D+0 + 2.D+0/RPGC(3,N)))
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
        ELSEIF( ISCHR(N).EQ.3 ) THEN
          SGPM = 1.D+0-ASL_M
          RKGM = (SGPM**2)*((1.D+0-ASL_M**(1.D+0/RPGC(3,N)))
     &      **RPGC(3,N))
          SGPF = 1.D+0-ASL_F
          RKGF = (SGPF**2)*((1.D+0-ASL_F**(1.D+0/RPGC(1,N)))
     &      **RPGC(1,N))
          RKGX = ( PERM(4,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKGF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
          RKGY = ( PERM(5,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(8,N)*RKGF*POR(4,N) )/
     &      ( PERM(5,N)*(1.D+0-POR(4,N)) + PERM(8,N)*POR(4,N)
     &      + SMALL )
          RKGZ = ( PERM(6,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(9,N)*RKGF*POR(4,N) )/
     &      ( PERM(6,N)*(1.D+0-POR(4,N)) + PERM(9,N)*POR(4,N)
     &      + SMALL )
          RKGX = MAX( RKGX,RKGY,RKGZ )
!
!---    Single-pressure dual-porosity Brooks and Corey  ---
!
        ELSEIF( ISCHR(N).EQ.4 ) THEN
          SGPM = 1.D+0-ASL_M
          RKGM = SGPM**(3.0D+0 + 2.0D+0/RPGC(3,N))
          SGPF = 1.D+0-ASL_F
          RKGF = SGPF**(3.0D+0 + 2.0D+0/RPGC(1,N))
          RKGX = ( PERM(4,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKGF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
          RKGY = ( PERM(5,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(8,N)*RKGF*POR(4,N) )/
     &      ( PERM(5,N)*(1.D+0-POR(4,N)) + PERM(8,N)*POR(4,N)
     &      + SMALL )
          RKGZ = ( PERM(6,N)*RKGM*(1.D+0-POR(4,N)) +
     &      PERM(9,N)*RKGF*POR(4,N) )/
     &      ( PERM(6,N)*(1.D+0-POR(4,N)) + PERM(9,N)*POR(4,N)
     &      + SMALL )
          RKGX = MAX( RKGX,RKGY,RKGZ )
        ENDIF
!
!---  Modified-Corey relative permeability function  ---
!
      ELSEIF( IRPG(N).EQ.3 ) THEN
        SLRX = RPGC(1,N)
        SGRX = RPGC(3,N)
        SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SGRX),0.D+0 ),1.D+0 )
        ASLPX = SLPX + ESGTX
        RKGX = RPGC(2,N)*((1.D+0-ASLPX)**2)*(1.D+0-ASLPX**2)
!
!---  Fatt and Klikoff relative permeability function  ---
!
      ELSEIF( IRPG(N).EQ.4 ) THEN
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        RKGX = (1.D+0-ASLZ)**3
!
!---  Free-Corey relative permeability function  ---
!
      ELSEIF( IRPG(N).EQ.7 ) THEN
        SLRX = RPGC(3,N)
        SGRX = RPGC(4,N)
        SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SGRX),0.D+0 ),1.D+0 )
        ASLPX = SLPX + ESGTX
        SGPX = 1.D+0 - ASLPX
        RKGX = RPGC(1,N)*(SGPX**(RPGC(2,N)))
!
!---  Classical-Corey relative permeability function  ---
!
      ELSEIF( IRPG(N).EQ.17 ) THEN
        SLRX = RPGC(3,N)
        SGRX = RPGC(4,N)
        SLPX = (SLX-SLRX)/(1.D+0-SLRX-SGRX)
        SLPX = MIN( MAX( SLPX,0.D+0 ),1.D+0 )
        ASLPX = SLPX + ESGTX
        RKGX = RPGC(1,N)*(1.D+0-(ASLPX**(RPGC(2,N)/2.D+0)))*
     &    ((1.D+0-ASLPX)**(RPGC(2,N)/2.D+0))
!
!---  van Genuchten gas relative permeability function  ---
!
      ELSEIF( IRPG(N).EQ.9 ) THEN
        SLRX = RPGC(1,N)
        SGRX = RPGC(3,N)
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        SGZ = (1.D+0-ASLZ)*(1.D+0-SLRX)
        SGPX = MIN( MAX( (SGZ-SGRX)/(1.D+0-SGRX),0.D+0 ),1.D+0 )
        RKGX = SQRT(SGPX)*((1.D+0-((1.D+0-SGPX)**(1.D+0/RPGC(2,N)))
     &      **RPGC(2,N))**2)
!
!---  Tabular gas relative permeability versus gas saturation
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( IRPG(N).EQ.10 ) THEN
        SGZ = MAX( SGX-SGTX,0.D+0 )
        ITBX = 0
        RKGX = FNTBLY( SGZ,IRGTBL(1,N),IRGTBL(2,N),ITBX )
!
!---  Tabular gas relative permeability versus gas saturation
!     w/ cubic spline interpolation  ---
!
      ELSEIF( IRPG(N).EQ.11 ) THEN
        SGZ = MAX( SGX-SGTX,0.D+0 )
        RKGX = FSPLNY( SGZ,IRGTBL(1,N),IRGTBL(2,N) )
!
!---  Tabular gas relative permeability versus ln(capillary head)
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( IRPG(N).EQ.12 ) THEN
        ITBX = 0
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKGX = FNTBLY( HDGLX,IRGTBL(1,N),IRGTBL(2,N),ITBX )
!
!---  Tabular gas relative permeability versus ln(capillary head)
!     w/ cubic spline interpolation  ---
!
      ELSEIF( IRPG(N).EQ.13 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKGX = FSPLNY( HDGLX,IRGTBL(1,N),IRGTBL(2,N) )
!
!---  Tabular gas relative permeability versus capillary head
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( IRPG(N).EQ.14 ) THEN
        ITBX = 0
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKGX = FNTBLY( HDGL,IRGTBL(1,N),IRGTBL(2,N),ITBX )
!
!---  Tabular gas relative permeability versus capillary head
!     w/ cubic spline interpolation  ---
!
      ELSEIF( IRPG(N).EQ.15 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKGX = FSPLNY( HDGL,IRGTBL(1,N),IRGTBL(2,N) )
!
!---  Doughty gas relative permeability function  ---
!
      ELSEIF( IRPG(N).EQ.18 ) THEN
        RKGMX = RPGC(1,N)
        GAMMAX = RPGC(2,N)
        SLRX = RPGC(3,N)
        CMX = RPGC(4,N)
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        ASLPX = MIN( MAX( (ESLX+ESGTX),0.D+0 ),1.D+0 )
        RKGX = RKGMX*((1.D+0-ASLPX)**GAMMAX)*
     &    (((1.D+0-ASLPX)**(1.D+0/CMX))**(2.D+0*CMX))
!
!---  Doughty drainage-imbibition gas relative permeability
!     function  ---
!
      ELSEIF( IRPG(N).EQ.19 ) THEN
        RKGMX = RPGC(1,N)
        GAMMAX = RPGC(2,N)
        SLRDX = RPGC(3,N)
        CMX = RPGC(4,N)
        SLRIX = RPGC(5,N)
        SGTMX = RPGC(6,N)
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
      ELSEIF( IRPG(N).EQ.20 ) THEN
        SLRX = RPGC(1,N)
        SGRX = RPGC(3,N)
        SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SGRX),0.D+0 ),1.D+0 )
        ASLPX = SLPX + ESGTX
        IF( SLRX.GT.EPSL ) THEN
          BX = (1.D+0-RPGC(1,N))/SLRX
          AX = 1.D+0 - BX
          DSLX = 5.D-2
          RKG1X = AX + BX*(1.D+0-SLX)
          RKG2X = RPGC(2,N)*((1.D+0-ASLPX)**2)*
     &      (1.D+0-ASLPX**RPGC(4,N))
          RKG1X = MIN( MAX( RKG1X,0.D+0 ),1.D+0 )
          RKG2X = MIN( MAX( RKG2X,0.D+0 ),1.D+0 )
          SIGMAX = 1.D+0/(1.D+0 + EXP((SLRX-SLX)/DSLX))
          RKGX = (1.D+0-SIGMAX)*RKG1X + SIGMAX*RKG2X
        ELSE
          RKGX = RPGC(2,N)*((1.D+0-ASLPX)**2)*
     &      (1.D+0-ASLPX**RPGC(4,N))
        ENDIF

!
!---  Extended power law relative permeability function  ---
!
      ELSEIF( IRPG(N).EQ.21 ) THEN
        SLRX = RPGC(3,N)
        SGRX = RPGC(4,N)
        SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SGRX),0.D+0 ),1.D+0 )
        ASLPX = SLPX + ESGTX
        SGPX = 1.D+0 - ASLPX
        IF( SLRX.GT.EPSL ) THEN
          BX = (1.D+0-RPGC(1,N))/SLRX
          AX = 1.D+0 - BX
          DSLX = 5.D-2
          RKG1X = AX + BX*(1.D+0-SLX)
          RKG2X = RPGC(1,N)*(SGPX**(RPGC(2,N)))
          RKG1X = MIN( MAX( RKG1X,0.D+0 ),1.D+0 )
          RKG2X = MIN( MAX( RKG2X,0.D+0 ),1.D+0 )
          SIGMAX = 1.D+0/(1.D+0 + EXP((SLRX-SLX)/DSLX))
          RKGX = (1.D+0-SIGMAX)*RKG1X + SIGMAX*RKG2X
        ELSE
          RKGX = RPGC(1,N)*(SGPX**(RPGC(2,N)))
        ENDIF
!
!---  Mualem porosity distribution function w/ van Genuchten 
!     saturation function  ---
!
      ELSEIF( IRPG(N).EQ.41 ) THEN
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        SGP = 1.D+0 - ASLZ
        RKGX = SQRT(SGP)*((1.D+0-ASLZ**(1.D+0/RPGC(3,N)))**RPGC(3,N))**2
!
!---  Mualem porosity distribution function w/ Brooks and Corey 
!     saturation function  ---
!
      ELSEIF( IRPG(N).EQ.42 ) THEN
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        SGP = 1.D+0 - ASLZ
        RKGX = SQRT(SGP)*(1.D+0-ASLZ**(1.D+0+1.D+0/RPGC(3,N)))**2
!
!---  Burdine porosity distribution function w/ van Genuchten 
!     saturation function  ---
!
      ELSEIF( IRPG(N).EQ.43 ) THEN
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        SGP = 1.D+0 - ASLZ
        RKGX = (SGP**2)*((1.D+0-ASLZ**(1.D+0/RPGC(3,N)))**RPGC(3,N))
!
!---  Burdine porosity distribution function w/ Brooks and Corey 
!     saturation function  ---
!
      ELSEIF( IRPG(N).EQ.44 ) THEN
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        SGP = 1.D+0 - ASLZ
        RKGX = (SGP**2)*(1.D+0-ASLZ**(1.D+0 + 2.D+0/RPGC(3,N)))
      ENDIF
      RKGX = MIN( MAX( RKGX,0.D+0 ),1.D+0 )
!
!---  Reset subroutine string sequence  ---
!
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
!     STOMPX-CO2
!
!     Aqueous relative permeability.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 June 2011.
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
      SUB_LOG(ISUB_LOG) = '/RKLS_CO2'
      ESLZ = MAX( ESLX,0.D+0 )
!
!---  Constant relative permeability function  ---
!
      IF( MOD( IRPL(N),100 ).EQ.0 ) THEN
        RKLX = RPLC(2,N)
!
!---    Single-pressure dual-porosity saturation functions  ---
!
        IF( ISCHR(N).EQ.3 .OR. ISCHR(N).EQ.4 ) THEN
          RKLM = RPLC(2,N)
          RKLF = RPLC(1,N)
          RKLX = ( PERM(4,N)*RKLM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKLF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
        ENDIF
!
!---  Mualem-irreducible porosity distribution function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.21 ) THEN
        SLRX = RPLC(1,N)
        SLPX = MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 )
        IF( ISCHR(N).EQ.1 .OR. ISCHR(N).EQ.101  ) THEN
          RKLX = SQRT(SLPX)*((1.D+0-(1.D+0-SLPX**(1.D+0/RPLC(2,N)))
     &      **RPLC(2,N))**2)
        ELSEIF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 .OR.
     &    ISCHR(N).EQ.202 ) THEN
          RKLX = SLPX**(2.5D+0 + 2.0D+0/RPLC(2,N))
        ENDIF
!
!---  Mualem porosity distribution function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.1 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(N).EQ.1 .OR. ISCHR(N).EQ.101 ) THEN
          RKLX = SQRT(ESLZ)*((1.D+0-(1.D+0-ESLZ**(1.D+0/RPLC(2,N)))
     &      **RPLC(2,N))**2)
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 .OR.
     &    ISCHR(N).EQ.202 ) THEN
          RKLX = ESLZ**(2.5D+0 + 2.0D+0/RPLC(2,N))
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
          ELSEIF( ISCHR(N).EQ.3  ) THEN
            RKLM = (1.D+0-(ASL_M**(1.D+0/RPLC(2,N))))
            RKLM = SQRT(ASL_M)*((1.D+0-RKLM**RPLC(2,N))**2)
            RKLF = (1.D+0-(ASL_F**(1.D+0/RPLC(1,N))))
            RKLF = SQRT(ASL_F)*((1.D+0-RKLF**RPLC(1,N))**2)
            RKLX = ( PERM(4,N)*RKLM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKLF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
!
!---    Single-pressure dual-porosity Brooks and Corey  ---
!
          ELSEIF( ISCHR(N).EQ.4  ) THEN
            RKLM = ASL_M**(2.5D+0 + 2.0D+0/RPLC(2,N))
            RKLF = ASL_F**(2.5D+0 + 2.0D+0/RPLC(1,N))
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
          RKLX = (ESLZ**2)*(1.D+0-(1.D+0-ESLZ**(1.D+0/RPLC(2,N)))
     &      **RPLC(2,N))
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 .OR.
     &    ISCHR(N).EQ.202 ) THEN
          RKLX = ESLZ**(3.0D+0 + 2.0D+0/RPLC(2,N))
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
        ELSEIF( ISCHR(N).EQ.3 ) THEN
          RKLM = (ASL_M**2)*(1.D+0-(1.D+0-ASL_M**(1.D+0/RPLC(2,N)))
     &      **RPLC(2,N))
          RKLF = (ASL_F**2)*(1.D+0-(1.D+0-ASL_F**(1.D+0/RPLC(1,N)))
     &      **RPLC(1,N))
          RKLX = ( PERM(4,N)*RKLM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKLF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
!
!---  Dual porosity Brooks and Corey  ---
!
        ELSEIF( ISCHR(N).EQ.4 ) THEN
          RKLM = ASL_M**(3.0D+0 + 2.0D+0/RPLC(2,N))
          RKLF = ASL_F**(3.0D+0 + 2.0D+0/RPLC(1,N))
          RKLX = ( PERM(4,N)*RKLM*(1.D+0-POR(4,N)) +
     &      PERM(7,N)*RKLF*POR(4,N) )/
     &      ( PERM(4,N)*(1.D+0-POR(4,N)) + PERM(7,N)*POR(4,N)
     &      + SMALL )
        ENDIF
!
!---  Corey relative permeability function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.3 ) THEN
        RKLX = ESLZ**4
!
!---  Fatt and Klikoff relative permeability function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.4 ) THEN
        RKLX = ESLZ**3
!
!---  Haverkamp relative permeability function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.5 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 .OR.
     &    ISCHR(N).EQ.202 .OR. ISCHR(N).EQ.5 ) THEN
          HDEN = SCHR(1,N)
        ELSE
          HDEN = ZERO
        ENDIF
        IF( HDGL.LE.HDEN ) THEN
          RKLX = 1.D+0
        ELSE
          RKLX = RPLC(1,N)/
     &      (RPLC(1,N) + (((HDGL-HDEN)/SCHR(5,N))**RPLC(2,N)))
        ENDIF
!
!---  Touma and Vauclin relative permeability function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.6 ) THEN
        RKLX = RPLC(1,N)*(ESLZ**RPLC(2,N))
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
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKLX = EXP( RPLC(1,N)*HDGL + RPLC(2,N) )
!
!---  Tabular aqueous relative permeability versus gas saturation
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.10 ) THEN
        ITBX = 0
        RKLX = FNTBLY( SLX,IRLTBL(1,N),IRLTBL(2,N),ITBX )
!
!---  Tabular aqueous relative permeability versus gas saturation
!     w/ cubic spline interpolation  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.11 ) THEN
        RKLX = FSPLNY( SLX,IRLTBL(1,N),IRLTBL(2,N) )
!
!---  Tabular aqueous relative permeability versus ln(capillary head)
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.12 ) THEN
        ITBX = 0
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKLX = FNTBLY( HDGLX,IRLTBL(1,N),IRLTBL(2,N),ITBX )
!
!---  Tabular aqueous relative permeability versus ln(capillary head)
!     w/ cubic spline interpolation  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.13 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKLX = FSPLNY( HDGLX,IRLTBL(1,N),IRLTBL(2,N) )
!
!---  Tabular aqueous relative permeability versus capillary head
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.14 ) THEN
        ITBX = 0
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKLX = FNTBLY( HDGL,IRLTBL(1,N),IRLTBL(2,N),ITBX )
!
!---  Tabular aqueous relative permeability versus capillary head
!     w/ cubic spline interpolation  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.15 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKLX = FSPLNY( HDGL,IRLTBL(1,N),IRLTBL(2,N) )
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
!
!---  Pruess anisotropy permeability function  ---
!
      ELSEIF( IRPL(N).GE.200 .AND. IRPL(N).LT.300 ) THEN
        HDGL_CM = HDGL*1.D+2
        ANISOX = RPLC(5,N)*(RPLC(6,N)**(RPLC(7,N)**HDGL_CM))
        ANISOX = MAX( ANISOX,1.D+0 )
        ANISOX = MIN( ANISOX,RPLC(5,N) )
        RKLX = RKLX*ANISOX
!
!---  Mualem porosity distribution function w van Genuchten 
!     saturation function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.41 ) THEN
        RKLX = SQRT(ESLZ)*((1.D+0-(1.D+0-ESLZ**(1.D+0/RPLC(2,N)))
     &      **RPLC(2,N))**2)
!
!---  Mualem porosity distribution function w Brooks and Corey 
!     saturation function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.42 ) THEN
        RKLX = ESLZ**(2.5D+0 + 2.0D+0/RPLC(2,N))
!
!---  Burdine porosity distribution function w van Genuchten 
!     saturation function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.43 ) THEN
        RKLX = (ESLZ**2)*(1.D+0-(1.D+0-ESLZ**(1.D+0/RPLC(2,N)))
     &      **RPLC(2,N))
!
!---  Burdine porosity distribution function w Brooks and Corey
!     saturation function  ---
!
      ELSEIF( MOD( IRPL(N),100 ).EQ.44 ) THEN
        RKLX = ESLZ**(3.0D+0 + 2.0D+0/RPLC(2,N))
      ENDIF
      RKLX = MIN( MAX( RKLX,1.D-24 ),1.D+0 )
!
!---  Reset subroutine string sequence  ---
!
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
!     STOMPX-CO2
!
!     Compute the maximum relative residuals
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 November 2021.
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
      USE FILES
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
      REAL*8, DIMENSION(LUK*(1+LWELL+LSPILL)) :: RSDLX	
      REAL*8, DIMENSION(5) :: VARX	
      INTEGER, DIMENSION(LUK*(1+LWELL+LSPILL)) :: NSDLX,NPHLX,NPHX
      INTEGER, DIMENSION(LUK*(1+LWELL+LSPILL)) :: IDLX,IDX
      INTEGER, DIMENSION(NFCGC(ID+1)) :: IRSDX
      INTEGER, DIMENSION(7) :: IVARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER*64 PH_CND(10)
!
!----------------------Data Statements---------------------------------!
!
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






        MPW = NMD + IEQW
        DPW = BLU(MPW)
        MPA = NMD + IEQA
        DPA = BLU(MPA)
        N_DB = ND(N)
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          MPS = NMD + IEQS
          DPS = BLU(MPS)
        ELSE
          DPS = 0.D+0
        ENDIF
!
!---    Increment equation counter for next active node  ---
!
        NMD = NMD + ISVC

!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE

!
!---    Skip selected nodes in the residual calculation
!       (not implemented)  ---
!
!        IF( ISKP(N).EQ.1 ) CYCLE
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
          RSDX = MIN( ABS(BLU(MPW))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSDLX(IEQW) ) THEN
            RSDLX(IEQW) = RSDX
            NSDLX(IEQW) = N
            NPHLX(IEQW) = NPHAZ(2,N)
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
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
            RSDX = MIN( ABS(BLU(MPA))/
     &        MAX( (PG(2,N)+PATM)/HCX,PATM/HCX ),
     &        ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSDLX(IEQA) ) THEN
              RSDLX(IEQA) = RSDX
              NSDLX(IEQA) = N
              NPHLX(IEQA) = NPHAZ(2,N)
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
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
            IF( RSDX.GT.RSDLX(IEQS) ) THEN
              RSDLX(IEQS) = RSDX
              NSDLX(IEQS) = N
              NPHLX(IEQS) = NPHAZ(2,N)
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
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
          RSDX = MIN( ABS(BLU(MPW))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSDLX(IEQW) ) THEN
            RSDLX(IEQW) = RSDX
            NSDLX(IEQW) = N
            NPHLX(IEQW) = NPHAZ(2,N)
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      CO2 mass equation  ---
!
          IF( SG(2,N).GT.1.D-3 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLA(2,N))*DTI*VOL(N)
            RSDX = MIN( ABS(BLU(MPA))/(ABS(PG(2,N))+PATM),
     &        ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSDLX(IEQA) ) THEN
              RSDLX(IEQA) = RSDX
              NSDLX(IEQA) = N
              NPHLX(IEQA) = NPHAZ(2,N)
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
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
            IF( RSDX.GT.RSDLX(IEQS) ) THEN
              RSDLX(IEQS) = RSDX
              NSDLX(IEQS) = N
              NPHLX(IEQS) = NPHAZ(2,N)
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
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
          RSDX = MIN( ABS(BLU(MPW))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSDLX(IEQW) ) THEN
            RSDLX(IEQW) = RSDX
            NSDLX(IEQW) = N
            NPHLX(IEQW) = NPHAZ(2,N)
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      CO2 mass equation  ---
!
          IF( SG(2,N).GT.1.D-5 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLA(2,N))*DTI*VOL(N)
            RSDX = MIN( ABS(BLU(MPA)/1.D+1),
     &        ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSDLX(IEQA) ) THEN
              RSDLX(IEQA) = RSDX
              NSDLX(IEQA) = N
              NPHLX(IEQA) = NPHAZ(2,N)
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
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
            IF( RSDX.GT.RSDLX(IEQS) ) THEN
              RSDLX(IEQS) = RSDX
              NSDLX(IEQS) = N
              NPHLX(IEQS) = NPHAZ(2,N)
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---    Fully unsaturated system w/ or w/o entrapped gas
!       Water mass - aqueous saturation
!       CO2 mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.8 .OR. NPHAZ(2,N).EQ.9
     &    .OR. NPHAZ(2,N).EQ.10 ) THEN
!
!---      Water mass equation  ---
!
          PGX = PG(2,N) + PATM
          CALL SP_B( T(2,N),XLS(2,N),PSBX )
          ACP = PORD(2,N)*RHOG(2,N)*SG(2,N)*XGW(2,N)*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPW))/PSBX,
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSDLX(IEQW) ) THEN
            RSDLX(IEQW) = RSDX
            NSDLX(IEQW) = N
            NPHLX(IEQW) = NPHAZ(2,N)
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      CO2 mass equation  ---
!
          ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &      RHOL(2,N)*SL(2,N)*XLA(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPA))/ABS(PGX),
     &      ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSDLX(IEQA) ) THEN
            RSDLX(IEQA) = RSDX
            NSDLX(IEQA) = N
            NPHLX(IEQA) = NPHAZ(2,N)
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      Salt mass equation, isobrine option  ---

          IF( ISLC(32).EQ.0 ) THEN
            CALL DENS_S( T(2,N),PGX,RHOSPX )
            ACP = TMS(2,N)*DTI*VOL(N)
            RSDX = MIN( (ABS(BLU(MPS))/RHOSPX),
     &        ABS(RSDL(IEQS,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSDLX(IEQS) ) THEN
              RSDLX(IEQS) = RSDX
              NSDLX(IEQS) = N
              NPHLX(IEQS) = NPHAZ(2,N)
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
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
        IF( ID.EQ.IDLX(M) .AND. NSDLX(M).GT.0 ) THEN
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
!
!---      Water mass equation  ---
!
          NX = NSD(IEQW)
          IF( NX.GT.0 ) THEN
            NPX = NPHX(IEQW)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            PRINT *,'  Water Mass Equation Maximum Residual = ',
     &        RSD(IEQW),': Node = ',NX,': Phase Condition = ',
     &        PH_CND(NPX)(1:NCHX)
          ENDIF
!
!---      CO2 mass equation  ---
!
          NX = NSD(IEQA)
          IF( NX.GT.0 ) THEN
            NPX = NPHX(IEQA)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            PRINT *,'  CO2 Mass Equation Maximum Residual = ',
     &        RSD(IEQA),': Node = ',NX,': Phase Condition = ',
     &        PH_CND(NPX)(1:NCHX)
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            NX = NSD(IEQS)
            IF( NX.GT.0 ) THEN
              NPX = NPHX(IEQS)
              NCHX = INDEX( PH_CND(NPX),'  ') - 1
              PRINT *,'  Salt Equation Maximum Residual = ',
     &          RSD(IEQS),': Node = ',NX,': Phase Condition = ',
     &          PH_CND(NPX)(1:NCHX)
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
          NCH = INDEX( UNTM(1:),'  ' ) - 1
          IF( ID.EQ.0 ) PRINT *,'  Time Step Reduced From ',
     &      DTX*CNVTM,UNTM(1:NCH),' to ',DT*CNVTM,UNTM(1:NCH)
          DO N = 1,NFCGC(ID+1)
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
          IVARX(1) = -1
          VARX(4) = DTX*CNVTM
          VARX(5) = DT*CNVTM
!
!---      Coupled-well pressure  ---
!
          IF( L_CW.EQ.1 ) THEN
            DO NCW = 1,N_CW
              P_CW(2,NCW) = P_CW(1,NCW)
            ENDDO
          ENDIF
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
     &      MPI_INTEGER,STATUS,IERR )
          OFFSET = OFFSET + NVAR*NBYTI
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
          RETURN
        ENDIF
!
!---    Write a convergence failure index in the NSTEP location
!       plus write the global node numbers and phase condition indices
!       for the location of maximum residuals for the water, CO2,
!       and salt equations to output.bin  ---
!
        NVAR = 7
        IVARX(2) = NSD(IEQW)
        IVARX(3) = NPHX(IEQW)
        IVARX(4) = NSD(IEQA)
        IVARX(5) = NPHX(IEQA)
        IF( ISLC(32).EQ.0 ) THEN
          IVARX(6) = NSD(IEQS)
          IVARX(7) = NPHX(IEQS)
        ELSE
          IVARX(6) = 0
          IVARX(7) = 0
        ENDIF
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &    MPI_INTEGER,STATUS,IERR )
        OFFSET = OFFSET + NVAR*NBYTI
        IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
!
!---    Write maximum residuals for the water, CO2, and salt
!       equations and time step reductions to output.bin  ---
!
        NVAR = 5
        VARX(1) = RSD(IEQW)
        VARX(2) = RSD(IEQA)
        IF( ISLC(32).EQ.0 ) THEN
          VARX(3) = RSD(IEQS)
        ELSE
          VARX(3) = 0.D+0
        ENDIF
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,VARX,NVAR,
     &    MPI_REAL8,STATUS,IERR )
        OFFSET = OFFSET + NVAR*NBYTR
        IOFFSET_REF = IOFFSET_REF + NVAR*NBYTR
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDL_CO2 group
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
!     STOMPX-CO2
!
!     Compute source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 November 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE REACT
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
      SUB_LOG(ISUB_LOG) = '/SORC_CO2'
!
!---  Load CO2 sources for linked aqueous CO2   ---
!
      IF( ISPLK(6).NE.0 ) THEN
!
!---    Load CO2 sources associated with chemical
!       reactions  ---
!
        DO NS = 1,NSR(ID+1)
          N = ISRN(NS)
          SRCAX = SRCA(1,N)*DTI
          DO M = 1,ISVF
            SRCA(M,N) = SRCAX
          ENDDO
        ENDDO
!
!---    Zero source terms  ---
!
        DO NS = 1,NSR(ID+1)
          N = ISRN(NS)
          DO M = 1,ISVF
            SRCW(M,N) = 0.D+0
            SRCS(M,N) = 0.D+0
          ENDDO
        ENDDO
      ELSE
!
!---    Zero source terms  ---
!
        DO NS = 1,NSR(ID+1)
          N = ISRN(NS)
          DO M = 1,ISVF
            SRCA(M,N) = 0.D+0
            SRCW(M,N) = 0.D+0
            SRCS(M,N) = 0.D+0
          ENDDO
        ENDDO
      ENDIF
!
!---  Loop over sources  ---
!
      DO NS = 1,NSR(ID+1)
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = ISRIN(NS)
        IF( TMZ.LE.SRC(1,1,MB) ) CYCLE
        IF( ISRM(NS).EQ.1 ) THEN
          DO N = 1,8
            SRX(N) = SRC(N,1,MB)
          ENDDO
        ELSE
          IFIND = 0
          DO M = 2,ISRM(NS)
            IF( TMZ.LE.SRC(1,M,MB) ) THEN
             DTSR = MIN( SRC(1,M,MB)-TMZ,DT )
             TFSR = (TMZ-0.5D+0*DTSR-SRC(1,M-1,MB))/
     &         (SRC(1,M,MB)-SRC(1,M-1,MB))
             DO N = 1,8
               SRX(N) = SRC(N,M-1,MB) + TFSR*(SRC(N,M,MB)-SRC(N,M-1,MB))
             ENDDO
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
          PGX = PG(M,N) + PATM
          PLX = PL(M,N) + PATM
          TX = T(M,N)
!
!---      Aqueous volumetric rate w/ aqueous-salt
!         and aqueous-CO2; limit aqueous components
!         to their solubility limit  ---
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
!---          Aqueous-salt concentration  ---
!
              IF( MOD(ISRT(NS),1000)/100.EQ.1 ) THEN
!
!---            Aqueous-CO2 concentration  ---
!
                IF( ISRT(NS)/1000.EQ.1 ) THEN
                  RHOLSX = SRX(8)
                  RHOLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Conc., Salt Aqu. Conc. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_11( TX,PX,PGX,RHOLSX,RHOLAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 relative saturation  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                  RHOLSX = SRX(8)
                  PHILAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Rel. Sat., Salt Aqu. Conc. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_12( TX,PX,PGX,RHOLSX,PHILAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 mass fraction  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                  RHOLSX = SRX(8)
                  XLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Conc. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_13( TX,PX,PGX,RHOLSX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            No aqueous-CO2  ---
!
                ELSE
                  RHOLSX = SRX(8)
                  XLAX = 0.D+0
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Conc. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_13( TX,PX,PGX,RHOLSX,
     &              YLSX,XLSX,XLAX,CHMSGX )
                ENDIF
!
!---          Aqueous-salt relative saturation  ---
!
              ELSEIF( MOD(ISRT(NS),1000)/100.EQ.2 ) THEN
!
!---            Aqueous-CO2 concentration  ---
!
                IF( ISRT(NS)/1000.EQ.1 ) THEN
                  PHILSX = SRX(8)
                  RHOLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'Salt Aqu. Rel. Conc., CO2 Aqu. Conc.  = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_21( TX,PX,PGX,PHILSX,RHOLAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 relative saturation  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                  PHILSX = SRX(8)
                  PHILAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'Salt Aqu. Rel. Sat., CO2 Aqu. Rel. Sat. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_22( TX,PX,PGX,PHILSX,PHILAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 mass fraction  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                  PHILSX = SRX(8)
                  XLAX = SRX(5)
                  CHMSGX(1) = 'null'
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_23( TX,PX,PGX,PHILSX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            No aqueous-CO2  ---
!
                ELSE
                  PHILSX = SRX(8)
                  XLAX = 0.D+0
                  CHMSGX(1) = 'null'
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_23( TX,PX,PGX,PHILSX,
     &              YLSX,XLSX,XLAX,CHMSGX )
                ENDIF
!
!---          Aqueous-salt mass fraction  ---
!
              ELSEIF( MOD(ISRT(NS),1000)/100.EQ.3 ) THEN
!
!---            Aqueous-CO2 concentration  ---
!
                IF( ISRT(NS)/1000.EQ.1 ) THEN
                  YLSX = SRX(8)
                  RHOLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Conc., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_31( TX,PX,PGX,RHOLAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 relative saturation  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                  YLSX = SRX(8)
                  PHILAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Rel. Sat., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_32( TX,PX,PGX,PHILAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 mass fraction  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                  YLSX = SRX(8)
                  XLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---            No aqueous-CO2  ---
!
                ELSE
                  YLSX = SRX(8)
                  XLAX = 0.D+0
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                ENDIF
!
!---          Aqueous-salt molality  ---
!
              ELSEIF( MOD(ISRT(NS),1000)/100.EQ.4 ) THEN
!
!---            Aqueous-CO2 concentration  ---
!
                IF( ISRT(NS)/1000.EQ.1 ) THEN
                  YLSX = SRX(8)
                  RHOLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Conc., Salt Aqu. Molality = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_41( TX,PX,PGX,RHOLAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 relative saturation  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                  YLSX = SRX(8)
                  PHILAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Rel. Sat., Salt Aqu. Molality = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_42( TX,PX,PGX,PHILAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 mass fraction  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                  YLSX = SRX(8)
                  XLAX = SRX(5)
                  CHMSGX(1) = 'null'
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_43( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---            No aqueous-CO2  ---
!
                ELSE
                  YLSX = SRX(8)
                  XLAX = 0.D+0
                  CHMSGX(1) = 'null'
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_43( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                ENDIF
!
!---          No aqueous-salt  ---
!
              ELSE
!
!---            Aqueous-CO2 concentration  ---
!
                IF( ISRT(NS)/1000.EQ.1 ) THEN
                  YLSX = 0.D+0
                  RHOLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Conc., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_31( TX,PX,PGX,RHOLAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 relative saturation  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                  YLSX = 0.D+0
                  PHILAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Rel. Sat., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_32( TX,PX,PGX,PHILAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 mass fraction  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                  YLSX = 0.D+0
                  XLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---            No aqueous-CO2  ---
!
                ELSE
                  YLSX = 0.D+0
                  XLAX = 0.D+0
                  CHMSGX(1) = 'Unconverged Source: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Aqu. Vol. Source: Transition ' //
     &              'Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                ENDIF
              ENDIF
!
!---          Aqueous component fractions and density  ---
!
              XLWX = 1.D+0 - XLSX - XLAX
              CALL DENS_B( TX,PX,XLSX,RHOBX )
              CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
              SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOLX*XLAX
              SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOLX*XLWX
              SRCS(M,N) = SRCS(M,N) + SRX(4)*RHOLX*XLSX
            ENDIF
!
!---      Gas volumetric rate w/ component relative humidity ---
!
          ELSEIF( MOD(ISRT(NS),100).EQ.4 ) THEN
            IF( SRX(4).GE.0.D+0 ) THEN
              PX = MAX( PGX,PLX,SRX(3) )
              CALL SP_W( TX,PSWX )
              PVWX = PSWX*SRX(5)
              PVAX = MAX( PX-PVWX,0.D+0 )
!
!---          Gas density and component fractions  ---
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
!---      Gas volumetric rate w/ component mass fractions ---
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
!---          Gas density and component fractions  ---
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
!---      Salt source  ---
!
          ELSEIF( MOD(ISRT(NS),100).EQ.12 ) THEN
            SRCS(M,N) = SRCS(M,N) + SRX(4)
!
!---      Aqueous mass rate w/ aqueous-CO2 relative
!         saturation, limit aqueous salt
!         sources to solubility limit  ---
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
!---          Aqueous-salt concentration  ---
!
              IF( MOD(ISRT(NS),1000)/100.EQ.1 ) THEN
!
!---            Aqueous-CO2 concentration  ---
!
                IF( ISRT(NS)/1000.EQ.1 ) THEN
                  RHOLSX = SRX(8)
                  RHOLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Conc., Salt Aqu. Conc. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_11( TX,PX,PGX,RHOLSX,RHOLAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 relative saturation  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                  RHOLSX = SRX(8)
                  PHILAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Rel. Sat., Salt Aqu. Conc. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_12( TX,PX,PGX,RHOLSX,PHILAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 mass fraction  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                  RHOLSX = SRX(8)
                  XLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Conc. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_13( TX,PX,PGX,RHOLSX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            No aqueous-CO2  ---
!
                ELSE
                  RHOLSX = SRX(8)
                  XLAX = 0.D+0
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Conc. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_13( TX,PX,PGX,RHOLSX,
     &              YLSX,XLSX,XLAX,CHMSGX )
                ENDIF
!
!---          Aqueous-salt relative saturation  ---
!
              ELSEIF( MOD(ISRT(NS),1000)/100.EQ.2 ) THEN
!
!---            Aqueous-CO2 concentration  ---
!
                IF( ISRT(NS)/1000.EQ.1 ) THEN
                  PHILSX = SRX(8)
                  RHOLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'Salt Aqu. Rel. Conc., CO2 Aqu. Conc.  = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_21( TX,PX,PGX,PHILSX,RHOLAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 relative saturation  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                  PHILSX = SRX(8)
                  PHILAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'Salt Aqu. Rel. Sat., CO2 Aqu. Rel. Sat. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_22( TX,PX,PGX,PHILSX,PHILAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 mass fraction  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                  PHILSX = SRX(8)
                  XLAX = SRX(5)
                  CHMSGX(1) = 'null'
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_23( TX,PX,PGX,PHILSX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            No aqueous-CO2  ---
!
                ELSE
                  PHILSX = SRX(8)
                  XLAX = 0.D+0
                  CHMSGX(1) = 'null'
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_23( TX,PX,PGX,PHILSX,
     &              YLSX,XLSX,XLAX,CHMSGX )
                ENDIF
!
!---          Aqueous-salt mass fraction  ---
!
              ELSEIF( MOD(ISRT(NS),1000)/100.EQ.3 ) THEN
!
!---            Aqueous-CO2 concentration  ---
!
                IF( ISRT(NS)/1000.EQ.1 ) THEN
                  YLSX = SRX(8)
                  RHOLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Conc., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_31( TX,PX,PGX,RHOLAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 relative saturation  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                  YLSX = SRX(8)
                  PHILAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Rel. Sat., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_32( TX,PX,PGX,PHILAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 mass fraction  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                  YLSX = SRX(8)
                  XLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---            No aqueous-CO2  ---
!
                ELSE
                  YLSX = SRX(8)
                  XLAX = 0.D+0
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                ENDIF
!
!---          Aqueous-salt mass fraction  ---
!
              ELSEIF( MOD(ISRT(NS),1000)/100.EQ.4 ) THEN
!
!---            Aqueous-CO2 concentration  ---
!
                IF( ISRT(NS)/1000.EQ.1 ) THEN
                  YLSX = SRX(8)
                  RHOLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Conc., Salt Aqu. Molality = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_41( TX,PX,PGX,RHOLAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 relative saturation  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                  YLSX = SRX(8)
                  PHILAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Rel. Sat., Salt Aqu. Molality = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_42( TX,PX,PGX,PHILAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 mass fraction  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                  YLSX = SRX(8)
                  XLAX = SRX(5)
                  CHMSGX(1) = 'null'
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_43( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---            No aqueous-CO2  ---
!
                ELSE
                  YLSX = SRX(8)
                  XLAX = 0.D+0
                  CHMSGX(1) = 'null'
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_43( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                ENDIF
!
!---          No aqueous-salt  ---
!
              ELSE
!
!---            Aqueous-CO2 concentration  ---
!
                IF( ISRT(NS)/1000.EQ.1 ) THEN
                  YLSX = 0.D+0
                  RHOLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Conc., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_31( TX,PX,PGX,RHOLAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 relative saturation  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.2 ) THEN
                  YLSX = 0.D+0
                  PHILAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Rel. Sat., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_32( TX,PX,PGX,PHILAX,
     &              YLSX,XLSX,XLAX,CHMSGX )
!
!---            Aqueous-CO2 mass fraction  ---
!
                ELSEIF( ISRT(NS)/1000.EQ.3 ) THEN
                  YLSX = 0.D+0
                  XLAX = SRX(5)
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
!
!---            No aqueous-CO2  ---
!
                ELSE
                  YLSX = 0.D+0
                  XLAX = 0.D+0
                  CHMSGX(1) = 'Unconverged Initial Conditions: ' //
     &              'CO2 Aqu. Mass Frac., Salt Aqu. Mass Frac. = '
                  CHMSGX(2) = 'Initial Condition Transition: ' //
     &              'Vapor Pressure > Gas Pressure: PVB + PVA = '
                  CALL FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
                ENDIF
              ENDIF
!
!---          CO2 vapor pressure ---
!
              CALL SP_B( TX,XLSX,PSBX )
              PVBX = PSBX
              CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &          XGAX,XGWX,XLAZ,XLSX,XLWX,XMGAX,XMGWX,XMLAX,
     &          XMLSX,XMLWX )
              PVAX = MIN( XLAX/XLAZ,1.D+0 )*PGAX
!
!---          Aqueous component fraction  ---
!
              XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
              SRCA(M,N) = SRCA(M,N) + SRX(4)*XLAX
              SRCW(M,N) = SRCW(M,N) + SRX(4)*XLWX
              SRCS(M,N) = SRCS(M,N) + SRX(4)*XLSX
            ENDIF
!
!---      Gas mass rate w/ component relative humidity ---
!
          ELSEIF( MOD(ISRT(NS),100).EQ.8 ) THEN
            IF( SRX(4).GE.0.D+0 ) THEN
              PX = MAX( PGX,PLX,SRX(3) )
              CALL SP_W( TX,PSWX )
              PVWX = PSWX*SRX(5)
              PVAX = MAX( PX-PVWX,0.D+0 )
!
!---          Gas density and component fractions  ---
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
!---      Gas mass rate w/ component mass fractions ---
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
!---          Gas density and component fractions  ---
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
!---      Salt density source  ---
!
          ELSEIF( MOD(ISRT(NS),100).EQ.11 ) THEN
            SRCS(M,N) = SRCS(M,N) + SRX(4)*VOL(N)
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     STOMPX-CO2
!
!     Compute coupled-equation source integrals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 17 December 2021
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
      REAL*8 SRX(8+LSOLU)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORIC_CO2'
!
!---  Integrate coupled-well rates  ---
!
      IF( L_CW.EQ.1 ) THEN
        DO NCW = 1,N_CW
          QM_CW(2,NCW) = QM_CW(2,NCW) + QM_CW(1,NCW)*DT
          QM_CW(4,NCW) = QM_CW(4,NCW) + QM_CW(3,NCW)*DT
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORIC_CO2 group  ---
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
!     STOMPX-CO2
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
      SUB_LOG(ISUB_LOG) = '/SP_CO2'
!
!---  van Genuchten saturation function w/o gas entrapment,
!     w/ or w/o Webb extension  ---
!
      IF( ISCHR(N).EQ.1 ) THEN
!
!---    van Genuchten saturation function w/o gas entrapment,
!       w/ Webb extension ASL = SL  ---
!
        IF( ISM(N).EQ.2 ) THEN
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(9,N)
          SLRX = SCHR(4,N)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(8,N)
            HDGL = MIN( HDGL,SCHR(12,N) )
            DMPX = SMPX/(LOG10(SCHR(12,N))-LOG10(HMPX))
            SLX = -(LOG10(HDGL)-LOG10(SCHR(12,N)))*DMPX
            ASLX = SLX
!
!---      Capillary head at or below the matching point head,
!         use van Genuchten function
!
          ELSE
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
          SLRX = SCHR(4,N)
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
      ELSEIF( ISCHR(N).EQ.2 ) THEN
!
!---    Brooks and Corey saturation function w/o gas entrapment,
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(N).EQ.2 ) THEN
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(9,N)
          SLRX = SCHR(4,N)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(8,N)
            HDGL = MIN( HDGL,SCHR(12,N) )
            DMPX = SMPX/(LOG10(SCHR(12,N))-LOG10(HMPX))
            SLX = -(LOG10(HDGL)-LOG10(SCHR(12,N)))*DMPX
            ASLX = SLX
!
!---      Capillary head at or below the matching point head,
!         use Brooks-Corey function
!
          ELSE
            CL = MAX( SCHR(3,N),SMALL )
            IF( HDGL.LE.SCHR(1,N) ) THEN
              ASLX = 1.D+0
            ELSE
              ASLX = (SCHR(1,N)/HDGL)**CL
            ENDIF
            SLRX = SCHR(4,N)
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
          CL = MAX( SCHR(3,N),SMALL )
          IF( HDGL.LE.SCHR(1,N) ) THEN
            ASLX = 1.D+0
          ELSE
            ASLX = (SCHR(1,N)/HDGL)**CL
          ENDIF
          SLRX = SCHR(4,N)
          SGRX = SCHR(15,N)
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
      ELSEIF( ISCHR(N).EQ.3 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
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
        ASL_M = (1.D+0/(1.D+0 + (SCHR(1,N)*HDGL)**CN))**CM
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,N)),ZERO )*REALX
        SLR_M = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        SLMX = ASL_M*(1.D+0-SLR_M) + SLR_M
        CN = MAX( SCHR(6,N),SMALL )
        IF( SCHR(15,N).LE.EPSL ) THEN
          IF( MOD( IRPL(N),100 ).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(15,N)
        ENDIF
        ASL_F = (1.D+0/(1.D+0 + (SCHR(5,N)*HDGL)**CN))**CM
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,N)),ZERO )*REALX
        SLR_F = MAX( (1.D+0-HSCL)*SCHR(7,N),ZERO )
        SLFX = ASL_F*(1.D+0-SLR_F) + SLR_F
        PORD_MX = (1.D+0-POR(4,N))*POR0(2,N)/
     &    ( POR(4,N) + (1.D+0-POR(4,N))*POR0(2,N) + SMALL )
        PORD_FX = POR(4,N)/
     &    ( POR(4,N) + (1.D+0-POR(4,N))*POR0(2,N) + SMALL )
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
      ELSEIF( ISCHR(N).EQ.4 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        CL = MAX( SCHR(3,N),SMALL )
        IF( HDGL.LE.SCHR(1,N) ) THEN
          ASL_M = 1.D+0
        ELSE
          ASL_M = (SCHR(1,N)/HDGL)**CL
        ENDIF
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,N)),ZERO )*REALX
        SLR_M = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        SLMX = ASL_M*(1.D+0-SLR_M) + SLR_M
        CL = MAX( SCHR(6,N),SMALL )
        IF( HDGL.LE.SCHR(5,N) ) THEN
          ASL_F = 1.D+0
        ELSE
          ASL_F = (SCHR(5,N)/HDGL)**CL
        ENDIF
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,N)),ZERO )*REALX
        SLR_F = MAX( (1.D+0-HSCL)*SCHR(7,N),ZERO )
        SLFX = ASL_F*(1.D+0-SLR_F) + SLR_F
        PORD_MX = (1.D+0-POR(4,N))*POR0(2,N)/
     &    ( POR(4,N) + (1.D+0-POR(4,N))*POR0(2,N) + SMALL )
        PORD_FX = POR(4,N)/
     &    ( POR(4,N) + (1.D+0-POR(4,N))*POR0(2,N) + SMALL )
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
      ELSEIF( ISCHR(N).EQ.5 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( HDGL.LE.SCHR(1,N) ) THEN
          ASLX = 1.D+0
        ELSE
          ASLX = SCHR(2,N)/(SCHR(2,N)
     &      + (((HDGL-SCHR(1,N))/SCHR(5,N))**SCHR(3,N)))
        ENDIF
        REALX = REAL(ISM(N))
        HSCL = MAX( LOG(HDGL)/LOG(SCHR(12,N)),ZERO )*REALX
        SLRX = MAX( (1.D+0-HSCL)*SCHR(4,N),ZERO )
        SLX = ASLX*(1.D+0-SLRX) + SLRX
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        ESGTX = 0.D+0
        SGTX = 0.D+0
!
!---  Linear or linear-log interpolation function  ---
!
      ELSEIF( ISCHR(N).EQ.10 .OR. ISCHR(N).EQ.11 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( ISCHR(N).EQ.11 ) HDGL = LOG(HDGL)
        ITBX = 0
        ASLX = FNTBLY( HDGL,ISLTBL(1,N),ISLTBL(2,N),ITBX )
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
      ELSEIF( ISCHR(N).EQ.12 .OR. ISCHR(N).EQ.13 ) THEN
        HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( ISCHR(N).EQ.13 ) HDGL = LOG(HDGL)
        ITBX = 0
        SLDX = FNTBLY( HDGL,ISLTBL(1,N),ISLTBL(2,N),ITBX )
        ITBX = 0
        SLIX = FNTBLY( HDGL,ISLTBL(3,N),ISLTBL(4,N),ITBX )
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
!
!---  van Genuchten saturation function w/ gas entrapment,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(N).EQ.101 ) THEN
!
!---    van Genuchten saturation function w/ gas entrapment,
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(N).EQ.2 ) THEN
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(9,N)
          SLRX = SCHR(4,N)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(8,N)
            HDGL = MIN( HDGL,SCHR(12,N) )
            DMPX = SMPX/(LOG10(SCHR(12,N))-LOG10(HMPX))
            SLX = -(LOG10(HDGL)-LOG10(SCHR(12,N)))*DMPX
            ASLX = SLX
!
!---      Capillary head at or below the matching point head,
!         use van Genuchten function
!
          ELSE
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
            PRINT *,'Unrecognized Trapped-Gas Option, ID = ',ID
            ICNV = 4
!
!---        Reset subroutine string sequence  ---
!
            ISUB_LOG = ISUB_LOG-1
            RETURN
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
          SLRX = SCHR(4,N)
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
            PRINT *,'Unrecognized Trapped-Gas Option, ID = ',ID
            ICNV = 4
!
!---        Reset subroutine string sequence  ---
!
            ISUB_LOG = ISUB_LOG-1
            RETURN
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
      ELSEIF( ISCHR(N).EQ.102 ) THEN
!
!---    Brooks and Corey saturation function w/ gas entrapment,
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(N).EQ.2 ) THEN
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(9,N)
          SLRX = SCHR(4,N)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(8,N)
            HDGL = MIN( HDGL,SCHR(12,N) )
            DMPX = SMPX/(LOG10(SCHR(12,N))-LOG10(HMPX))
            SLX = -(LOG10(HDGL)-LOG10(SCHR(12,N)))*DMPX
            ASLX = SLX
!
!---      Capillary head at or below the matching point head,
!         use van Genuchten function
!
          ELSE
            CL = MAX( SCHR(3,N),SMALL )
            IF( HDGL.LE.SCHR(1,N) ) THEN
              ASLX = 1.D+0
            ELSE
              ASLX = (SCHR(1,N)/HDGL)**CL
            ENDIF
            SLRX = SCHR(4,N)
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
            PRINT *,'Unrecognized Trapped-Gas Option, ID = ',ID
            ICNV = 4
!
!---        Reset subroutine string sequence  ---
!
            ISUB_LOG = ISUB_LOG-1
            RETURN
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
          CL = MAX( SCHR(3,N),SMALL )
          IF( HDGL.LE.SCHR(1,N) ) THEN
            ASLX = 1.D+0
          ELSE
            ASLX = (SCHR(1,N)/HDGL)**CL
          ENDIF
          SLRX = SCHR(4,N)
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
            PRINT *,'Unrecognized Trapped-Gas Option, ID = ',ID
            ICNV = 4
!
!---        Reset subroutine string sequence  ---
!
            ISUB_LOG = ISUB_LOG-1
            RETURN
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
      ELSEIF( ISCHR(N).EQ.201 ) THEN
!
!---    van Genuchten saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(N).EQ.2 ) THEN
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
          HMPDX = SCHR(9,N)
          HMPIX = SCHR(11,N)
          SLRX = SCHR(4,N)
          HDGL = MIN( HDGL,SCHR(12,N) )
!
!---      Capillary head above the drainage matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPDX ) THEN
              SMPDX = SCHR(8,N)
              DMPDX = SMPDX/(LOG10(SCHR(12,N))-LOG10(HMPDX))
              SLDX = -(LOG10(HDGL)-LOG10(SCHR(12,N)))*DMPDX
              ASLDX = SLDX
              ASLMX = MAX( MIN( ASLDX,ASLMINX ),0.D+0 )
!
!---        Capillary head above the imbibition matching point head,
!           use Webb extension  ---
!
            IF( HDGL.GT.HMPIX ) THEN
              SMPIX = SCHR(10,N)
              DMPIX = SMPIX/(LOG10(SCHR(12,N))-LOG10(HMPIX))
              SLIX = -(LOG10(HDGL)-LOG10(SCHR(12,N)))*DMPIX
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
              ASLIX = (1.D+0/(1.D+0 + (SCHR(2,N)*HDGL)**CNI))**CMI
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
            ASLDX = (1.D+0/(1.D+0 + (SCHR(1,N)*HDGL)**CND))**CMD
            ASLDX = ASLDX*(1.D+0-SLRX) + SLRX
            ASLMX = MAX( MIN( ASLDX,ASLMINX ),0.D+0 )
!
!---        Capillary head above the imbibition matching point head,
!           use Webb extension  ---
!
            IF( HDGL.GT.HMPIX ) THEN
              SMPIX = SCHR(10,N)
              DMPIX = SMPIX/(LOG10(SCHR(12,N))-LOG10(HMPIX))
              SLIX = -(LOG10(HDGL)-LOG10(SCHR(12,N)))*DMPIX
              ASLIX = SLIX
              IF( ASLDX.GE.ASLIX ) THEN
                ASLX = 5.D-1*(ASLIX + SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ELSE
                ASLX = 5.D-1*(ASLIX - SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ENDIF
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
              ASLIX = (1.D+0/(1.D+0 + (SCHR(2,N)*HDGL)**CNI))**CMI
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
            PRINT *,'Unrecognized Trapped-Gas Option, ID = ',ID
            ICNV = 4
!
!---        Reset subroutine string sequence  ---
!
            ISUB_LOG = ISUB_LOG-1
            RETURN
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
          ASLDX = (1.D+0/(1.D+0 + (SCHR(1,N)*HDGL)**CND))**CMD
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
          ASLIX = (1.D+0/(1.D+0 + (SCHR(2,N)*HDGL)**CNI))**CMI
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
            PRINT *,'Unrecognized Trapped-Gas Option, ID = ',ID
            ICNV = 4
!
!---        Reset subroutine string sequence  ---
!
            ISUB_LOG = ISUB_LOG-1
            RETURN
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
      ELSEIF( ISCHR(N).EQ.202 ) THEN
!
!---    Brooks and Corey drainage-imbibition saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(N).EQ.2 ) THEN
          HDGL = MAX( BTGLX*(PGX-PLX)/RHORL/GRAV,1.D-14 )
          HMPDX = SCHR(9,N)
          HMPIX = SCHR(11,N)
          SLRX = SCHR(4,N)
          HDGL = MIN( HDGL,SCHR(12,N) )
!
!---      Capillary head above the drainage matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPDX ) THEN
            SMPDX = SCHR(8,N)
            DMPDX = SMPDX/(LOG10(SCHR(12,N))-LOG10(HMPDX))
            SLDX = -(LOG10(HDGL)-LOG10(SCHR(12,N)))*DMPDX
            ASLDX = SLDX
            ASLMX = MAX( MIN( ASLDX,ASLMINX ),0.D+0 )
!
!---        Capillary head above the imbibition matching point head,
!           use Webb extension  ---
!
            IF( HDGL.GT.HMPIX ) THEN
              SMPIX = SCHR(10,N)
              DMPIX = SMPIX/(LOG10(SCHR(12,N))-LOG10(HMPIX))
              SLIX = -(LOG10(HDGL)-LOG10(SCHR(12,N)))*DMPIX
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
              CLI = MAX( SCHR(6,N),SMALL )
              IF( HDGL.LE.SCHR(5,N) ) THEN
                ASLIX = 1.D+0
              ELSE
                ASLIX = (SCHR(5,N)/HDGL)**CLI
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
            CLD = MAX( SCHR(3,N),SMALL )
            IF( HDGL.LE.SCHR(1,N) ) THEN
              ASLDX = 1.D+0
            ELSE
              ASLDX = (SCHR(1,N)/HDGL)**CLD
            ENDIF
            ASLDX = ASLDX*(1.D+0-SLRX) + SLRX
            ASLMX = MAX( MIN( ASLDX,ASLMINX ),0.D+0 )
!
!---        Capillary head above the imbibition matching point head,
!           use Webb extension  ---
!
            IF( HDGL.GT.HMPIX ) THEN
              SMPIX = SCHR(10,N)
              DMPIX = SMPIX/(LOG10(SCHR(12,N))-LOG10(HMPIX))
              SLIX = -(LOG10(HDGL)-LOG10(SCHR(12,N)))*DMPIX
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
              CLI = MAX( SCHR(6,N),SMALL )
              IF( HDGL.LE.SCHR(5,N) ) THEN
                ASLIX = 1.D+0
              ELSE
                ASLIX = (SCHR(5,N)/HDGL)**CLI
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
            PRINT *,'Unrecognized Trapped-Gas Option, ID = ',ID
            ICNV = 4
!
!---        Reset subroutine string sequence  ---
!
            ISUB_LOG = ISUB_LOG-1
            RETURN
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
          CLD = MAX( SCHR(3,N),SMALL )
          IF( HDGL.LE.SCHR(1,N) ) THEN
            ASLDX = 1.D+0
          ELSE
            ASLDX = (SCHR(1,N)/HDGL)**CLD
          ENDIF
          CLI = MAX( SCHR(6,N),SMALL )
          IF( HDGL.LE.SCHR(5,N) ) THEN
            ASLIX = 1.D+0
          ELSE
            ASLIX = (SCHR(5,N)/HDGL)**CLI
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
            PRINT *,'Unrecognized Trapped-Gas Option, ID = ',ID
            ICNV = 4
!
!---        Reset subroutine string sequence  ---
!
            ISUB_LOG = ISUB_LOG-1
            RETURN
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
!---  End of SP_CO2 group  ---
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
!     STOMPX-CO2
!
!     Update the primary variables on field cells w/o ghost cells.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 November 2021.
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
      REAL*8, DIMENSION(12) :: VARZ	
      REAL*8, DIMENSION(12) :: VARX	
      INTEGER, DIMENSION(12) :: IVARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
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
     &      MPI_INTEGER,STATUS,IERR )
          OFFSET = OFFSET + NVAR*NBYTI
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
          NVAR = 2
          VARX(1) = DTX*CNVTM
          VARX(2) = DT*CNVTM
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,VARX,NVAR,
     &      MPI_REAL8,STATUS,IERR )
          OFFSET = OFFSET + NVAR*NBYTR
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTR
!
!---      Reset principal variables to old time step values
!         on field and ghost cells  ---
!
          DO N = 1,NFCGC(ID+1)
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
!---      Coupled-well pressure  ---
!
          IF( L_CW.EQ.1 ) THEN
            DO NCW = 1,N_CW
              P_CW(2,NCW) = P_CW(1,NCW)
            ENDDO
          ENDIF
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
!---      Coupled-well pressure  ---
!
          IF( L_CW.EQ.1 ) THEN
            DO NCW = 1,N_CW
              P_CW(2,NCW) = P_CW(1,NCW)
            ENDDO
          ENDIF
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
     &      MPI_INTEGER,STATUS,IERR )
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
!---    Skip for inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE

        IF( IERRL.NE.(NFLD+1) ) CYCLE
        MPW = NMD + IEQW
        DPW = BLU(MPW)
        MPA = NMD + IEQA
        DPA = BLU(MPA)
        N_DB = ND(N)
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          MPS = NMD + IEQS
          DPS = BLU(MPS)
        ELSE
          DPS = 0.D+0
        ENDIF
!
!---    Increment equation counter for next active node  ---
!
        NMD = NMD + ISVC
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
          DPX = 1.D+6
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          PL(2,N) = PL(2,N) + DPW
          PL(2,N) = MIN( PL(2,N),5.D+8 )
!
!---      Zero negative corrections for zero aqueous CO2  ---
!
          IF( XLA(2,N)/EPSL.LT.EPSL .AND. BLU(MPA)/EPSL.LT.EPSL ) THEN
            BLU(MPA) = 0.D+0
            DPA = 0.D+0
          ENDIF
          XLA(2,N) = MAX( (XLA(2,N)+DPA),0.D+0 )
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
          IF( ISCHR(N).EQ.2 .OR. ISCHR(N).EQ.102 .OR.
     &        ISCHR(N).EQ.202 ) THEN
            ENPR = SCHR(1,N)*RHORL*GRAV
          ELSEIF( ISCHR(N).EQ.4 ) THEN
            ENPR = MIN( SCHR(1,N),SCHR(5,N) )*RHORL*GRAV
          ELSE
            ENPR = 0.D+0
          ENDIF
!
!---      Limit changes in pressure, and trap large changes in
!         gas pressure  ---
!
          DPX = MAX( 1.D+6,2.5D-1*(PG(2,N)-PL(2,N)) )
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          DPA = SIGN( MIN(ABS(DPX),ABS(DPA)),DPA )
!
!---      Relax pressure updates when transitioning to unsaturated
!         conditions  ---
!
          IF( (PG(2,N)+DPA)-(PL(2,N)+DPW).LT.ENPR ) THEN
            DPA = 6.D-1*DPA
            DPW = 6.D-1*DPW
          ENDIF
          PG(2,N) = PG(2,N) + DPA
          PG(2,N) = MIN( PG(2,N),5.D+8 )
          PL(2,N) = PL(2,N) + DPW
          PL(2,N) = MAX( PL(2,N),(PG(2,N)-SCHR(12,N)*RHORL*GRAV) )
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
!---      Maintain the gas pressure above or at the water vapor
!         pressure  ---
!
          CALL SP_B( T(2,N),XLS(2,N),PSBX )
          CALL DENS_B( T(2,N),PSBX,XLS(2,N),RHOBX )
          PCX = MAX( PSBX-PL(2,N),0.D+0 )
          CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLS(2,N),N )
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
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          PL(2,N) = PL(2,N) + DPW
          PL(2,N) = MIN( PL(2,N),5.D+8 )
!
!---      Limit changes in trapped gas  ---
!
          DPX = 1.D-1*SCHR(15,N)
          DPA = SIGN( MIN(ABS(DPX),ABS(DPA)),DPA )
          SG(2,N) = SG(2,N) + DPA
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
          DPA = SIGN( MIN(ABS(DPX),ABS(DPA)),DPA )
          PG(2,N) = PG(2,N) + DPA
!
!---      Limit changes in water vapor pressure  ---
!
          IF( PVW(2,N)/EPSL.LT.EPSL .AND. BLU(MPW)/EPSL.LT.EPSL ) THEN
            BLU(MPW) = 0.D+0
            DPW = 0.D+0
          ENDIF
          CALL SOL_LS( T(2,N),XLSMX )
          XLSX = MIN( XLS(2,N),XLSMX )
          CALL SP_B( T(2,N),XLSX,PSBX )
          DPX = 2.5D-1*PSBX
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          PVW(2,N) = PVW(2,N) + DPW
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
          CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX,XLSX,N )
          PG(2,N) = MAX( PG(2,N),(PVBX-PATM) )
        ENDIF

!
!---    Skip checks for excessive pressure or temperature for
!       ghost cells   ---
!
        IF( IGHC(N).EQ.1 ) CYCLE

!
!---    Check for excessive pressure or temperature   ---
!
        PGX = PG(2,N)+PATM
        PLX = PL(2,N)+PATM
        TKX = T(2,N)+TABS
        IF( PGX.GT.8.D+8 .OR. PGX.LT.0.D+0 ) IERR = 1
        IF( PLX.GT.8.D+8 ) IERR = 1
        IF( TKX.GT.TCRW .OR. TKX.LT.TABS ) IERR = 1
!
!---    Excess changes in primary variables, skip updates  ---
!
        IF( IERR.EQ.1 .AND. IERRL.EQ.(NFLD+1) ) THEN
          IERRL = ND(N)
          VARZ(1) = REAL( N )
          VARZ(2) = REAL( NPHAZ(2,N) )
          VARZ(3) = DPW
          VARZ(4) = DPA
          VARZ(5) = DPS
          VARZ(6) = PL(2,N)+PATM
          VARZ(7) = PG(2,N)+PATM
          VARZ(8) = SG(2,N)
          VARZ(9) = XLA(2,N)
          VARZ(10) = PVW(2,N)
          VARZ(11) = YLS(2,N)
          VARZ(12) = TMS(2,N)
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
          CALL MPI_RECV( VARZ,12,MPI_REAL8,IDX,18,MPI_COMM_WORLD,
     &      STATUS,IERR )
          N = INT(VARZ(1))
          NPHZX = MOD(INT(VARZ(2)),100)
          PRINT *,'---  Excessive Primary Variable Change  ---'
          PRINT *,'  Node = ',ND(N)
          PRINT *,'  Phase Condition = ',PH_CND(NPHZX)
          PRINT *,'  DPW = ',VARZ(3)
          PRINT *,'  DPA = ',VARZ(4)
          IF( ISLC(32).EQ.0 ) THEN
            PRINT *,'  DPS = ',VARZ(5)
          ENDIF
          PRINT *,'  Aqueous Pressure, Pa = ',VARZ(6)
          PRINT *,'  Gas Pressure, Pa = ',VARZ(7)
          PRINT *,'  Gas Saturation = ',VARZ(8)
          PRINT *,'  Aqueous-CO2 Mass Fraction = ',VARZ(9)
          PRINT *,'  Water Vapor Pressure = ',VARZ(10)
          PRINT *,'  Total-Salt Aqu. Mass Fraction = ',VARZ(11)
          PRINT *,'  Salt Volumetric Concentration = ',VARZ(12)
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
!---      Write a excessive primary variable index in the NSTEP location
!         plus write the global node numbers and phase condition indices
!         for the location of maximum residuals for the water, CO2,
!         and salt equations to output.bin  ---
!
          NVAR = 3
          IVARX(1) = -4
          IVARX(2) = ND(IERRG)
          IVARX(3) = NPHZX
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &      MPI_INTEGER,STATUS,IERR )
          OFFSET = OFFSET + NVAR*NBYTI
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
!
!---      Write maximum residuals for the water, CO2, and salt
!         equations and time step reductions to output.bin  ---
!
          NVAR = 9
          VARX(1) = VARZ(6)
          VARX(2) = VARZ(7)
          VARX(3) = VARZ(8)
          VARX(4) = VARZ(9)
          VARX(5) = VARZ(10)
          VARX(6) = VARZ(11)
          VARX(7) = VARZ(12)
          VARX(8) = DTX*CNVTM
          VARX(9) = DT*CNVTM
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,VARX,NVAR,
     &      MPI_REAL8,STATUS,IERR )
          OFFSET = OFFSET + NVAR*NBYTR
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTR
!
!---    Number of time step reductions failure: stop simulation  ---
!
        ELSE
          DO N = 1,NFCGC(ID+1)
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
     &      MPI_INTEGER,STATUS,IERR )
          OFFSET = OFFSET + NVAR*NBYTI
          IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_GC_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     STOMPX-CO2
!
!     Update the primary variables on ghost cells
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/UPDT_GC_CO2'
      NPVX = 12
!      IF( ID.EQ.0 ) THEN
!        DO N = 1,NP
!          DO M = 1,6
!            PRINT *,'NCGC(',M,',',N,') = ',NCGC(M,N)
!          ENDDO
!        ENDDO
!      ENDIF
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
          SBFB(NCS+1) = PG(2,NLSGC(M+MCS))
          SBFB(NCS+2) = PL(2,NLSGC(M+MCS))
          SBFB(NCS+3) = PVA(2,NLSGC(M+MCS))
          SBFB(NCS+4) = PVW(2,NLSGC(M+MCS))
          SBFB(NCS+5) = SG(2,NLSGC(M+MCS))
          SBFB(NCS+6) = SGT(2,NLSGC(M+MCS))
          SBFB(NCS+7) = SL(2,NLSGC(M+MCS))
          SBFB(NCS+8) = T(2,NLSGC(M+MCS))
          SBFB(NCS+9) = TMS(2,NLSGC(M+MCS))
          SBFB(NCS+10) = XLA(2,NLSGC(M+MCS))
          SBFB(NCS+11) = XLS(2,NLSGC(M+MCS))
          SBFB(NCS+12) = YLS(2,NLSGC(M+MCS))
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
          PG(2,NLRGC(M+MCR)) = RBFB(NCR+1)
          PL(2,NLRGC(M+MCR)) = RBFB(NCR+2)
          PVA(2,NLRGC(M+MCR)) = RBFB(NCR+3)
          PVW(2,NLRGC(M+MCR)) = RBFB(NCR+4)
          SG(2,NLRGC(M+MCR)) = RBFB(NCR+5)
          SGT(2,NLRGC(M+MCR)) = RBFB(NCR+6)
          SL(2,NLRGC(M+MCR)) = RBFB(NCR+7)
          T(2,NLRGC(M+MCR)) = RBFB(NCR+8)
          TMS(2,NLRGC(M+MCR)) = RBFB(NCR+9)
          XLA(2,NLRGC(M+MCR)) = RBFB(NCR+10)
          XLS(2,NLRGC(M+MCR)) = RBFB(NCR+11)
          YLS(2,NLRGC(M+MCR)) = RBFB(NCR+12)
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
          SBFS(NCS+1) = PG(2,NLSGC(M+MCS))
          SBFS(NCS+2) = PL(2,NLSGC(M+MCS))
          SBFS(NCS+3) = PVA(2,NLSGC(M+MCS))
          SBFS(NCS+4) = PVW(2,NLSGC(M+MCS))
          SBFS(NCS+5) = SG(2,NLSGC(M+MCS))
          SBFS(NCS+6) = SGT(2,NLSGC(M+MCS))
          SBFS(NCS+7) = SL(2,NLSGC(M+MCS))
          SBFS(NCS+8) = T(2,NLSGC(M+MCS))
          SBFS(NCS+9) = TMS(2,NLSGC(M+MCS))
          SBFS(NCS+10) = XLA(2,NLSGC(M+MCS))
          SBFS(NCS+11) = XLS(2,NLSGC(M+MCS))
          SBFS(NCS+12) = YLS(2,NLSGC(M+MCS))
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
          PG(2,NLRGC(M+MCR)) = RBFS(NCR+1)
          PL(2,NLRGC(M+MCR)) = RBFS(NCR+2)
          PVA(2,NLRGC(M+MCR)) = RBFS(NCR+3)
          PVW(2,NLRGC(M+MCR)) = RBFS(NCR+4)
          SG(2,NLRGC(M+MCR)) = RBFS(NCR+5)
          SGT(2,NLRGC(M+MCR)) = RBFS(NCR+6)
          SL(2,NLRGC(M+MCR)) = RBFS(NCR+7)
          T(2,NLRGC(M+MCR)) = RBFS(NCR+8)
          TMS(2,NLRGC(M+MCR)) = RBFS(NCR+9)
          XLA(2,NLRGC(M+MCR)) = RBFS(NCR+10)
          XLS(2,NLRGC(M+MCR)) = RBFS(NCR+11)
          YLS(2,NLRGC(M+MCR)) = RBFS(NCR+12)
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
          SBFW(NCS+1) = PG(2,NLSGC(M+MCS))
          SBFW(NCS+2) = PL(2,NLSGC(M+MCS))
          SBFW(NCS+3) = PVA(2,NLSGC(M+MCS))
          SBFW(NCS+4) = PVW(2,NLSGC(M+MCS))
          SBFW(NCS+5) = SG(2,NLSGC(M+MCS))
          SBFW(NCS+6) = SGT(2,NLSGC(M+MCS))
          SBFW(NCS+7) = SL(2,NLSGC(M+MCS))
          SBFW(NCS+8) = T(2,NLSGC(M+MCS))
          SBFW(NCS+9) = TMS(2,NLSGC(M+MCS))
          SBFW(NCS+10) = XLA(2,NLSGC(M+MCS))
          SBFW(NCS+11) = XLS(2,NLSGC(M+MCS))
          SBFW(NCS+12) = YLS(2,NLSGC(M+MCS))
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
          PG(2,NLRGC(M+MCR)) = RBFW(NCR+1)
          PL(2,NLRGC(M+MCR)) = RBFW(NCR+2)
          PVA(2,NLRGC(M+MCR)) = RBFW(NCR+3)
          PVW(2,NLRGC(M+MCR)) = RBFW(NCR+4)
          SG(2,NLRGC(M+MCR)) = RBFW(NCR+5)
          SGT(2,NLRGC(M+MCR)) = RBFW(NCR+6)
          SL(2,NLRGC(M+MCR)) = RBFW(NCR+7)
          T(2,NLRGC(M+MCR)) = RBFW(NCR+8)
          TMS(2,NLRGC(M+MCR)) = RBFW(NCR+9)
          XLA(2,NLRGC(M+MCR)) = RBFW(NCR+10)
          XLS(2,NLRGC(M+MCR)) = RBFW(NCR+11)
          YLS(2,NLRGC(M+MCR)) = RBFW(NCR+12)
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
          SBFE(NCS+1) = PG(2,NLSGC(M+MCS))
          SBFE(NCS+2) = PL(2,NLSGC(M+MCS))
          SBFE(NCS+3) = PVA(2,NLSGC(M+MCS))
          SBFE(NCS+4) = PVW(2,NLSGC(M+MCS))
          SBFE(NCS+5) = SG(2,NLSGC(M+MCS))
          SBFE(NCS+6) = SGT(2,NLSGC(M+MCS))
          SBFE(NCS+7) = SL(2,NLSGC(M+MCS))
          SBFE(NCS+8) = T(2,NLSGC(M+MCS))
          SBFE(NCS+9) = TMS(2,NLSGC(M+MCS))
          SBFE(NCS+10) = XLA(2,NLSGC(M+MCS))
          SBFE(NCS+11) = XLS(2,NLSGC(M+MCS))
          SBFE(NCS+12) = YLS(2,NLSGC(M+MCS))
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
          PG(2,NLRGC(M+MCR)) = RBFE(NCR+1)
          PL(2,NLRGC(M+MCR)) = RBFE(NCR+2)
          PVA(2,NLRGC(M+MCR)) = RBFE(NCR+3)
          PVW(2,NLRGC(M+MCR)) = RBFE(NCR+4)
          SG(2,NLRGC(M+MCR)) = RBFE(NCR+5)
          SGT(2,NLRGC(M+MCR)) = RBFE(NCR+6)
          SL(2,NLRGC(M+MCR)) = RBFE(NCR+7)
          T(2,NLRGC(M+MCR)) = RBFE(NCR+8)
          TMS(2,NLRGC(M+MCR)) = RBFE(NCR+9)
          XLA(2,NLRGC(M+MCR)) = RBFE(NCR+10)
          XLS(2,NLRGC(M+MCR)) = RBFE(NCR+11)
          YLS(2,NLRGC(M+MCR)) = RBFE(NCR+12)
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
          SBFN(NCS+1) = PG(2,NLSGC(M+MCS))
          SBFN(NCS+2) = PL(2,NLSGC(M+MCS))
          SBFN(NCS+3) = PVA(2,NLSGC(M+MCS))
          SBFN(NCS+4) = PVW(2,NLSGC(M+MCS))
          SBFN(NCS+5) = SG(2,NLSGC(M+MCS))
          SBFN(NCS+6) = SGT(2,NLSGC(M+MCS))
          SBFN(NCS+7) = SL(2,NLSGC(M+MCS))
          SBFN(NCS+8) = T(2,NLSGC(M+MCS))
          SBFN(NCS+9) = TMS(2,NLSGC(M+MCS))
          SBFN(NCS+10) = XLA(2,NLSGC(M+MCS))
          SBFN(NCS+11) = XLS(2,NLSGC(M+MCS))
          SBFN(NCS+12) = YLS(2,NLSGC(M+MCS))
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
          PG(2,NLRGC(M+MCR)) = RBFN(NCR+1)
          PL(2,NLRGC(M+MCR)) = RBFN(NCR+2)
          PVA(2,NLRGC(M+MCR)) = RBFN(NCR+3)
          PVW(2,NLRGC(M+MCR)) = RBFN(NCR+4)
          SG(2,NLRGC(M+MCR)) = RBFN(NCR+5)
          SGT(2,NLRGC(M+MCR)) = RBFN(NCR+6)
          SL(2,NLRGC(M+MCR)) = RBFN(NCR+7)
          T(2,NLRGC(M+MCR)) = RBFN(NCR+8)
          TMS(2,NLRGC(M+MCR)) = RBFN(NCR+9)
          XLA(2,NLRGC(M+MCR)) = RBFN(NCR+10)
          XLS(2,NLRGC(M+MCR)) = RBFN(NCR+11)
          YLS(2,NLRGC(M+MCR)) = RBFN(NCR+12)
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
          SBFT(NCS+1) = PG(2,NLSGC(M+MCS))
          SBFT(NCS+2) = PL(2,NLSGC(M+MCS))
          SBFT(NCS+3) = PVA(2,NLSGC(M+MCS))
          SBFT(NCS+4) = PVW(2,NLSGC(M+MCS))
          SBFT(NCS+5) = SG(2,NLSGC(M+MCS))
          SBFT(NCS+6) = SGT(2,NLSGC(M+MCS))
          SBFT(NCS+7) = SL(2,NLSGC(M+MCS))
          SBFT(NCS+8) = T(2,NLSGC(M+MCS))
          SBFT(NCS+9) = TMS(2,NLSGC(M+MCS))
          SBFT(NCS+10) = XLA(2,NLSGC(M+MCS))
          SBFT(NCS+11) = XLS(2,NLSGC(M+MCS))
          SBFT(NCS+12) = YLS(2,NLSGC(M+MCS))
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
          PG(2,NLRGC(M+MCR)) = RBFT(NCR+1)
          PL(2,NLRGC(M+MCR)) = RBFT(NCR+2)
          PVA(2,NLRGC(M+MCR)) = RBFT(NCR+3)
          PVW(2,NLRGC(M+MCR)) = RBFT(NCR+4)
          SG(2,NLRGC(M+MCR)) = RBFT(NCR+5)
          SGT(2,NLRGC(M+MCR)) = RBFT(NCR+6)
          SL(2,NLRGC(M+MCR)) = RBFT(NCR+7)
          T(2,NLRGC(M+MCR)) = RBFT(NCR+8)
          TMS(2,NLRGC(M+MCR)) = RBFT(NCR+9)
          XLA(2,NLRGC(M+MCR)) = RBFT(NCR+10)
          XLS(2,NLRGC(M+MCR)) = RBFT(NCR+11)
          YLS(2,NLRGC(M+MCR)) = RBFT(NCR+12)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_GC_CO2 group
!
      RETURN
      END




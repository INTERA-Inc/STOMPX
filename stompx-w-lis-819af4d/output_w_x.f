!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REFNOD_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Prints convergence and variable information to the screen and
!     output file.
!
!     1 aqueous pressure (Pa), ' PL '
!     2 gas pressure (Pa), ' PG '
!     4 temperature (C), ' T  '
!     5 phase condition (null) 'PHCN'
!     6 aqueous gauge pressure (Pa), 'GPL '
!     7 gas gauge pressure (Pa), 'GPG '
!     9 apparent aqueous saturation (null), 'ASL '
!     11 aqueous saturation (null), ' SL '
!     12 gas saturation (null), ' SG '
!     15 aqueous moisture content (null), 'MCL '
!     19 effective trapped gas (null), 'ESGT'
!     20 diffusive porosity (null), 'PORD'
!     24 aqueous water mass fraction (null), 'XLW '
!     27 aqueous hydraulic head (m), 'HHL '
!     30 rock/soil type (null), 'RSZN'
!     31 aqueous relative permeability (null), 'RKL '
!     34 aqueous density (kg/m^3), 'RHOL'
!     37 total water mass (kg), 'TMW '
!     40 water mass source integral (kg), 'SRIW'
!     49 aqueous courant (null), 'CRNL'
!     63 Webb matching point head (m), 'WMPH'
!     80 stress-xx (Pa), 'SIGXX'
!     83 aqueous matrix (null), 'DSLM'
!     84 aqueous fracture (null), 'DSLF'
!     87 xnc aqueous volumetric flux (m/s), 'ULNC'
!     88 ync aqueous volumetric flux (m/s), 'VLNC'
!     89 znc aqueous volumetric flux (m/s), 'WLNC'
!     100 node number (null), 'NODE'
!     101 stress-yy (Pa), 'SIGXX'
!     103 aqueous compressibility (1/Pa), CMPL
!     105 trapped gas saturation (null), 'SGT '
!     127 minimum effect aqueous saturation (null), 'ESLM'
!     130 stress-zz (Pa), 'SIGZZ'
!     132 X displacement variance (m^2), 'SIGXX'
!     133 Y displacement variance (m^2), 'SIGYY'
!     134 Z displacement variance (m^2), 'SIGZZ'
!     137 water relative humidity (null), 'RHW '
!     140 water mass source rate (kg/s), 'SRCW'
!     150 Webb matching point saturation (null), 'WMPS'
!     176 aqueous viscosity 'VISL'
!     191 integrated water mass (kg), 'IMW '
!     194 integrated aqueous water mass (kg), 'IMLW'
!     200 reserved to control plot file output
!     214 stress-yz (Pa), 'SIGYZ'
!     224 stress-xz (Pa), 'SIGXZ'
!     226 stress-xy (Pa), 'SIGXY'
!     247 x-direction intrinsic permeability (m^2), ' UK '
!     248 y-direction intrinsic permeability (m^2), ' VK '
!     249 z-direction intrinsic permeability (m^2), ' WK '
!     261 integrated water mass source (kg), 'IMSW'
!     275 MCStress (Pa), 'MCSTR'
!     291 x node position, 'XP'
!     292 y node position, 'YP'
!     293 z node position, 'ZP'
!     299 similitude parameter (m^2/s), 'SIMV'
!     369 strain xx, EPSXX
!     370 strain yy, EPSYY
!     371 strain zz, EPSZZ
!     372 strain yz, EPSYZ
!     373 strain xz, EPSXZ
!     374 strain xy, EPSXY
!     375 x displacement, m, DISPX
!     376 y displacement, m, DISPY
!     377 z displacement, m, DISPZ
!     401 solute volumetric concentration 'C   '
!     401 total_species volumetric concentration, mol/m^3, 'SP  '
!     402 solute aqueous concentration 'CL  '
!     402 total_species aqueous concentration 'SPL  '
!     405 solute aqueous mole fraction 'YL  '
!     408 xnc solute flux 'UC  '
!     409 ync solute flux 'VC  '
!     410 znc solute flux 'WC  '
!     411 solute source 'SRC '
!     411 total_species source 'SPSR'
!     412 exchange total_species volumetric conc., mol/m^3, 'SPX '
!     413 exchange total_species aqueous concentration mol/m^3, 'SPLX'
!     423 solute integrated mass 'ICM'
!     423 total_species integrated mass 'SPIM'
!     426 exchange total_species volumetric conc., mol/m^3, 'SPS '
!     427 exchange total_species aqueous concentration mol/m^3, 'SPLS'
!     430 solute mass concentration, 1/kg soil, 'CM'
!     431 solute activity, Bq, 'CA'
!     432 integrated solute activity, Bq, 'CIA'
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 17 November 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE REACT
      USE PROP
      USE OUTPU
      USE HYST
      USE GRID
      USE GEO_MECH
      USE FLUX
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
      REAL*8, DIMENSION(1:NVREF) :: VARX,VARIX
      INTEGER, DIMENSION(1:7) :: IVARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REFNOD_W'
!
!---  Set the offset point to the integrated offset for all
!     processors  ---
!
      OFFSET = IOFFSET_REF
!
!---  Write header of output.bin  ---
!
      IF( (NSTEP-NRST).EQ.0 ) THEN
!
!---    Write number of reference nodes, number of reference
!       node variables, number of equations, and number of
!       significant digits  ---
!
        NVAR = 7
        IVARX(1) = NREF
        IVARX(2) = NVREF
        IVARX(3) = ISVC
        IVARX(4) = ISGNO
        IVARX(5) = IFLD
        IVARX(6) = JFLD
        IVARX(7) = KFLD
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &     MPI_INTEGER,STATUS,IERR)
        OFFSET = OFFSET + NVAR*NBYTI
        IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
!        IF( ID.EQ.0 ) PRINT *,'NREF = ',NREF,' NVREF = ',NVREF,
!     &    ' ISVC = ',ISVC
!
!---    Write reference node variable header strings  ---
!
        NVAR = NVREF*64
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,CHREF,NVAR,
     &     MPI_CHAR,STATUS,IERR)
        OFFSET = OFFSET + NVAR
        IOFFSET_REF = IOFFSET_REF + NVAR
!        IF( ID.EQ.0 ) THEN
!          DO M = 1,NVREF
!            PRINT *,'CHREF(',M,') = ',CHREF(M)
!          ENDDO
!        ENDIF
!
!---    Write reference node unit header strings  ---
!
        NVAR = NVREF*64
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,UNREF,NVAR,
     &     MPI_CHAR,STATUS,IERR)
        OFFSET = OFFSET + NVAR
        IOFFSET_REF = IOFFSET_REF + NVAR
!        IF( ID.EQ.0 ) THEN
!          DO M = 1,NVREF
!            PRINT *,'UNREF(',M,') = ',UNREF(M)
!          ENDDO
!        ENDIF
!
!---    Write time unit string  ---
!
        NVAR = 64
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,UNTM,NVAR,
     &     MPI_CHAR,STATUS,IERR)
        OFFSET = OFFSET + NVAR
        IOFFSET_REF = IOFFSET_REF + NVAR
!        IF( ID.EQ.0 ) THEN
!          PRINT *,'UNTM = ',UNTM
!        ENDIF
!
!---    Write number of reference nodes, number of reference
!       node variables, number of equations, and number of
!       significant digits  ---
!
        NVAR = NREF
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,NDREF,NVAR,
     &     MPI_INTEGER,STATUS,IERR)
        OFFSET = OFFSET + NVAR*NBYTI
        IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
!        IF( ID.EQ.0 ) THEN
!          DO M = 1,NREF
!            PRINT *,'NDREF(',M,') = ',NDREF(M)
!          ENDDO
!        ENDIF
      ENDIF
!
!---  Write current time step and number of iterations to output.bin  --
!
      NVAR = 2
      IVARX(1) = NSTEP
      IVARX(2) = NITER
      IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &   MPI_INTEGER,STATUS,IERR)
      OFFSET = OFFSET + NVAR*NBYTI
!      PRINT *,'3: OFFSET = ',OFFSET,' ID = ',ID
!      IF( ID.EQ.0 ) PRINT *,'NSTEP = ',NSTEP,' NITER = ',NITER
!
!---  Write current time and timestep to output.bin  --
!
      NVAR = 2
      VARX(1) = TM*CNVTM
      VARX(2) = DT*CNVTM
      IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,VARX,NVAR,
     &   MPI_REAL8,STATUS,IERR)
      OFFSET = OFFSET + NVAR*NBYTR
!      PRINT *,'4: OFFSET = ',OFFSET,' ID = ',ID
!      IF( ID.EQ.0 ) PRINT *,'DT = ',VARX(1),' TM = ',VARX(2)
!
!---  Compute integrated outputs  ---
!
      NVAR = 1
      DO NV = 1,NVREF
        IRNV = IREF(NV)
!
!---    Integrated water mass  ---
!
        IF( IRNV.EQ.191 ) THEN
          VARLX = 0.D+0
!
!---      Local summation, skiping inactive and ghost cells  ---
!
          DO N = 1,NFCGC(ID+1)
            IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
            VARLX = VARLX + PORD(2,N)*VOL(N)*XLW(2,N)*SL(2,N)*RHOL(2,N)
          ENDDO
!
!---      Global summation  ---
!
          CALL MPI_ALLREDUCE( VARLX,VARIX(NV),NVAR,MPI_REAL8,MPI_SUM,
     &      MPI_COMM_WORLD,IERR )
!
!---    Integrated aqueous water mass  ---
!
        ELSEIF( IRNV.EQ.194 ) THEN
          VARLX = 0.D+0
!
!---      Local summation, skiping inactive and ghost cells  ---
!
          DO N = 1,NFCGC(ID+1)
            IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
            VARLX = VARLX + PORD(2,N)*VOL(N)*XLW(2,N)*SL(2,N)*RHOL(2,N)
          ENDDO
!
!---      Global summation  ---
!
          CALL MPI_ALLREDUCE( VARLX,VARIX(NV),NVAR,MPI_REAL8,MPI_SUM,
     &      MPI_COMM_WORLD,IERR )
!
!---    Integrated water mass source  ---
!
        ELSEIF( IRNV.EQ.261 ) THEN
          VARLX = 0.D+0
!
!---      Local summation, skiping inactive and ghost cells  ---
!
          DO N = 1,NFCGC(ID+1)
            IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
            VARLX = VARLX + SRCIW(N)
          ENDDO
!
!---      Global summation  ---
!
          CALL MPI_ALLREDUCE( VARLX,VARIX(NV),NVAR,MPI_REAL8,MPI_SUM,
     &      MPI_COMM_WORLD,IERR )
        ENDIF
      ENDDO
      IOFFSET = OFFSET
!
!---  Loop over number of local reference nodes  ---
!
      DO M = 1,NREFL
        OFFSET = IOFFSET
!
!---    Compute offset  ---
!
        NLINE = NDREFI(M) - 1
        OFFSET = OFFSET + NLINE*(NBYTI + NVREF*NBYTR)
!        PRINT *,'5: OFFSET = ',OFFSET,' ID = ',ID
!
!---    Write global node number  ---
!
        NVAR = 1
        CALL MPI_FILE_WRITE_AT( IWR,OFFSET,ND(NDREFL(M)),NVAR,
     &    MPI_INTEGER,STATUS,IERR)
        OFFSET = OFFSET + NBYTI
!        PRINT *,'6: OFFSET = ',OFFSET,' ID = ',ID
!
!---    Write reference node variable values  ---
!
        N = NDREFL(M)
        DO NV = 1,NVREF
          IRNV = IREF(NV)
!
!---      Aqueous pressure (absolute)  ---
!
          IF( IRNV.EQ.1 ) THEN
            VARX(NV) = PL(2,N) + PATM
!
!---      Gas pressure (absolute)  ---
!
          ELSEIF( IRNV.EQ.2 ) THEN
            VARX(NV) = PG(2,N) + PATM
!
!---      Temperature  ---
!
          ELSEIF( IRNV.EQ.4 ) THEN
            VARX(NV) = T(2,N)
!
!---      Phase condition  ---
!
          ELSEIF( IRNV.EQ.5 ) THEN
            VARX(NV) = REAL( NPHAZ(2,N) )
!
!---      Aqueous pressure (gauge)  ---
!
          ELSEIF( IRNV.EQ.6 ) THEN
            VARX(NV) = PL(2,N)
!
!---      Gas pressure (gauge)  ---
!
          ELSEIF( IRNV.EQ.7 ) THEN
            VARX(NV) = PG(2,N)
!
!---      Apparent aqueous saturation or
!         hydrate specific area (m^2/m^3)  ---
!
          ELSEIF( IRNV.EQ.9 ) THEN
            VARX(NV) = ASL(N)
!
!---      Actual aqueous saturation  ---
!
          ELSEIF( IRNV.EQ.11 ) THEN
            VARX(NV) = SL(2,N)
!
!---      Actual gas saturation  ---
!
          ELSEIF( IRNV.EQ.12 ) THEN
            VARX(NV) = SG(2,N)
!
!---      Aqueous moisture content  ---
!
          ELSEIF( IRNV.EQ.15 ) THEN
            VARX(NV) = SL(2,N)*PORD(2,N)
!
!---      Apparent trapped-gas saturation  ---
!
          ELSEIF( IRNV.EQ.19 ) THEN
            VARX(NV) = ASGT(N)
!
!---      Diffusive porosity  ---
!
          ELSEIF( IRNV.EQ.20 ) THEN
            VARX(NV) = PORD(2,N)
!
!---      Aqueous water mass fraction  ---
!
          ELSEIF( IRNV.EQ.24 ) THEN
            VARX(NV) = XLW(2,N)
!
!---      Aqueous hydraulic head  ---
!
          ELSEIF( IRNV.EQ.27 ) THEN
            VARX(NV) = PL(2,N)/RHORL/GRAV + ZP(N)
!
!---      Rock/soil type  ---
!
          ELSEIF( IRNV.EQ.30 ) THEN
            VARX(NV) = REAL( IZ(N) )
!
!---      Aqueous relative permeability ---
!
          ELSEIF( IRNV.EQ.31 ) THEN
            VARX(NV) = RKL(2,N)
!
!---      Aqueous density ---
!
          ELSEIF( IRNV.EQ.34 ) THEN
            VARX(NV) = RHOL(2,N)
!
!---      Total water mass ---
!
          ELSEIF( IRNV.EQ.37 ) THEN
            VARX(NV) = PORD(2,N)*VOL(N)*XLW(2,N)*SL(2,N)*RHOL(2,N)
!
!---      Water mass source integral ---
!
          ELSEIF( IRNV.EQ.40 ) THEN
            VARX(NV) = SRCIW(N)
!
!---      Matric potential head or Webb matching point head  ---
!
          ELSEIF( IRNV.EQ.63 ) THEN
            VARX(NV) = SCHR(9,N)
!
!---      Stress-xx, Pa  ---
!
          ELSEIF( IRNV.EQ.80 ) THEN
!
!---        Geomechanics  ---
!
            IF( ISLC(50).NE.0 ) THEN
              ISLC50X = ABS(ISLC(50))
              IMODX = MOD(ISLC50X,10)
!
!---          Effective stress  ---
!
              VARX(NV) = SIG_GM(1,N)
!
!---          Total stress  ---
!
              IF( IRNV_GC.EQ.0 ) THEN
!
!---            Poroelasticity ---
!
                IF( IMODX.NE.2 .AND. IMODX.NE.4 ) THEN
                  PX = MAX( PG(2,N),PL(2,N) ) + PATM
                  VARX(NV) = VARX(NV) + PROP_GM(3,N)*(PX-PCMP(N))
                ENDIF
!
!---            Thermoelasticity ---
!
                IF( IMODX.NE.3 .AND. IMODX.NE.4 ) THEN
                  BLK3X = PROP_GM(1,N)/
     &              (1.0D+0-2.0D+0*PROP_GM(2,N))
                  VARX(NV) = VARX(NV) + 
     &              PROP_GM(4,N)*BLK3X*(T(2,N)-TCMP(N))
                ENDIF
              ENDIF
            ELSE
              VARX(NV) = 0.D+0
            ENDIF
!
!---      X node-centered aqueous volumetric flux  ---
!
          ELSEIF( IRNV.EQ.87 ) THEN
            VARX(NV) = 0.5D+0*(UL(1,1,N)+UL(1,2,N))
!
!---      Y node-centered aqueous volumetric flux  ---
!
          ELSEIF( IRNV.EQ.88 ) THEN
            VARX(NV) = 0.5D+0*(VL(1,1,N)+VL(1,2,N))
!
!---      Z node-centered aqueous volumetric flux  ---
!
          ELSEIF( IRNV.EQ.89 ) THEN
            VARX(NV) = 0.5D+0*(WL(1,1,N)+WL(1,2,N))
!
!---      Node number  ---
!
          ELSEIF( IRNV.EQ.100 ) THEN
            VARX(NV) = REAL(N)
!
!---      Stress-yy, Pa  ---
!
          ELSEIF( IRNV.EQ.101 ) THEN
!
!---        Geomechanics  ---
!
            IF( ISLC(50).NE.0 ) THEN
              ISLC50X = ABS(ISLC(50))
              IMODX = MOD(ISLC50X,10)
!
!---          Effective stress  ---
!
              VARX(NV) = SIG_GM(2,N)
!
!---          Total stress  ---
!
              IF( IRNV_GC.EQ.0 ) THEN
!
!---            Poroelasticity ---
!
                IF( IMODX.NE.2 .AND. IMODX.NE.4 ) THEN
                  PX = MAX( PG(2,N),PL(2,N) ) + PATM
                  VARX(NV) = VARX(NV) + PROP_GM(3,N)*(PX-PCMP(N))
                ENDIF
!
!---            Thermoelasticity ---
!
                IF( IMODX.NE.3 .AND. IMODX.NE.4 ) THEN
                  BLK3X = PROP_GM(1,N)/
     &              (1.0D+0-2.0D+0*PROP_GM(2,N))
                  VARX(NV) = VARX(NV) + 
     &              PROP_GM(4,N)*BLK3X*(T(2,N)-TCMP(N))
                ENDIF
              ENDIF
            ENDIF
!
!---      Aqueous compressibility, 1/Pa  ---
!
          ELSEIF( IRNV.EQ.103 ) THEN
            VARX(NV) = TORL(1,N)
!
!---      Trapped-gas saturation ---
!
          ELSEIF( IRNV.EQ.105 ) THEN
            VARX(NV) = SGT(2,N)
!
!---      Minimum effective aqueous saturation  ---
!
          ELSEIF( IRNV.EQ.127 ) THEN
            VARX(NV) = ASLMIN(2,N)
!!
!!---      Stress-zz, Pa  ---
!!
!          ELSEIF( IRNV.EQ.130 ) THEN
!!
!!---        Geomechanics  ---
!!
!            IF( ISLC(50).NE.0 ) THEN
!              ISLC50X = ABS(ISLC(50))
!              IMODX = MOD(ISLC50X,10)
!!
!!---          Effective stress  ---
!!
!              VARX(NV) = SIG_GM(3,N)
!!
!!---          Total stress  ---
!!
!              IF( IRNV_GC.EQ.0 ) THEN
!!
!!---            Poroelasticity ---
!!
!                IF( IMODX.NE.2 .AND. IMODX.NE.4 ) THEN
!                  PX = MAX( PG(2,N),PL(2,N) ) + PATM
!                  VARX(NV) = VARX(NV) + PROP_GM(3,N)*(PX-PCMP(N))
!                ENDIF
!!
!!---            Thermoelasticity ---
!!
!                IF( IMODX.NE.3 .AND. IMODX.NE.4 ) THEN
!                  BLK3X = PROP_GM(1,N)/
!     &              (1.0D+0-2.0D+0*PROP_GM(2,N))
!                  VARX(NV) = VARX(NV) + 
!     &              PROP_GM(4,N)*BLK3X*(T(2,N)-TCMP(N))
!                ENDIF
!              ENDIF
!            ENDIF
!
!---      Water mass source rate (kg/s)  ---
!
          ELSEIF( IRNV.EQ.140 ) THEN
            VARX(NV) = SRCW(2,N)
!
!---      Webb matching point saturation  ---
!
          ELSEIF( IRNV.EQ.150 ) THEN
            VARX(NV) = SCHR(8,N)
!
!---      Aqueous viscosity  ---
!
          ELSEIF( IRNV.EQ.176 ) THEN
            VARX(NV) = VISL(2,N)
!
!---      Integrated water mass  ---
!
          ELSEIF( IRNV.EQ.191 ) THEN
            IUNM = 1
            VARX(NV) = VARIX(NV)
!
!---      Integrated aqueous water mass  ---
!
          ELSEIF( IRNV.EQ.194 ) THEN
            IUNM = 1
            VARX(NV) = VARIX(NV)
!
!---      Stress-yz, Pa  ---
!
          ELSEIF( IRNV.EQ.214 ) THEN
            VARX(NV) = SIG_GM(4,N)
!
!---      Stress-xz, Pa  ---
!
          ELSEIF( IRNV.EQ.224 ) THEN
            VARX(NV) = SIG_GM(5,N)
!
!---      Stress-xy, Pa  ---
!
          ELSEIF( IRNV.EQ.226 ) THEN
            VARX(NV) = SIG_GM(5,N)
!
!---      X-Direction Intrinsic Permeability  ---
!
          ELSEIF( IRNV.EQ.247 ) THEN
            VARX(NV) = PERM(1,N)*PERMRF(2,N)
!
!---      Y-Direction Intrinsic Permeability
!         or Fracture Permeability  ---
!
          ELSEIF( IRNV.EQ.248 ) THEN
            VARX(NV) = PERM(2,N)*PERMRF(2,N)
!
!---      Z-Direction Intrinsic Permeability  ---
!
          ELSEIF( IRNV.EQ.249 ) THEN
            VARX(NV) = PERM(3,N)*PERMRF(2,N)
!
!---      Integrated water mass source  ---
!
          ELSEIF( IRNV.EQ.261 ) THEN
            IUNM = 1
            VARX(NV) = VARIX(NV)
!
!---      X-Node Centroid Position  ---
!
          ELSEIF( IRNV.EQ.291 ) THEN
            VARX(NV) = XP(N)
!
!---      Y-Node Centroid Position  ---
!
          ELSEIF( IRNV.EQ.292 ) THEN
            VARX(NV) = YP(N)
!
!---      Z-Node Centroid Position  ---
!
          ELSEIF( IRNV.EQ.293 ) THEN
            VARX(NV) = ZP(N)
!
!---      Similarity variable  ---
!
          ELSEIF( IRNV.EQ.299 ) THEN
            VARX(NV) = (XP(N)**2)/(TM+SMALL)
!
!---      Strain-xx  ---
!
          ELSEIF( IRNV.EQ.369 ) THEN
            VARX(NV) = EPS_GM(1,N)
!
!---      Strain-yy  ---
!
          ELSEIF( IRNV.EQ.370 ) THEN
            VARX(NV) = EPS_GM(2,N)
!
!---      Strain-zz  ---
!
          ELSEIF( IRNV.EQ.371 ) THEN
            VARX(NV) = EPS_GM(3,N)
!
!---      Reference volumetric strain  ---
!
          ELSEIF( IRNV.EQ.372 ) THEN
            VARX(NV) = EPSV_CMP(N)
!
!---      Volumetric strain  ---
!
          ELSEIF( IRNV.EQ.373 ) THEN
            VARX(NV) = EPSV_GM(2,N)
!
!---      Strain-yz  ---
!
          ELSEIF( IRNV.EQ.372 ) THEN
            VARX(NV) = EPS_GM(4,N)
!
!---      Strain-xz  ---
!
          ELSEIF( IRNV.EQ.373 ) THEN
            VARX(NV) = EPS_GM(5,N)
!
!---      Strain-xy  ---
!
          ELSEIF( IRNV.EQ.374 ) THEN
            VARX(NV) = EPS_GM(6,N)
!!
!!---      X Displacement  ---
!!
!          ELSEIF( IRNV.EQ.375 ) THEN
!            IF( IRNV_GC.EQ.0 ) THEN
!              VARX(NV) = 0.D+0
!              DO L = 1,8
!                NFEN = ND_GM(L,N)
!                VARX(NV) = VARX(NV) + U_GM(2,NFEN) - U_GM(1,NFEN)
!              ENDDO
!              VARX(NV) = 1.25D-1*VARX(NV)
!            ELSE
!              L = -IRNV_GC
!              NFEN = ND_GM(L,N)
!              VARX(NV) = U_GM(2,NFEN) - U_GM(1,NFEN)
!            ENDIF
!!
!!---      Y Displacement  ---
!!
!          ELSEIF( IRNV.EQ.376 ) THEN
!            IF( IRNV_GC.EQ.0 ) THEN
!              VARX(NV) = 0.D+0
!              DO L = 1,8
!                NFEN = ND_GM(L,N)
!                VARX(NV) = VARX(NV) + V_GM(2,NFEN) - V_GM(1,NFEN)
!              ENDDO
!              VARX(NV) = 1.25D-1*VARX(NV)
!            ELSE
!              L = -IRNV_GC
!              NFEN = ND_GM(L,N)
!              VARX(NV) = V_GM(2,NFEN) - V_GM(1,NFEN)
!            ENDIF
!!
!!---      Z Displacement  ---
!!
!          ELSEIF( IRNV.EQ.377 ) THEN
!            IF( IRNV_GC.EQ.0 ) THEN
!              VARX(NV) = 0.D+0
!              DO L = 1,8
!                NFEN = ND_GM(L,N)
!                VARX(NV) = VARX(NV) + W_GM(2,NFEN) - W_GM(1,NFEN)
!              ENDDO
!              VARX(NV) = 1.25D-1*VARX(NV)
!            ELSE
!              L = -IRNV_GC
!              NFEN = ND_GM(L,N)
!              VARX(NV) = W_GM(2,NFEN) - W_GM(1,NFEN)
!            ENDIF
          ELSE
            VARX(NV) = 0.D+0
          ENDIF
!
!---      Solute, conservation-component species, and
!         kinetic-component species reference-node output ---
!
          INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
          IF( IRNV.GT.400 .AND. IRNV.LE.INDX ) THEN
            IF( MOD((IRNV-400),33).EQ.0 ) THEN
              NSL = ((IRNV-400)/33)
              IRNVX = 33
            ELSE
              NSL = ((IRNV-400)/33) + 1
              IRNVX = MOD((IRNV-400),33)
            ENDIF
            IF( NSL.GT.NSOLU ) IUNMOL = 1
!            PRINT *,'IRNVX = ',IRNVX,' NSOLU = ',NSOLU,' ID = ',ID
!
!           401 solute volumetric concentration 'C   '
!           401 total_species volumetric concentration, mol/m^3, 'SP  '
!
            IF( IRNVX.EQ.1 ) THEN
              IUNM = -3
              IF( NSL.GT.NSOLU )THEN
                NEQ = NSL-NSOLU
                VARX(NV) = 0.D+0
                IUNMOL = 1
                DO MX = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(MX+1,NEQ)
                  IF( NSP.LE.NSPL ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NV) = VARX(NV) + EQ_C(MX,NEQ)*SP_CX
                  ENDIF
                ENDDO
              ELSE
                VARX(NV) = C(N,NSL)
              ENDIF
!
!           402 solute aqueous concentration 'CL  '
!           402 total_species aqueous concentration 'SPL  '
!
            ELSEIF( IRNVX.EQ.2 ) THEN
              IUNM = -3
              IF( NSL.GT.NSOLU )THEN
                NEQ = NSL-NSOLU
                VARX(NV) = 0.D+0
                IUNMOL = 1
                DO MX = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(MX+1,NEQ)
                  IF( NSP.LE.NSPL ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NV) = VARX(NV) + EQ_C(MX,NEQ)*SP_CX
                  ENDIF
                ENDDO
              ELSE
                VARX(NV) = C(N,NSL)
              ENDIF
              IF( SL(2,N)*PORD(2,N).GT.SMALL ) THEN
                VARX(NV) = VARX(NV)*YL(N,NSL)/(SL(2,N)*PORD(2,N))
              ELSE
                VARX(NV) = 0.D+0
              ENDIF
!
!           405 solute aqueous mole fraction 'YL  '
!
            ELSEIF( IRNVX.EQ.5 ) THEN
              VARX(NV) = YL(N,NSL)
!
!           408 x-direction node-centered solute/species flux '
!
            ELSEIF( IRNVX.EQ.8 ) THEN
              VARX(NV) = 0.5D+0*(UC(1,N,NSL)+UC(2,N,NSL))
!
!           409 y-direction node-centered solute/species flux '
!
            ELSEIF( IRNVX.EQ.9 ) THEN
              VARX(NV) = 0.5D+0*(VC(1,N,NSL)+VC(2,N,NSL))
!
!           410 z-direction node-centered solute/species flux '
!
            ELSEIF( IRNVX.EQ.10 ) THEN
              VARX(NV) = 0.5D+0*(WC(1,N,NSL)+WC(2,N,NSL))
!
!           411 solute source 'SRC '
!           411 total_species source 'SPSR'
!
            ELSEIF( IRNVX.EQ.11 ) THEN
              VARX(NV) = SRCIC(N,NSL)
!
!           412 exchange total_species volumetric conc., mol/m^3, 'SPX '
!
            ELSEIF( IRNVX.EQ.12 ) THEN
              IUNM = -3
              IF( NSL.GT.NSOLU )THEN
                NEQ = NSL-NSOLU
                VARX(NV) = 0.D+0
                IUNMOL = 1
                DO MX = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(MX+1,NEQ)
                  NSP1 = NSPL + NSPS
                  NSP2 = NSP1 + NSPE
                  IF( NSP.GT.NSP1 .AND. NSP.LE.NSP2 ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NV) = VARX(NV) + EQ_C(MX,NEQ)*SP_CX
                  ENDIF
                ENDDO
              ELSE
                VARX(NV) = YL(N,NSL) + YG(N,NSL)
              ENDIF    
! 
!           413 exchange total_species aqueous conc. mol/m^3, 'SPLX'
!
            ELSEIF( IRNVX.EQ.13 ) THEN
              IUNM = -3
              IF( NSL.GT.NSOLU )THEN
                NEQ = NSL-NSOLU
                VARX(NV) = 0.D+0
                IUNMOL = 1
                DO MX = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(MX+1,NEQ)
                  NSP1 = NSPL + NSPS
                  NSP2 = NSP1 + NSPE
                  IF( NSP.GT.NSP1 .AND. NSP.LE.NSP2 ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NV) = VARX(NV) + EQ_C(MX,NEQ)*SP_CX
                  ENDIF
                ENDDO
                IF( SL(2,N)*PORD(2,N).GT.SMALL ) THEN
                  VARX(NV) = VARX(NV)*YL(N,NSL)/(SL(2,N)*PORD(2,N))
                ELSE
                  VARX(NV) = 0.D+0
                ENDIF
              ELSE
                VARX(NV) = C(N,NSL)
              ENDIF   
!
!           426 exchange total_species volumetric conc., mol/m^3, 'SPS '
!
!           ELSEIF( IRNVX.EQ.26 ) THEN
              IUNM = -3
              IF( NSL.GT.NSOLU )THEN
                NEQ = NSL-NSOLU
                VARX(NV) = 0.D+0
                IUNMOL = 1
                DO MX = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(MX+1,NEQ)
                  NSP1 = NSPL + NSPS
                  IF( NSP.GT.NSPL .AND. NSP.LE.NSP1 ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NV) = VARX(NV) + EQ_C(MX,NEQ)*SP_CX
                  ENDIF
                ENDDO
              ENDIF     
!
!           427 exchange total_species aqueous conc. mol/m^3, 'SPLS'
!
            ELSEIF( IRNVX.EQ.27 ) THEN
              IUNM = -3
              IF( NSL.GT.NSOLU )THEN
                NEQ = NSL-NSOLU
                VARX(NV) = 0.D+0
                IUNMOL = 1
                DO MX = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(MX+1,NEQ)
                  NSP1 = NSPL + NSPS
                  IF( NSP.GT.NSPL .AND. NSP.LE.NSP1 ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NV) = VARX(NV) + EQ_C(MX,NEQ)*SP_CX
                  ENDIF
                ENDDO
                IF( SL(2,N)*PORD(2,N).GT.SMALL ) THEN
                  VARX(NV) = VARX(NV)*YL(N,NSL)/(SL(2,N)*PORD(2,N))
                ELSE
                  VARX(NV) = 0.D+0
                ENDIF
              ENDIF     
            ENDIF    
          ENDIF
!
!---      Unit conversion  ---
!
          VARX(NV) = VARX(NV)*CNVREF(NV)
        ENDDO
!        DO NV = 1,NVREF
!          PRINT *,'VARX(',NV,') = ',VARX(NV),' ID = ',ID
!        ENDDO
!
!---    Write time step, time and variables to output.bin file  ---
!
        NVAR = NVREF
        CALL MPI_FILE_WRITE_AT( IWR,OFFSET,VARX,NVAR,MPI_REAL8,
     &    STATUS,IERR)
      ENDDO
!
!---  Increment integrated offset by 2 x integer (time step number and 
!     number of iterations), 2 x real (time and time step) + 
!     number of reference nodes x ( 1 x integer (global node number)
!     + number of reference node variables x real ) ---
!
      IOFFSET_REF = IOFFSET_REF + 2*NBYTI + 2*NBYTR + 
     &  NREF*(NBYTI + NVREF*NBYTR)
!      PRINT *,'IOFFSET_REF = ',IOFFSET_REF,' ID = ',ID
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REFNOD_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFIN_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Write a surface.bin file
!
!     Surface flux and source term integrator
!
!   2   aqueous volumetric flux - ULV, VLV, WLV
!   5   aqueous mass flux - ULM, VLM, WLM
! 101 to 100+NSOLU
!       solute flux - UC, VC, WC
! 101+NSOLU to 100+NSOLU+NEQC   
!       conservation-component species flux - UCC, VCC, WCC
! 101+NSOLU+NEQC to 100+NSOLU+NEQC+NEQK
!       kinetic-component species flux - UKC, VKC, WKC
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 February 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE REACT
      USE PROP
      USE OUTPU
      USE HYST
      USE GRID
      USE GEO_MECH
      USE FLUX
      USE FILES
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
      REAL, DIMENSION(:,:), ALLOCATABLE :: SF_L,SF_G
      REAL, DIMENSION(2) :: VARZ
      INTEGER NCTX(LSF),ISFCX(6)
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      LOGICAL IFLAG
      CHARACTER(32) :: CHMSG
      CHARACTER(3) :: CHX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFIN_W'
!
!---  Integrate source terms, looping over all nodes, 
!     skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        SRCIW(N) = SRCIW(N) + SRCW(2,N)*DT
        SRCIA(N) = SRCIA(N) + SRCA(2,N)*DT
        SRCIS(N) = SRCIS(N) + SRCS(2,N)*DT
      ENDDO
!
!---  Return if there are no surface flux outputs  ---
!
      IF( NSF.EQ.0 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Allocate temporary memory for surface flux data  ---
!
      ALLOCATE( SF_L(1:2,1:NSF),STAT=ISTAT )
      CHMSG = 'SF_L'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SF_G(1:2,1:NSF),STAT=ISTAT )
      CHMSG = 'SF_G'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      DO NS = 1,NSF
        DO M = 1,2
          SF_L(M,NS) = 0.D+0
          SF_G(M,NS) = 0.D+0
        ENDDO
      ENDDO
!
!---  Create an surface.bin file and write header information  ---
!
      IF( IHSF.EQ.0 ) THEN
        IHSF = 1
!
!---    Set the boundary condition pointer array,
!       looping over the local number of surface-flux nodes  ---
!
        DO NC = 1,NSFN(ID+1)
          N = ISFN(NC)
          NS = ISFS(NC)
          ISFB(NC) = 0
!
!---      Bottom surface  ---
!
          IF( ISFD(NS).EQ.-3 ) THEN
            NBX = ICM(1,N)
            IF( NBX.EQ.0 ) THEN
              DO NB = 1,NBC(ID+1)
                IF( IBCN(NB).EQ.N .AND. IBCD(NB).EQ.ISFD(NS) ) THEN
                  ISFB(NC) = NB
                  EXIT
                ENDIF
              ENDDO
            ENDIF
!
!---      South surface  ---
!
          ELSEIF( ISFD(NS).EQ.-2 ) THEN
            NSX = ICM(2,N)
            IF( NSX.EQ.0 ) THEN
              DO NB = 1,NBC(ID+1)
                IF( IBCN(NB).EQ.N .AND. IBCD(NB).EQ.ISFD(NS) ) THEN
                  ISFB(NC) = NB
                  EXIT
                ENDIF
              ENDDO
            ENDIF
!
!---      West surface  ---
!
          ELSEIF( ISFD(NS).EQ.-1 ) THEN
            NWX = ICM(3,N)
            IF( NWX.EQ.0 ) THEN
              DO NB = 1,NBC(ID+1)
                IF( IBCN(NB).EQ.N .AND. IBCD(NB).EQ.ISFD(NS) ) THEN
                  ISFB(NC) = NB
                  EXIT
                ENDIF
              ENDDO
            ENDIF
!
!---      East surface  ---
!
          ELSEIF( ISFD(NS).EQ.1 ) THEN
            NEX = ICM(4,N)
            IF( NEX.EQ.0 ) THEN
              DO NB = 1,NBC(ID+1)
                IF( IBCN(NB).EQ.N .AND. IBCD(NB).EQ.ISFD(NS) ) THEN
                  ISFB(NC) = NB
                  EXIT
                ENDIF
              ENDDO
            ENDIF
!
!---      North surface  ---
!
          ELSEIF( ISFD(NS).EQ.2 ) THEN
            NNX = ICM(5,N)
            IF( NNX.EQ.0 ) THEN
              DO NB = 1,NBC(ID+1)
                IF( IBCN(NB).EQ.N .AND. IBCD(NB).EQ.ISFD(NS) ) THEN
                  ISFB(NC) = NB
                  EXIT
                ENDIF
              ENDDO
            ENDIF
!
!---      Top surface  ---
!
          ELSEIF( ISFD(NS).EQ.3 ) THEN
            NTX = ICM(6,N)
            IF( NTX.EQ.0 ) THEN
              DO NB = 1,NBC(ID+1)
                IF( IBCN(NB).EQ.N .AND. IBCD(NB).EQ.ISFD(NS) ) THEN
                  ISFB(NC) = NB
                  EXIT
                ENDIF
              ENDDO
            ENDIF
          ENDIF        
        ENDDO 
!
!---    Skip the default surface file if its unused  ---
!
        NSTART = 1
        IF( ISFGP(1).EQ.0 ) NSTART = 2
!
!---    Loop over the number of surface-flux files  ---
!
        DO NSG = NSTART,NSFGP
!
!---      Set the offset point to the integrated offset for all
!         processors  ---
!
          OFFSET = IOFFSET_SF(NSG)
!
!---      Open the surface file  ---
!
          NCH = INDEX(FNSF(NSG)(1:),'  ')-1
          CALL MPI_FILE_OPEN( MPI_COMM_WORLD,FNSF(NSG)(1:NCH),
     &      MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,
     &      ISF(NSG),IERR )
!
!---      Write operational mode index  ---
!
          NVAR = 1
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,
     &       IOM,NVAR,MPI_INTEGER,STATUS,IERR)
          OFFSET = OFFSET + NVAR*NBYTI
          IOFFSET_SF(NSG) = IOFFSET_SF(NSG) + NVAR*NBYTI
!
!---      Write number of surfaces  ---
!
          NVAR = 1
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,
     &       ISFGP(NSG),NVAR,MPI_INTEGER,STATUS,IERR)
          OFFSET = OFFSET + NVAR*NBYTI
          IOFFSET_SF(NSG) = IOFFSET_SF(NSG) + NVAR*NBYTI
!
!---      Write time unit string  ---
!
          NVAR = 64
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,UNTM,
     &       NVAR,MPI_CHAR,STATUS,IERR)
          OFFSET = OFFSET + NVAR
          IOFFSET_SF(NSG) = IOFFSET_SF(NSG) + NVAR
!
!---      Loop over surfaces  ---
!
          DO NS = 1,NSF
!
!---        Surface not associated with surface file,
!           skip output  ---
!
            IF( NSG.NE.ISFF(NS) ) CYCLE
!
!---        Write surface type  ---
!
            NVAR = 1
            IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,
     &        ISFT(NS),NVAR,MPI_INTEGER,STATUS,IERR)
            OFFSET = OFFSET + NVAR*NBYTI
            IOFFSET_SF(NSG) = IOFFSET_SF(NSG) + NVAR*NBYTI
!
!---        Write surface flux header character strings  ---
!
            NCH1 = INDEX(CHSF(1,NS),'  ')-1
            NCH2 = INDEX(CHSF(2,NS),'  ')-1
            NVAR = 64
            IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,
     &        CHSF(1,NS),NVAR,MPI_CHAR,STATUS,IERR)
            OFFSET = OFFSET + NVAR
            IOFFSET_SF(NSG) = IOFFSET_SF(NSG) + NVAR
            NVAR = 64
            IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,
     &        CHSF(2,NS),NVAR,MPI_CHAR,STATUS,IERR)
            OFFSET = OFFSET + NVAR
            IOFFSET_SF(NSG) = IOFFSET_SF(NSG) + NVAR
!
!---        Write header rate and integral units  ---
!
            NCH1 = INDEX(UNSF(1,NS),'  ')-1
            NCH2 = INDEX(UNSF(2,NS),'  ')-1
            NVAR = 64
            IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,
     &        UNSF(1,NS),NVAR,MPI_CHAR,STATUS,IERR)
            OFFSET = OFFSET + NVAR
            IOFFSET_SF(NSG) = IOFFSET_SF(NSG) + NVAR
            NVAR = 64
            IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,
     &        UNSF(2,NS),NVAR,MPI_CHAR,STATUS,IERR)
            OFFSET = OFFSET + NVAR
            IOFFSET_SF(NSG) = IOFFSET_SF(NSG) + NVAR
          ENDDO
        ENDDO
      ENDIF
!
!---  Loop over the local number of surface-flux nodes  ---
!
      DO NC = 1,NSFN(ID+1)
        N = ISFN(NC)
        NS = ISFS(NC)
!
!---    Bottom surface  ---
!
        IF( ISFD(NS).EQ.-3 ) THEN
!
!---      Aqueous volumetric flux, m^3/s  ---
!
          IF( ISFT(NS).EQ.2 ) THEN
            SFX = WL(1,1,N)*AFZ(1,N)
!
!---      Aqueous mass flux, kg/s  ---
!
          ELSEIF( ISFT(NS).EQ.5 ) THEN
            RHOLX = RHOL(2,N)
            IF( WL(1,1,N).GT.ZERO ) THEN
              NBX = ICM(1,N)
              IF( NBX.EQ.0 ) THEN
                NB = ISFB(NC)
                RHOLX = RHOLB(2,NB)
              ELSE
                RHOLX = RHOL(2,NBX)
              ENDIF
            ENDIF
            SFX = WL(1,1,N)*AFZ(1,N)*RHOLX
          ELSEIF( ISFT(NS).GT.100 .AND. ISFT(NS).LE.(100+NSOLU) ) THEN
            NSL = ISFT(NS)-100
            SFX = WC(1,N,NSL)*AFZ(1,N)
          ELSEIF( ISFT(NS).GT.(100+NSOLU) .AND. 
     &      ISFT(NS).LE.(100+NSOLU+NEQC+NEQK) ) THEN
            NSL = ISFT(NS)-100
            SFX = WC(1,N,NSL)*AFZ(1,N)*1.D-3
          ENDIF
!
!---    South surface  ---
!
        ELSEIF( ISFD(NS).EQ.-2 ) THEN
!
!---      Aqueous volumetric flux, m^3/s  ---
!
          IF( ISFT(NS).EQ.2 ) THEN
            SFX = VL(1,1,N)*AFY(1,N)
!
!---      Aqueous mass flux, kg/s  ---
!
          ELSEIF( ISFT(NS).EQ.5 ) THEN
            RHOLX = RHOL(2,N)
            IF( VL(1,1,N).GT.ZERO ) THEN
              NSX = ICM(2,N)
              IF( NSX.EQ.0 ) THEN
                NB = ISFB(NC)
                RHOLX = RHOLB(2,NB)
              ELSE
                RHOLX = RHOL(2,NSX)
              ENDIF
            ENDIF
            SFX = VL(1,1,N)*AFY(1,N)*RHOLX
          ELSEIF( ISFT(NS).GT.100 .AND. ISFT(NS).LE.(100+NSOLU) ) THEN
            NSL = ISFT(NS)-100
            SFX = VC(1,N,NSL)*AFY(1,N)
          ELSEIF( ISFT(NS).GT.(100+NSOLU) .AND. 
     &      ISFT(NS).LE.(100+NSOLU+NEQC+NEQK) ) THEN
            NSL = ISFT(NS)-100
            SFX = VC(1,N,NSL)*AFY(1,N)*1.D-3
          ENDIF
!
!---    West surface  ---
!
        ELSEIF( ISFD(NS).EQ.-1 ) THEN
!
!---      Aqueous volumetric flux, m^3/s  ---
!
          IF( ISFT(NS).EQ.2 ) THEN
            SFX = UL(1,1,N)*AFX(1,N)
!
!---      Aqueous mass flux, kg/s  ---
!
          ELSEIF( ISFT(NS).EQ.5 ) THEN
            RHOLX = RHOL(2,N)
            IF( UL(1,1,N).GT.ZERO ) THEN
              NWX = ICM(3,N)
              IF( NWX.EQ.0 ) THEN
                NB = ISFB(NC)
                RHOLX = RHOLB(2,NB)
              ELSE
                RHOLX = RHOL(2,NWX)
              ENDIF
            ENDIF
            SFX = UL(1,1,N)*AFX(1,N)*RHOLX
          ELSEIF( ISFT(NS).GT.100 .AND. ISFT(NS).LE.(100+NSOLU) ) THEN
            NSL = ISFT(NS)-100
            SFX = UC(1,N,NSL)*AFX(1,N)
          ELSEIF( ISFT(NS).GT.(100+NSOLU) .AND. 
     &      ISFT(NS).LE.(100+NSOLU+NEQC+NEQK) ) THEN
            NSL = ISFT(NS)-100
            SFX = UC(1,N,NSL)*AFX(1,N)*1.D-3
          ENDIF
!
!---    East surface  ---
!
        ELSEIF( ISFD(NS).EQ.1 ) THEN
!
!---      Aqueous volumetric flux, m^3/s  ---
!
          IF( ISFT(NS).EQ.2 ) THEN
            SFX = UL(1,2,N)*AFX(2,N)
!
!---      Aqueous mass flux, kg/s  ---
!
          ELSEIF( ISFT(NS).EQ.5 ) THEN
            RHOLX = RHOL(2,N)
            IF( UL(1,2,N).LT.ZERO ) THEN
              NEX = ICM(4,N)
              IF( NEX.EQ.0 ) THEN
                NB = ISFB(NC)
                RHOLX = RHOLB(2,NB)
              ELSE
                RHOLX = RHOL(2,NEX)
              ENDIF
            ENDIF
            SFX = UL(1,2,N)*AFX(2,N)*RHOLX
          ELSEIF( ISFT(NS).GT.100 .AND. ISFT(NS).LE.(100+NSOLU) ) THEN
            NSL = ISFT(NS)-100
            SFX = UC(2,N,NSL)*AFX(2,N)
          ELSEIF( ISFT(NS).GT.(100+NSOLU) .AND. 
     &      ISFT(NS).LE.(100+NSOLU+NEQC+NEQK) ) THEN
            NSL = ISFT(NS)-100
            SFX = UC(2,N,NSL)*AFX(2,N)*1.D-3
          ENDIF
!
!---    North surface  ---
!
        ELSEIF( ISFD(NS).EQ.2 ) THEN
!
!---      Aqueous volumetric flux, m^3/s  ---
!
          IF( ISFT(NS).EQ.2 ) THEN
            SFX = VL(1,2,N)*AFY(2,N)
!
!---      Aqueous mass flux, kg/s  ---
!
          ELSEIF( ISFT(NS).EQ.5 ) THEN
            RHOLX = RHOL(2,N)
            IF( VL(1,2,N).LT.ZERO ) THEN
              NNX = ICM(4,N)
              IF( NNX.EQ.0 ) THEN
                NB = ISFB(NC)
                RHOLX = RHOLB(2,NB)
              ELSE
                RHOLX = RHOL(2,NNX)
              ENDIF
            ENDIF
            SFX = VL(1,2,N)*AFY(2,N)*RHOLX
          ELSEIF( ISFT(NS).GT.100 .AND. ISFT(NS).LE.(100+NSOLU) ) THEN
            NSL = ISFT(NS)-100
            SFX = VC(2,N,NSL)*AFY(2,N)
          ELSEIF( ISFT(NS).GT.(100+NSOLU) .AND. 
     &      ISFT(NS).LE.(100+NSOLU+NEQC+NEQK) ) THEN
            NSL = ISFT(NS)-100
            SFX = VC(2,N,NSL)*AFY(2,N)*1.D-3
          ENDIF
!
!---    Top surface  ---
!
        ELSEIF( ISFD(NS).EQ.3 ) THEN
!
!---      Aqueous volumetric flux, m^3/s  ---
!
          IF( ISFT(NS).EQ.2 ) THEN
            SFX = WL(1,2,N)*AFZ(2,N)
!
!---      Aqueous mass flux, kg/s  ---
!
          ELSEIF( ISFT(NS).EQ.5 ) THEN
            RHOLX = RHOL(2,N)
            IF( WL(1,2,N).LT.ZERO ) THEN
              NTX = ICM(6,N)
              IF( NTX.EQ.0 ) THEN
                NB = ISFB(NC)
                RHOLX = RHOLB(2,NB)
              ELSE
                RHOLX = RHOL(2,NTX)
              ENDIF
            ENDIF
            SFX = WL(1,2,N)*AFZ(2,N)*RHOLX
          ELSEIF( ISFT(NS).GT.100 .AND. ISFT(NS).LE.(100+NSOLU) ) THEN
            NSL = ISFT(NS)-100
            SFX = WC(2,N,NSL)*AFZ(2,N)
          ELSEIF( ISFT(NS).GT.(100+NSOLU) .AND. 
     &      ISFT(NS).LE.(100+NSOLU+NEQC+NEQK) ) THEN
            NSL = ISFT(NS)-100
            SFX = WC(2,N,NSL)*AFZ(2,N)*1.D-3
          ENDIF
        ENDIF        
        RSFDX = REAL(ISFD(NS))
        SF_L(1,NS) = SF_L(1,NS) + SFX*SIGN(1.D+0,RSFDX)
        SF_L(2,NS) = SF_L(2,NS) + SFX*SIGN(1.D+0,RSFDX)*DT
      ENDDO 
!
!---  Global values of surface fluxes  ---
!
      NVAR = 1
      DO NS = 1,NSF
        DO M = 1,2
          VARLX = SF_L(M,NS)
          CALL MPI_ALLREDUCE( VARLX,VARGX,NVAR,MPI_REAL8,
     &      MPI_SUM,MPI_COMM_WORLD,IERR )
          SF_G(M,NS) = VARGX
        ENDDO
      ENDDO
!
!---  Skip the default surface file if its unused  ---
!
      NSTART = 1
      IF( ISFGP(1).EQ.0 ) NSTART = 2
!
!---  Loop over the number of surface-flux files  ---
!
      DO NSG = NSTART,NSFGP
!
!---    Set the offset point to the integrated offset for all
!       processors  ---
!
        OFFSET = IOFFSET_SF(NSG)
!
!---    Write current time to binary surface files  --
!
        NVAR = 1
        VARX = TM*CNVTM
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,VARX,
     &     NVAR,MPI_REAL8,STATUS,IERR )
        OFFSET = OFFSET + NVAR*NBYTR
!
!---    Loop over surfaces  ---
!
        DO NS = 1,NSF
!
!---      Surface not associated with surface file  ---
!
          IF( NSG.NE.ISFF(NS) ) EXIT
          NVAR = 1
          VARX = SF_G(1,NS)*CNVSF(1,NS)
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,VARX,
     &       NVAR,MPI_REAL8,STATUS,IERR )
          OFFSET = OFFSET + NVAR*NBYTR
          NVAR = 1
          VARX = SF_G(2,NS)*CNVSF(2,NS)
          IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( ISF(NSG),OFFSET,VARX,
     &       NVAR,MPI_REAL8,STATUS,IERR )
          OFFSET = OFFSET + NVAR*NBYTR
        ENDDO
!
!---    Reset the integrated offset  ---
!
        IOFFSET_SF(NSG) = OFFSET
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFIN_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRPLOT_W
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Write a plot_xxx.bin file.
!
!     1 aqueous pressure (Pa)
!     2 gas pressure (Pa)
!     3 surface vapor pressure (Pa)
!     4 temperature (C)
!     5 phase condition (null) 'PHCN'
!     6 aqueous gauge pressure (Pa)
!     7 gas gauge pressure (Pa)
!     9 apparent aqueous saturation (null), 'ASL '
!     11 aqueous saturation (null), ' SL '
!     12 gas saturation (null), ' SG '
!     15 aqueous moisture content (null), 'MCL '
!     19 effective trapped gas (null), 'ESGT'
!     20 diffusive porosity (null), 'PORD'
!     24 aqueous water mass fraction (null), 'XLW '
!     27 aqueous hydraulic head (m), 'HHL '
!     28 gas hydraulic head (m), 'HHG '
!     30 rock/soil type (null), 'RSZN'
!     31 aqueous relative permeability (null), 'RKL '
!     34 aqueous density (kg/m^3), 'RHOL'
!     37 total water mass (kg), 'TMW '
!     40 water mass source integral (kg), 'SRIW'
!     49 aqueous courant (null), 'CRNL'
!     63 Webb matching point head (m), 'WMPH'
!     80 stress-xx (Pa), 'SIGXX'
!     83 aqueous matrix (null), 'DSLM'
!     84 aqueous fracture (null), 'DSLF'
!     87 xnc aqueous volumetric flux (m/s), 'ULNC'
!     88 ync aqueous volumetric flux (m/s), 'VLNC'
!     89 znc aqueous volumetric flux (m/s), 'WLNC'
!     100 node number (null), 'NODE'
!     101 stress-yy (Pa), 'SIGXX'
!     103 aqueous compressibility (1/Pa), CMPL
!     105 trapped gas saturation (null), 'SGT '
!     127 minimum effect aqueous saturation (null), 'ESLM'
!     130 stress-zz (Pa), 'SIGZZ'
!     132 X displacement variance (m^2), 'SIGXX'
!     133 Y displacement variance (m^2), 'SIGYY'
!     134 Z displacement variance (m^2), 'SIGZZ'
!     137 water relative humidity (null), 'RHW '
!     140 water mass source rate (kg/s), 'SRCW'
!     150 Webb matching point saturation (null), 'WMPS'
!     176 aqueous viscosity 'VISL'
!     200 reserved to control plot file output
!     201 x aqueous relative permeability 'RKLX'
!     202 y aqueous relative permeability 'RKLY'
!     203 z aqueous relative permeability 'RKLZ'
!     214 stress-yz (Pa), 'SIGYZ'
!     224 stress-xz (Pa), 'SIGXZ'
!     226 stress-xy (Pa), 'SIGXY'
!     247 x-direction intrinsic permeability (m^2), ' UK '
!     248 y-direction intrinsic permeability (m^2), ' VK '
!     249 z-direction intrinsic permeability (m^2), ' WK '
!     291 x node position, 'XP'
!     292 y node position, 'YP'
!     293 z node position, 'ZP'
!     299 similitude parameter (m^2/s), 'SIMV'
!     369 strain xx, EPSXX
!     370 strain yy, EPSYY
!     371 strain zz, EPSZZ
!     372 strain yz, EPSYZ
!     373 strain xz, EPSXZ
!     374 strain xy, EPSXY
!     401 solute volumetric concentration 'C   '
!     401 total_species volumetric concentration, mol/m^3, 'SP  '
!     402 solute aqueous concentration 'CL  '
!     402 total_species aqueous concentration 'SPL  '
!     405 solute aqueous mole fraction 'YL  '
!     408 xnc solute flux 'UC  '
!     409 ync solute flux 'VC  '
!     410 znc solute flux 'WC  '
!     411 solute source 'SRC '
!     411 total_species source 'SPSR'
!     412 exchange total_species volumetric conc., mol/m^3, 'SPX '
!     413 exchange total_species aqueous concentration mol/m^3, 'SPLX'
!     423 solute integrated mass 'ICM'
!     423 total_species integrated mass 'SPIM'
!     426 exchange total_species volumetric conc., mol/m^3, 'SPS '
!     427 exchange total_species aqueous concentration mol/m^3, 'SPLS'
!     430 solute mass concentration, 1/kg soil, 'CM'
!     431 solute activity, Bq, 'CA'
!     432 integrated solute activity, Bq, 'CIA'
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 November 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE REACT
      USE PROP
      USE OUTPU
      USE HYST
      USE GRID
      USE GEO_MECH
      USE FLUX
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
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARX,VVARX
      CHARACTER*16 FN
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: PLTVNMX
      CHARACTER(64) :: SPNMX
      CHARACTER(32) :: CHMSG
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      INTEGER MXZ(4),MXY(4),MYZ(4),MX(2),MY(2),MZ(2)
!
!----------------------Data Statements---------------------------------!
!
      DATA MXZ / 1,2,5,6 /
      DATA MXY / 1,2,3,4 /
      DATA MYZ / 1,3,5,7 /
      DATA MX / 1,2 /
      DATA MY / 1,3 /
      DATA MZ / 1,5 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRPLOT_W'
!
!---  Write vertex files once  ---
!
      IF( IPLSKP.EQ.0 ) THEN
!
!---    Open plot_xgrid.bin file  ---
!
        CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'plot_xgrid.bin',
     &    MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,IPLX,IERR )
!
!---    Open plot_ygrid.bin file  ---
!
        CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'plot_ygrid.bin',
     &    MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,IPLY,IERR )
!
!---    Open plot_zgrid.bin file  ---
!
        CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'plot_zgrid.bin',
     &    MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,IPLZ,IERR )
      ENDIF
!
!---  Create a plot file name  ---
!
      FN = 'plot_0000000.bin'
      IF( NSTEP.LE.9 ) THEN
        FN(12:12) = CHAR(48+NSTEP)
      ELSEIF( NSTEP.LE.99 ) THEN
        I1 = INT(NSTEP/10)
        I2 = MOD( NSTEP,10 )
        FN(11:11) = CHAR(48+I1)
        FN(12:12) = CHAR(48+I2)
      ELSEIF( NSTEP.LE.999 ) THEN
        I1 = INT(NSTEP/100)
        I2 = INT( MOD( NSTEP,100 )/10 )
        I3 = MOD( NSTEP,10 )
        FN(10:10) = CHAR(48+I1)
        FN(11:11) = CHAR(48+I2)
        FN(12:12) = CHAR(48+I3)
      ELSEIF( NSTEP.LE.9999 ) THEN
        I1 = INT(NSTEP/1000)
        I2 = INT( MOD( NSTEP,1000 )/100 )
        I3 = INT( MOD( NSTEP,100 )/10 )
        I4 = MOD( NSTEP,10 )
        FN(9:9) = CHAR(48+I1)
        FN(10:10) = CHAR(48+I2)
        FN(11:11) = CHAR(48+I3)
        FN(12:12) = CHAR(48+I4)
      ELSEIF( NSTEP.LE.99999 ) THEN
        I1 = INT(NSTEP/10000)
        I2 = INT( MOD( NSTEP,10000 )/1000 )
        I3 = INT( MOD( NSTEP,1000 )/100 )
        I4 = INT( MOD( NSTEP,100 )/10 )
        I5 = MOD( NSTEP,10 )
        FN(8:8) = CHAR(48+I1)
        FN(9:9) = CHAR(48+I2)
        FN(10:10) = CHAR(48+I3)
        FN(11:11) = CHAR(48+I4)
        FN(12:12) = CHAR(48+I5)
      ELSEIF( NSTEP.LE.999999 ) THEN
        I1 = INT(NSTEP/100000)
        I2 = INT( MOD( NSTEP,100000 )/10000 )
        I3 = INT( MOD( NSTEP,10000 )/1000 )
        I4 = INT( MOD( NSTEP,1000 )/100 )
        I5 = INT( MOD( NSTEP,100 )/10 )
        I6 = MOD( NSTEP,10 )
        FN(7:7) = CHAR(48+I1)
        FN(8:8) = CHAR(48+I2)
        FN(9:9) = CHAR(48+I3)
        FN(10:10) = CHAR(48+I4)
        FN(11:11) = CHAR(48+I5)
        FN(12:12) = CHAR(48+I6)
      ELSEIF( NSTEP.LE.9999999 ) THEN
        I1 = INT(NSTEP/1000000)
        I2 = INT( MOD( NSTEP,1000000 )/100000 )
        I3 = INT( MOD( NSTEP,100000 )/10000 )
        I4 = INT( MOD( NSTEP,10000 )/1000 )
        I5 = INT( MOD( NSTEP,1000 )/100 )
        I6 = INT( MOD( NSTEP,100 )/10 )
        I7 = MOD( NSTEP,10 )
        FN(6:6) = CHAR(48+I1)
        FN(7:7) = CHAR(48+I2)
        FN(8:8) = CHAR(48+I3)
        FN(9:9) = CHAR(48+I4)
        FN(10:10) = CHAR(48+I5)
        FN(11:11) = CHAR(48+I6)
        FN(12:12) = CHAR(48+I7)
      ELSE
        FN = 'plot_xxxxxxx.bin'
      ENDIF
!
!---  Open plot_xxxxxxx.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,FN,
     &  MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,IPL,IERR )
!
!---  Allocate temporary memory for dummy integer array  ---
!
      ALLOCATE( IVARX(1:14),STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Write time step, global number of nodes, number of 
!     plot variables and time (sec)  ---
!
      NVPLOTX = 0
      DO NV = 1,NVPLOT
        IPNV = IPLOT(NV)
        IF( IPNV.EQ.200 ) CYCLE
        NVPLOTX = NVPLOTX + 1
      ENDDO
      IOFFSET = 0
      NVAR = 14
      IVARX(1) = NSTEP
      IVARX(2) = NFLD
      IVARX(3) = IFLD
      IVARX(4) = JFLD
      IVARX(5) = KFLD
      IVARX(6) = NPFLD
      IVARX(7) = IPFLD
      IVARX(8) = JPFLD
      IVARX(9) = KPFLD
      IVARX(10) = NXP
      IVARX(11) = NVPLOTX
      IVARX(12) = ISGNP
      IVARX(13) = ICS
      IVARX(14) = ISLC(63)
      OFFSET = IOFFSET
      IOFFSET = IOFFSET + NVAR*NBYTI
!      PRINT *,'1: OFFSET = ',OFFSET,' IOFFSET = ',IOFFSET,' ID = ',ID
      IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IPL,OFFSET,IVARX,NVAR,
     &   MPI_INTEGER,STATUS,IERR )
      NVAR = 1
      OFFSET = IOFFSET
      IOFFSET = IOFFSET + NBYTR
!      PRINT *,'2: OFFSET = ',OFFSET,' IOFFSET = ',IOFFSET,' ID = ',ID
      IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IPL,OFFSET,TM,NVAR,
     &   MPI_REAL8,STATUS,IERR )
!
!---  Deallocate temporary memory for dummy integer array  ---
!
      DEALLOCATE( IVARX,STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Write time units, length units, and plot variable units  ---
!
      NVAR = 64
      OFFSET = IOFFSET
      IOFFSET = IOFFSET + NVAR
!      PRINT *,'3: OFFSET = ',OFFSET,' IOFFSET = ',IOFFSET,' ID = ',ID
      IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IPL,OFFSET,UNTM,NVAR,
     &   MPI_CHAR,STATUS,IERR )
      NVAR = 64
      OFFSET = IOFFSET
      IOFFSET = IOFFSET + NVAR
 !     PRINT *,'4: OFFSET = ',OFFSET,' IOFFSET = ',IOFFSET,' ID = ',ID
      IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IPL,OFFSET,UNLN,NVAR,
     &   MPI_CHAR,STATUS,IERR )
!
!---  Plane node counts  ---
!
      IJFLD = IFLD*JFLD
      JKFLD = JFLD*KFLD
      KIFLD = KFLD*IFLD
!
!---  Number of vertices: cylindrical coordinates  ---
!
      IF( ICS.EQ.2 ) THEN
        IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1) .OR.
     &    IJFLD.EQ.1 .OR. JKFLD.EQ.1 .OR. KIFLD.EQ.1 .OR. 
     &    ISLC(63).EQ.1 ) THEN
          NVTC = 8
        ELSEIF( IFLD.GT.1 .AND. JFLD.GT.1 ) THEN
          NVTC = 4
        ELSEIF( IFLD.GT.1 .AND. KFLD.GT.1 ) THEN
          NVTC = 4
        ELSEIF( IFLD.GT.1 ) THEN
          NVTC = 2
        ENDIF
!
!---  Number of vertices: non-cylindrical coordinates  ---
!
      ELSE
        IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1) .OR.
     &    ICS.EQ.3 .OR. ISLC(63).EQ.1 .OR.
     &    IJFLD.EQ.1 .OR. JKFLD.EQ.1 .OR. KIFLD.EQ.1 ) THEN
          NVTC = 8
        ELSEIF( IFLD.GT.1 .AND. JFLD.GT.1 ) THEN
          NVTC = 4
        ELSEIF( IFLD.GT.1 .AND. KFLD.GT.1 ) THEN
          NVTC = 4
        ELSEIF( IFLD.GT.1 ) THEN
          NVTC = 2
        ENDIF
      ENDIF
!
!---  Write number of vertices to binary plot file  ---
!
      NVAR = 1
      OFFSET = IOFFSET
      IOFFSET = IOFFSET + NVAR*NBYTI
      IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IPL,OFFSET,NVTC,NVAR,
     &  MPI_INTEGER,STATUS,IERR )
!
!---  Write vertex files once  ---
!
      IF( IPLSKP.EQ.0 ) THEN
!
!---    Allocate temporary memory for vertices  ---
!
        ALLOCATE( VARX(1:NVTC*NFC(ID+1)),STAT=ISTAT )
        CHMSG = 'VARX'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---    X-direction vertices: cylindrical coordinates  ---
!
        IF( ICS.EQ.2 ) THEN
          IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1) .OR.
     &      IJFLD.EQ.1 .OR. JKFLD.EQ.1 .OR. KIFLD.EQ.1 .OR. 
     &      ISLC(63).EQ.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(M,N)*COS(YE(M,N))
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 .AND. JFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MXY(M),N)*COS(YE(MXY(M),N))
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 .AND. KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MXZ(M),N)*COS(YE(MXZ(M),N))
              ENDDO
            ENDDO
          ELSEIF( JFLD.GT.1 .AND. KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MYZ(M),N)*COS(YE(MYZ(M),N))
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MX(M),N)*COS(YE(MX(M),N))
              ENDDO
            ENDDO
          ELSEIF( JFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MY(M),N)*COS(YE(MY(M),N))
              ENDDO
            ENDDO
          ELSEIF( KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MZ(M),N)*COS(YE(MZ(M),N))
              ENDDO
            ENDDO
          ENDIF
!
!---    X-direction vertices: non-cylindrical coordinates  ---
!
        ELSE
          IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1) .OR.
     &      ICS.EQ.3 .OR. ISLC(63).EQ.1 .OR.
     &      IJFLD.EQ.1 .OR. JKFLD.EQ.1 .OR. KIFLD.EQ.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(M,N)
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 .AND. JFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MXY(M),N)
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 .AND. KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MXZ(M),N)
              ENDDO
            ENDDO
          ELSEIF( JFLD.GT.1 .AND. KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MYZ(M),N)
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MX(M),N)
              ENDDO
            ENDDO
          ELSEIF( JFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MY(M),N)
              ENDDO
            ENDDO
          ELSEIF( KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MZ(M),N)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
!
!---    Write x-direction vertices to binary plot file  ---
!
!        OFFSET = IOFFSET + NFC_L*NVTC*NBYTR
!        IOFFSET = IOFFSET + NFC_G*NVTC*NBYTR
        OFFSET = NFC_L*NVTC*NBYTR
        CALL MPI_FILE_WRITE_AT( IPLX,OFFSET,VARX,NC,MPI_REAL8,
     &    STATUS,IERR )
!
!---    Y-direction vertices: cylindrical coordinates  ---
!
        IF( ICS.EQ.2 ) THEN
          IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1) .OR.
     &      IJFLD.EQ.1 .OR. JKFLD.EQ.1 .OR. KIFLD.EQ.1 .OR. 
     &      ISLC(63).EQ.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(M,N)*SIN(YE(M,N))
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 .AND. JFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MXY(M),N)*SIN(YE(MXY(M),N))
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 .AND. KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MXZ(M),N)*SIN(YE(MXZ(M),N))
              ENDDO
            ENDDO
          ELSEIF( JFLD.GT.1 .AND. KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MYZ(M),N)*SIN(YE(MYZ(M),N))
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MX(M),N)*SIN(YE(MX(M),N))
              ENDDO
            ENDDO
          ELSEIF( JFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MY(M),N)*SIN(YE(MY(M),N))
              ENDDO
            ENDDO
          ELSEIF( KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*XE(MZ(M),N)*SIN(YE(MZ(M),N))
              ENDDO
            ENDDO
          ENDIF
!
!---    Y-direction vertices: non-cylindrical coordinates  ---
!
        ELSE
          IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1) .OR.
     &      ICS.EQ.3 .OR. ISLC(63).EQ.1 .OR.
     &      IJFLD.EQ.1 .OR. JKFLD.EQ.1 .OR. KIFLD.EQ.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*YE(M,N)
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 .AND. JFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*YE(MXY(M),N)
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 .AND. KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*YE(MXZ(M),N)
              ENDDO
            ENDDO
          ELSEIF( JFLD.GT.1 .AND. KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*YE(MYZ(M),N)
              ENDDO
            ENDDO
          ELSEIF( IFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*YE(MX(M),N)
              ENDDO
            ENDDO
          ELSEIF( JFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*YE(MY(M),N)
              ENDDO
            ENDDO
          ELSEIF( KFLD.GT.1 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              DO M = 1,NVTC
                NC = NC + 1
                VARX(NC) = CNVLN*YE(MZ(M),N)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
!
!---    Write y-direction vertices to binary plot file  ---
!
!        OFFSET = IOFFSET + NFC_L*NVTC*NBYTR
!        IOFFSET = IOFFSET + NFC_G*NVTC*NBYTR
        OFFSET = NFC_L*NVTC*NBYTR
        CALL MPI_FILE_WRITE_AT( IPLY,OFFSET,VARX,NC,MPI_REAL8,
     &    STATUS,IERR )
!
!---    Z-direction vertices  ---
!
        IF( IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1 .OR.
     &      ICS.EQ.3 .OR. ISLC(63).EQ.1 .OR.
     &      IJFLD.EQ.1 .OR. JKFLD.EQ.1 .OR. KIFLD.EQ.1 ) THEN
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            DO M = 1,NVTC
              NC = NC + 1
              VARX(NC) = CNVLN*ZE(M,N)
            ENDDO
          ENDDO
        ELSEIF( IFLD.GT.1 .AND. JFLD.GT.1 ) THEN
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            DO M = 1,NVTC
              NC = NC + 1
              VARX(NC) = CNVLN*ZE(MXY(M),N)
            ENDDO
          ENDDO
        ELSEIF( IFLD.GT.1 .AND. KFLD.GT.1 ) THEN
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            DO M = 1,NVTC
              NC = NC + 1
              VARX(NC) = CNVLN*ZE(MXZ(M),N)
            ENDDO
          ENDDO
        ELSEIF( JFLD.GT.1 .AND. KFLD.GT.1 ) THEN
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            DO M = 1,NVTC
              NC = NC + 1
              VARX(NC) = CNVLN*ZE(MYZ(M),N)
            ENDDO
          ENDDO
        ELSEIF( IFLD.GT.1 ) THEN
          OFFSET = IOFFSET + NFC_L*NVTC*NBYTR
          IOFFSET = IOFFSET + NFC_G*NVTC*NBYTR
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            DO M = 1,NVTC
              NC = NC + 1
              VARX(NC) = CNVLN*ZE(MX(M),N)
            ENDDO
          ENDDO
        ELSEIF( JFLD.GT.1 ) THEN
          OFFSET = IOFFSET + NFC_L*NVTC*NBYTR
          IOFFSET = IOFFSET + NFC_G*NVTC*NBYTR
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            DO M = 1,NVTC
              NC = NC + 1
              VARX(NC) = CNVLN*ZE(MY(M),N)
            ENDDO
          ENDDO
        ELSEIF( KFLD.GT.1 ) THEN
          OFFSET = IOFFSET + NFC_L*NVTC*NBYTR
          IOFFSET = IOFFSET + NFC_G*NVTC*NBYTR
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            DO M = 1,NVTC
              NC = NC + 1
              VARX(NC) = CNVLN*ZE(MZ(M),N)
            ENDDO
          ENDDO
        ENDIF
!
!---    Write z-direction vertices to binary plot file  ---
!
!        OFFSET = IOFFSET + NFC_L*NVTC*NBYTR
!        IOFFSET = IOFFSET + NFC_G*NVTC*NBYTR
        OFFSET = NFC_L*NVTC*NBYTR
        CALL MPI_FILE_WRITE_AT( IPLZ,OFFSET,VARX,NC,MPI_REAL8,
     &    STATUS,IERR )
!
!---    Deallocate temporary memory for vertices  ---
!
        DEALLOCATE( VARX,STAT=ISTAT )
        CHMSG = 'VARX'
        CALL DEALLOC_ERROR( CHMSG,ISTAT )
      ENDIF
!
!---  Allocate temporary memory for plot file variables  ---
!
      ALLOCATE( VARX(NFC(ID+1)),STAT=ISTAT )
      CHMSG = 'VARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Grid-cell volume  ---
!
      OFFSET = IOFFSET + NFC_L*NBYTR
      IOFFSET = IOFFSET + NFC_G*NBYTR
      NC = 0
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for ghost cells only  ---
!
        IF( IGHC(N).EQ.1 ) CYCLE
        NC = NC + 1
        VARX(NC) = (CNVLN**3)*VOL(N)
      ENDDO
      CALL MPI_FILE_WRITE_AT( IPL,OFFSET,VARX,NC,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Allocate temporary memory for dummy integer array  ---
!
      ALLOCATE( IVARX(1:NFC(ID+1)),STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Node map  ---
!
      OFFSET = IOFFSET + NFC_L*NBYTI
      IOFFSET = IOFFSET + NFC_G*NBYTI
      NC = 0
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for ghost cells only  ---
!
        IF( IGHC(N).EQ.1 ) CYCLE
        NC = NC + 1
        IVARX(NC) = IXP(N)
      ENDDO
      CALL MPI_FILE_WRITE_AT( IPL,OFFSET,IVARX,NC,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Deallocate temporary memory for dummy integer array  ---
!
      DEALLOCATE( IVARX,STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for plot variable names  ---
!
      ALLOCATE( PLTVNMX(1:NVPLOT),STAT=ISTAT )
      CHMSG = 'PLTVNMX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Index for solute, conservation-component species, and
!     kinetic-component species plot output ---
!
      INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)
!
!---  Loop over number of plot file variables  ---
!
      DO NV = 1,NVPLOT
        IPNV = IPLOT(NV)
        IF( IPNV.EQ.200 ) CYCLE
!
!---    Aqueous pressure, Pa  ---
!
        IF( IPNV.EQ.1 ) THEN
          PLTVNMX = 'Aqueous Pressure, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*(PL(2,N)+PATM)
          ENDDO
!          PRINT *,'Aqueous pressure, Pa: ID = ',ID
!
!---    Gas pressure, Pa  ---
!
        ELSEIF( IPNV.EQ.2 ) THEN
          PLTVNMX = 'Gas Pressure, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*(PG(2,N)+PATM)
          ENDDO
!          PRINT *,'Gas pressure, Pa: ID = ',ID
!
!---    Temperature, ˚C  ---
!
        ELSEIF( IPNV.EQ.4 ) THEN
          PLTVNMX = 'Temperature, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            IF( INT(CNVPLOT(NV)).EQ.-2 ) THEN
              VARX(NC) = (T(2,N)+273.15D+0)
            ELSEIF( INT(CNVPLOT(NV)).EQ.-3 ) THEN
              VARX(NC) = (T(2,N)*1.8D+0+3.2D+1)
            ELSEIF( INT(CNVPLOT(NV)).EQ.-4 ) THEN
              VARX(NC) = (T(2,N)*1.8D+0+4.92D+2)
            ELSE
              VARX(NC) = T(2,N)
            ENDIF
          ENDDO
!
!---    Phase condition  ---
!
        ELSEIF( IPNV.EQ.5 ) THEN
          PLTVNMX = 'Phase Condition'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = REAL(NPHAZ(2,N))
          ENDDO
!
!---    Aqueous gauge pressure, Pa  ---
!
        ELSEIF( IPNV.EQ.6 ) THEN
          PLTVNMX = 'Aqueous Gauge Pressure, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*PL(2,N)
          ENDDO
!
!---    Gas gauge pressure, Pa  ---
!
        ELSEIF( IPNV.EQ.7 ) THEN
          PLTVNMX = 'Gas Gauge Pressure, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*PG(2,N)
          ENDDO
!
!---    Apparent aqueous saturation  ---
!
        ELSEIF( IPNV.EQ.9 ) THEN
          PLTVNMX = 'Apparent Aqueous Saturation'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = ASL(N)
          ENDDO
!
!---    Aqueous saturation  ---
!
        ELSEIF( IPNV.EQ.11 ) THEN
          PLTVNMX = 'Aqueous Saturation'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = SL(2,N)
          ENDDO
!
!---    Gas saturation  ---
!
        ELSEIF( IPNV.EQ.12 ) THEN
          PLTVNMX = 'Gas Saturation'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = SG(2,N)
          ENDDO
!
!---    Aqueous moisture content  ---
!
        ELSEIF( IPNV.EQ.15 ) THEN
          PLTVNMX = 'Aqueous Moisture Content'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = PORD(2,N)*SL(2,N)
          ENDDO
!
!---    Effective trapped gas saturation  ---
!
        ELSEIF( IPNV.EQ.19 ) THEN
          PLTVNMX = 'Effective Trapped Gas Saturation'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = ASGT(N)
          ENDDO
!
!---    Diffusive porosity  ---
!
        ELSEIF( IPNV.EQ.20 ) THEN
          PLTVNMX = 'Diffusive Porosity'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = PORD(2,N)
          ENDDO
!          PRINT *,'Diffusive porosity, Pa: ID = ',ID
!
!---    Aqueous water mass fraction  ---
!
        ELSEIF( IPNV.EQ.24 ) THEN
          PLTVNMX = 'Aqueous Water Mass Fraction'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = XLW(2,N)
          ENDDO
!
!---    Aqueous hydraulic head, m  ---
!
        ELSEIF( IPNV.EQ.27 ) THEN
          PLTVNMX = 'Aqueous Hydraulic Head ' //
     &    '(Fresh Water),' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*(PL(2,N)/RHORL/GRAV + ZP(N))
          ENDDO
!
!---    Gas hydraulic head, m  ---
!
        ELSEIF( IPNV.EQ.28 ) THEN
          PLTVNMX = 'Gas Hydraulic Head ' //
     &    '(Fresh Water),' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*(PG(2,N)/RHORL/GRAV + ZP(N))
          ENDDO
!
!---    Rock/soil type  ---
!
        ELSEIF( IPNV.EQ.30 ) THEN
          PLTVNMX = 'Rock/Soil Type'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = REAL(IZ(N))
          ENDDO
!
!---    Aqueous relative permeability  ---
!
        ELSEIF( IPNV.EQ.31 ) THEN
          PLTVNMX = 'Aqueous Relative Permeability'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = RKL(2,N)
          ENDDO
!
!---    Aqueous density, kg/m^3 aqu  ---
!
        ELSEIF( IPNV.EQ.34 ) THEN
          PLTVNMX = 'Aqueous Density, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*RHOL(2,N)
          ENDDO
!
!---    Total water mass, kg  ---
!
        ELSEIF( IPNV.EQ.37 ) THEN
          PLTVNMX = 'Total Water Mass, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*PORD(2,N)*VOL(N)*
     &        XLW(2,N)*SL(2,N)*RHOL(2,N)
          ENDDO
!
!---    Water mass source integral, kg  ---
!
        ELSEIF( IPNV.EQ.40 ) THEN
          PLTVNMX = 'Water Mass Source Integral, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*SRCIW(N)
          ENDDO
!
!---    Aqueous Courant number  ---
!
        ELSEIF( IPNV.EQ.49 ) THEN
          PLTVNMX = 'Maximum Local Aqueous Courant Number '
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CRNTL(N)
          ENDDO
!
!---    Webb Matching Point Head, m  ---
!
        ELSEIF( IPNV.EQ.63 ) THEN
          PLTVNMX = 'Webb Matching Point Head, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*SCHR(9,N)
          ENDDO
!
!---    Total Stress-xx, Effective Stress-xx, or ice pressure Pa  ---
!
        ELSEIF( IPNV.EQ.80 ) THEN
!
!---    Geomechanics  ---
!
          IF( ISLC(50).NE.0 ) THEN
            ISLC50X = ABS(ISLC(50))
            IMODX = MOD(ISLC50X,10)
            IF( IPNVGC.EQ.0 ) THEN
              PLTVNMX = 'Total Stress-xx, ' // UNPLOT(NV)
            ELSE
              PLTVNMX = 'Effective Stress-xx, ' // UNPLOT(NV)
            ENDIF
          ELSE
            PLTVNMX = 'Ice Pressure, ' // UNPLOT(NV)
          ENDIF
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
!
!---        Geomechanics  ---
!
            IF( ISLC(50).NE.0 ) THEN
              IF( IPNVGC.EQ.0 ) THEN
                VARX(NC) = SIG_GM(1,N)
!
!---            Poroelasticity ---
!
                IF( IMODX.NE.2 .AND. IMODX.NE.4 ) THEN
                  PX = MAX( PG(2,N),PL(2,N) ) + PATM
                  VARX(NC) = VARX(NC) + PROP_GM(3,N)*(PX-PCMP(N))
                ENDIF
!
!---            Thermoelasticity ---
!
                IF( IMODX.NE.3 .AND. IMODX.NE.4 ) THEN
                 BLK3X = PROP_GM(1,N)/(1.D+0-2.D+0*PROP_GM(2,N))
                 VARX(NC) = VARX(NC) + PROP_GM(4,N)*
     &             BLK3X*(T(2,N)-TCMP(N))
                ENDIF
                VARX(NC) = CNVPLOT(NV)*VARX(NC)
              ELSE
                VARX(NC) = CNVPLOT(NV)*SIG_GM(1,N)
              ENDIF
            ELSE
              VARX(NC) = 0.D+0
            ENDIF
          ENDDO
!
!---    X-direction node-centered aqueous flux, m/s  ---
!
        ELSEIF( IPNV.EQ.87 ) THEN
          PLTVNMX = 'X-Dir. Aqueous Darcy Velocity ' //
     &    '(Node Centered), ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*(UL(1,1,N)+UL(1,2,N))
          ENDDO
!          PRINT *,'X-direction node-centered aqueous flux: ID = ',ID
!
!---    Y-direction node-centered aqueous flux, m/s  ---
!
        ELSEIF( IPNV.EQ.88 ) THEN
          PLTVNMX = 'Y-Dir. Aqueous Darcy Velocity ' //
     &    '(Node Centered), ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*(VL(1,1,N)+VL(1,2,N))
          ENDDO
!          PRINT *,'Y-direction node-centered aqueous flux: ID = ',ID
!
!---    Z-direction node-centered aqueous flux, m/s  ---
!
        ELSEIF( IPNV.EQ.89 ) THEN
          PLTVNMX = 'Z-Dir. Aqueous Darcy Velocity ' //
     &    '(Node Centered), ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*(WL(1,1,N)+WL(1,2,N))
          ENDDO
!          PRINT *,'Z-direction node-centered aqueous flux: ID = ',ID
!
!---    Node number  ---
!
        ELSEIF( IPNV.EQ.100 ) THEN
          PLTVNMX = 'Node Number'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = REAL(ND(N))
          ENDDO
!
!---    Total Stress-yy or Effective Stress-yy, Pa  ---
!
        ELSEIF( IPNV.EQ.101 ) THEN
          ISLC50X = ABS(ISLC(50))
          IMODX = MOD(ISLC50X,10)
          IF( IPNVGC.EQ.0 ) THEN
            PLTVNMX = 'Total Stress-yy, ' // UNPLOT(NV)
          ELSE
            PLTVNMX = 'Effective Stress-yy, ' // UNPLOT(NV)
          ENDIF
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            IF( IPNVGC.EQ.0 ) THEN
              VARX(NC) = SIG_GM(2,N)
!
!---          Poroelasticity ---
!
              IF( IMODX.NE.2 .AND. IMODX.NE.4 ) THEN
                PX = MAX( PG(2,N),PL(2,N) ) + PATM
                VARX(NC) = VARX(NC) + PROP_GM(3,N)*(PX-PCMP(N))
              ENDIF
!
!---          Thermoelasticity ---
!
              IF( IMODX.NE.3 .AND. IMODX.NE.4 ) THEN
               BLK3X = PROP_GM(1,N)/(1.D+0-2.D+0*PROP_GM(2,N))
               VARX(NC) = VARX(NC) + PROP_GM(4,N)*BLK3X*(T(2,N)-TCMP(N))
              ENDIF
              VARX(NC) = CNVPLOT(NV)*VARX(NC)
            ELSE
              VARX(NC) = CNVPLOT(NV)*SIG_GM(2,N)
            ENDIF
          ENDDO
!
!---    Trapped gas saturation  ---
!
        ELSEIF( IPNV.EQ.105 ) THEN
          PLTVNMX = 'Trapped Gas Saturation'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = SGT(2,N)
          ENDDO
!
!---    Minimum effective aqueous saturation  ---
!
        ELSEIF( IPNV.EQ.127 ) THEN
          PLTVNMX = 'Minimum Effective Aqueous Saturation'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = ASLMIN(2,N)
          ENDDO
!
!---    Total Stress-zz, Effective Stress-zz, or CH4 vapor partial 
!       pressure, Pa  ---
!
        ELSEIF( IPNV.EQ.130 ) THEN
!
!---      Geomechanics  ---
!
          IF( ISLC(50).NE.0 ) THEN
            ISLC50X = ABS(ISLC(50))
            IMODX = MOD(ISLC50X,10)
            IF( IPNVGC.EQ.0 ) THEN
              PLTVNMX = 'Total Stress-zz, ' // UNPLOT(NV)
            ELSE
              PLTVNMX = 'Effective Stress-zz, ' // UNPLOT(NV)
            ENDIF
            PLTVNMX = 'CH4 Vapor Partial Pressure, ' // UNPLOT(NV)
          ENDIF
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
!
!---        Geomechanics  ---
!
            IF( ISLC(50).NE.0 ) THEN
              IF( IPNVGC.EQ.0 ) THEN
                VARX(NC) = SIG_GM(3,N)
!
!---            Poroelasticity ---
!
                IF( IMODX.NE.2 .AND. IMODX.NE.4 ) THEN
                  PX = MAX( PG(2,N),PL(2,N) ) + PATM
                  VARX(NC) = VARX(NC) + PROP_GM(3,N)*(PX-PCMP(N))
                ENDIF
!
!---            Thermoelasticity ---
!
                IF( IMODX.NE.3 .AND. IMODX.NE.4 ) THEN
                 BLK3X = PROP_GM(1,N)/(1.D+0-2.D+0*PROP_GM(2,N))
                 VARX(NC) = VARX(NC) + 
     &             PROP_GM(4,N)*BLK3X*(T(2,N)-TCMP(N))
                ENDIF
                VARX(NC) = CNVPLOT(NV)*VARX(NC)
              ELSE
                VARX(NC) = CNVPLOT(NV)*SIG_GM(3,N)
              ENDIF
            ELSE
              VARX(NC) = 0.D+0
            ENDIF
          ENDDO
!
!---    Water mass source rate, kg/s  ---
!
        ELSEIF( IPNV.EQ.140 ) THEN
          PLTVNMX = 'Water Mass Source Rate, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*SRCW(2,N)
          ENDDO
!
!---    Webb matching point saturation  ---
!
        ELSEIF( IPNV.EQ.150 ) THEN
          PLTVNMX = 'Webb Matching Point Saturation'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = SCHR(8,N)
          ENDDO
!
!---    Aqueous viscosity, Pa s  ---
!
        ELSEIF( IPNV.EQ.176 ) THEN
          PLTVNMX = 'Aqueous Viscosity, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*VISL(2,N)
          ENDDO
!
!---    X-direction aqueous relative permeability  ---
!
        ELSEIF( IPNV.EQ.201 ) THEN
          PLTVNMX = 'X-Direction Aqueous Relative Permeability'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = RKL(2,N)
          ENDDO
!
!---    Y-direction aqueous relative permeability  ---
!
        ELSEIF( IPNV.EQ.201 ) THEN
          PLTVNMX = 'Y-Direction Aqueous Relative Permeability'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = RKL(2,N)
          ENDDO
!
!---    Z-direction aqueous relative permeability  ---
!
        ELSEIF( IPNV.EQ.201 ) THEN
          PLTVNMX = 'Z-Direction Aqueous Relative Permeability'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = RKL(2,N)
          ENDDO
!
!---    Stress-yz, Pa  ---
!
        ELSEIF( IPNV.EQ.214 ) THEN
          PLTVNMX = 'Stress-yz, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*SIG_GM(4,N)
          ENDDO
!
!---    Stress-xz, Pa  ---
!
        ELSEIF( IPNV.EQ.224 ) THEN
          PLTVNMX = 'Stress-xz, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*SIG_GM(5,N)
          ENDDO
!
!---    Stress-xy, Pa  ---
!
        ELSEIF( IPNV.EQ.226 ) THEN
          PLTVNMX = 'Stress-xy, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*SIG_GM(6,N)
          ENDDO
!
!---    X-direction permeability, m^2  ---
!
        ELSEIF( IPNV.EQ.247 ) THEN
          PLTVNMX = 'X-Direction Intrinsic Permeability, ' // 
     &      UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*PERM(1,N)*PERMRF(2,N)
          ENDDO
!
!---    Y-direction permeability, m^2  ---
!
        ELSEIF( IPNV.EQ.248 ) THEN
          PLTVNMX = 'Y-Direction Intrinsic Permeability, ' // 
     &      UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*PERM(2,N)*PERMRF(2,N)
          ENDDO
!
!---    Z-direction permeability, m^2  ---
!
        ELSEIF( IPNV.EQ.249 ) THEN
          PLTVNMX = 'Z-Direction Intrinsic Permeability, ' // 
     &      UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*PERM(3,N)*PERMRF(2,N)
          ENDDO
!
!---    X-direction grid-cell centroid, m ---
!
        ELSEIF( IPNV.EQ.291 ) THEN
          PLTVNMX = 'X Node-Centroid Position, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*XP(N)
          ENDDO
!
!---    Y-direction grid-cell centroid, m ---
!
        ELSEIF( IPNV.EQ.292 ) THEN
          PLTVNMX = 'Y Node-Centroid Position, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*YP(N)
          ENDDO
!
!---    Z-direction grid-cell centroid, m ---
!
        ELSEIF( IPNV.EQ.293 ) THEN
          PLTVNMX = 'Z Node-Centroid Position, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*ZP(N)
          ENDDO
!
!---    Similarity variable, m^2/s  ---
!
        ELSEIF( IPNV.EQ.299 ) THEN
          PLTVNMX = 'Similarity Variable, ' // UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = CNVPLOT(NV)*(XP(N)**2)/(TM+SMALL)
          ENDDO
          CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---    Strain-xx  ---
!
        ELSEIF( IPNV.EQ.369 ) THEN
          PLTVNMX = 'Strain-xx'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = EPS_GM(1,N)
          ENDDO
!
!---    Strain-yy  ---
!
        ELSEIF( IPNV.EQ.370 ) THEN
          PLTVNMX = 'Strain-yy'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = EPS_GM(2,N)
          ENDDO
!
!---    Strain-zz  ---
!
        ELSEIF( IPNV.EQ.371 ) THEN
          PLTVNMX = 'Strain-zz'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = EPS_GM(3,N)
          ENDDO
!
!---    Strain-yz  ---
!
        ELSEIF( IPNV.EQ.372 ) THEN
          PLTVNMX = 'Strain-yz'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = EPS_GM(4,N)
          ENDDO
!
!---    Strain-xz  ---
!
        ELSEIF( IPNV.EQ.373 ) THEN
          PLTVNMX = 'Strain-xz'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = EPS_GM(5,N)
          ENDDO
!
!---    Strain-xy  ---
!
        ELSEIF( IPNV.EQ.374 ) THEN
          PLTVNMX = 'Strain-xy'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            VARX(NC) = EPS_GM(6,N)
          ENDDO
!!
!!---    X displacement, m  ---
!!
!        ELSEIF( IPNV.EQ.375 ) THEN
!          IF( IPNVGC.EQ.0 ) THEN
!            PLTVNMX = 'X-Displacement, ' // UNPLOT(NV)
!          ELSE
!            PLTVNMX = 'X-Displacement: FE-Node: ' // UNPLOT(NV)
!          ENDIF
!          NC = 0
!          DO N = 1,NFCGC(ID+1)
!!
!!---        Skip for ghost cells only  ---
!!
!            IF( IGHC(N).EQ.1 ) CYCLE
!            NC = NC + 1
!            IF( IPNVGC.EQ.0 ) THEN
!              VARX(NC) = 0.D+0
!              DO L = 1,8
!                NFEN = ND_GM(L,N)
!                VARX(NC) = VARX(NC) + CNVPLOT(NV)*
!     &            (U_GM(2,NFEN)-U_GM(1,NFEN))
!              ENDDO
!            ELSE
!              NFEN = ND_GM(IPNVGC,N)
!              VARX(NC) = CNVPLOT(NV)*(U_GM(2,NFEN)-U_GM(1,NFEN))
!            ENDIF
!          ENDDO
!!
!!---    Y displacement, m  ---
!!
!        ELSEIF( IPNV.EQ.376 ) THEN
!          IF( IPNVGC.EQ.0 ) THEN
!            PLTVNMX = 'Y-Displacement, ' // UNPLOT(NV)
!          ELSE
!            PLTVNMX = 'Y-Displacement: FE-Node: ' // UNPLOT(NV)
!          ENDIF
!          NC = 0
!          DO N = 1,NFCGC(ID+1)
!!
!!---        Skip for ghost cells only  ---
!!
!            IF( IGHC(N).EQ.1 ) CYCLE
!            NC = NC + 1
!            IF( IPNVGC.EQ.0 ) THEN
!              VARX(NC) = 0.D+0
!              DO L = 1,8
!                NFEN = ND_GM(L,N)
!                VARX(NC) = VARX(NC) + CNVPLOT(NV)*
!     &            (V_GM(2,NFEN)-V_GM(1,NFEN))
!              ENDDO
!            ELSE
!              NFEN = ND_GM(IPNVGC,N)
!              VARX(NC) = CNVPLOT(NV)*(V_GM(2,NFEN)-V_GM(1,NFEN))
!            ENDIF
!          ENDDO
!!
!!---    Z displacement, m  ---
!!
!        ELSEIF( IPNV.EQ.377 ) THEN
!          IF( IPNVGC.EQ.0 ) THEN
!            PLTVNMX = 'Z-Displacement, ' // UNPLOT(NV)
!          ELSE
!            PLTVNMX = 'Z-Displacement: FE-Node: ' // UNPLOT(NV)
!          ENDIF
!          NC = 0
!          DO N = 1,NFCGC(ID+1)
!!
!!---        Skip for ghost cells only  ---
!!
!            IF( IGHC(N).EQ.1 ) CYCLE
!            NC = NC + 1
!            IF( IPNVGC.EQ.0 ) THEN
!              VARX(NC) = 0.D+0
!              DO L = 1,8
!                NFEN = ND_GM(L,N)
!                VARX(NC) = VARX(NC) + CNVPLOT(NV)*
!     &            (W_GM(2,NFEN)-W_GM(1,NFEN))
!              ENDDO
!            ELSE
!              NFEN = ND_GM(IPNVGC,N)
!              VARX(NC) = CNVPLOT(NV)*(W_GM(2,NFEN)-W_GM(1,NFEN))
!            ENDIF
!          ENDDO
!
!---    X-direction node index  ---
!
        ELSEIF( IPNV.EQ.388 ) THEN
          PLTVNMX = 'I Node Index'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            IDX = MOD( ND(N),IFLD )
            VARX(NC) = REAL(IDX)
          ENDDO
!
!---    Y-direction node index  ---
!
        ELSEIF( IPNV.EQ.389 ) THEN
          PLTVNMX = 'J Node Index'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            JDX = MOD( ND(N),IFLD*JFLD )/IFLD + 1
            VARX(NC) = REAL(JDX)
          ENDDO
!
!---    Z-direction node index  ---
!
        ELSEIF( IPNV.EQ.390 ) THEN
          PLTVNMX = 'K Node Index'
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            KDX = ND(N)/(IFLD*JFLD) + 1
            VARX(NC) = REAL(KDX)
          ENDDO
!
!---    X-direction node-centered surface area  ---
!
        ELSEIF( IPNV.EQ.391 ) THEN
          PLTVNMX = 'X-Dir. Surface Area (Node Centered), ' // 
     &      UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            IDX = MOD( ND(N),IFLD )
            VARX(NC) = CNVPLOT(NV)*0.5D+0*(AFX(1,N)+AFX(2,N))
          ENDDO
!
!---    Y-direction node-centered surface area  ---
!
        ELSEIF( IPNV.EQ.392 ) THEN
          PLTVNMX = 'Y-Dir. Surface Area (Node Centered), ' // 
     &      UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            IDX = MOD( ND(N),IFLD )
            VARX(NC) = CNVPLOT(NV)*0.5D+0*(AFY(1,N)+AFY(2,N))
          ENDDO
!
!---    Z-direction node-centered surface area  ---
!
        ELSEIF( IPNV.EQ.393 ) THEN
          PLTVNMX = 'Z-Dir. Surface Area (Node Centered), ' // 
     &      UNPLOT(NV)
          NC = 0
          DO N = 1,NFCGC(ID+1)
!
!---        Skip for ghost cells only  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
            NC = NC + 1
            IDX = MOD( ND(N),IFLD )
            VARX(NC) = CNVPLOT(NV)*0.5D+0*(AFZ(1,N)+AFZ(2,N))
          ENDDO
!
!---    Solute, conservation-component species, and
!       kinetic-component species plot output ---
!
        ELSEIF( IPNV.GT.400 .AND. IPNV.LE.INDX ) THEN
          NSL = ((IPNV-400)/33) + 1
          IF( NSL.GT.NSOLU ) IUNMOL = 1
          IPNVX = MOD((IPNV-400),33)
          IDB = INDEX( SOLUT(NSL)(1:),'  ') - 1
!
!---      401 Solute volumetric concentration, 1/m^3  ---
!---      401 Total_species volumetric concentration, mol/m^3  ---
!
          IF( IPNVX.EQ.1 ) THEN
!
!---        401 Total_species volumetric concentration, mol/m^3  ---
!
            IF( NSL.GT.NSOLU )THEN
              NEQ = NSL-NSOLU
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = 0.D+0
                VARX(NC) = CNVPLOT(NV)*0.5D+0*(AFZ(1,N)+AFZ(2,N))
                DO M = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(M+1,NEQ)
                  IF( NSP.LE.NSPL ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NC) = VARX(NC) + EQ_C(M,NEQ)*SP_CX
                  ENDIF
                ENDDO
                VARX(NC) = CNVPLOT(NV)*VARX(NC)*1.D-3
              ENDDO
!
!---        401 Solute volumetric concentration, 1/m^3  ---
!
            ELSE
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = CNVPLOT(NV)*C(N,NSL)
              ENDDO
            ENDIF
!
!---      402 Solute aqueous concentration, 1/m^3  ---
!---      402 Total_species aqueous concentration, mol/m^3  ---
!
          ELSEIF( IPNVX.EQ.2 ) THEN
!
!---        402 Total_species aqueous concentration, mol/m^3  ---
!
            IF( NSL.GT.NSOLU )THEN
              NEQ = NSL-NSOLU
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = 0.D+0
                VARX(NC) = CNVPLOT(NV)*0.5D+0*(AFZ(1,N)+AFZ(2,N))
                DO M = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(M+1,NEQ)
                  IF( NSP.LE.NSPL ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NC) = VARX(NC) + EQ_C(M,NEQ)*SP_CX
                  ENDIF
                ENDDO
                VARX(NC) = CNVPLOT(NV)*VARX(NC)*YL(N,NSL)*1.D-3/
     &            (SL(2,N)*PORD(2,N)+SMALL)
              ENDDO
!
!---        402 Solute aqueous concentration, 1/m^3  ---
!
            ELSE
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = CNVPLOT(NV)*C(N,NSL)*YL(N,NSL)/
     &            (SL(2,N)*PORD(2,N)+SMALL)
              ENDDO
            ENDIF
!
!---      405 Solute/species aqueous mole fraction  ---
!
          ELSEIF( IPNVX.EQ.5 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              NC = NC + 1
              VARX(NC) = YL(N,NSL)
            ENDDO
!
!---      407 Solute/species inventory  ---
!
          ELSEIF( IPNVX.EQ.6 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              NC = NC + 1
              VARX(NC) = YN(N,NSL)
            ENDDO
!
!---      408 X-dir. solute flux (node centered) 1/m^2 s  ---
!---      408 X-dir. total_species flux (node centered) mol/m^2 s  ---
!
          ELSEIF( IPNVX.EQ.8 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              NC = NC + 1
              VARX(NC) = CNVPLOT(NV)*5.D-1*(UC(1,N,NSL)+UC(2,N,NSL))
            ENDDO
!
!---      409 Y-dir. solute flux (node centered) 1/m^2 s  ---
!---      409 Y-dir. total_species flux (node centered) mol/m^2 s  ---
!
          ELSEIF( IPNVX.EQ.9 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              NC = NC + 1
              VARX(NC) = CNVPLOT(NV)*5.D-1*(VC(1,N,NSL)+VC(2,N,NSL))
            ENDDO
!
!---      410 Z-dir. solute flux (node centered) 1/m^2 s  ---
!---      410 Z-dir. total_species flux (node centered) mol/m^2 s  ---
!
          ELSEIF( IPNVX.EQ.10 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              NC = NC + 1
              VARX(NC) = CNVPLOT(NV)*5.D-1*(WC(1,N,NSL)+WC(2,N,NSL))
            ENDDO
!
!---      411 Solute source integral, 1/s  ---
!---      411 Total-species source integral, mol/s  ---
!
          ELSEIF( IPNVX.EQ.11 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              NC = NC + 1
              VARX(NC) = CNVPLOT(NV)*SRCIC(N,NSL)
            ENDDO
!
!---      412 Solute (unused)  ---
!---      412 Total-species exchange volumetric conc., mol/m^3  ---
!
          ELSEIF( IPNVX.EQ.12 ) THEN
!
!---        412 Total-species exchange volumetric conc., mol/m^3  ---
!
            IF( NSL.GT.NSOLU )THEN
              NEQ = NSL-NSOLU
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = 0.D+0
                DO M = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(M+1,NEQ)
                  IF( NSP.GT.NSPL+NSPS.AND.NSP.LE.NSPL+NSPS+NSPE ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NC) = VARX(NC) + EQ_C(M,NEQ)*SP_CX
                  ENDIF
                ENDDO
                VARX(NC) = CNVPLOT(NV)*VARX(NC)*1.D-3
              ENDDO
!
!---        412 Solute (unused)  ---
!
            ELSE
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = 0.D+0
              ENDDO
            ENDIF
!
!---      413 Solute (unused)  ---
!---      413 Total-species exchange aqueous conc., mol/m^3  ---
!
          ELSEIF( IPNVX.EQ.13 ) THEN
!
!---        413 Total-species exchange aqueous conc., mol/m^3  ---
!
            IF( NSL.GT.NSOLU )THEN
              NEQ = NSL-NSOLU
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = 0.D+0
                DO M = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(M+1,NEQ)
                  IF( NSP.GT.NSPL+NSPS.AND.NSP.LE.NSPL+NSPS+NSPE ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NC) = VARX(NC) + EQ_C(M,NEQ)*SP_CX
                  ENDIF
                ENDDO
                VARX(NC) = CNVPLOT(NV)*VARX(NC)*1.D-3*YL(N,NSL)/
     &            (SL(2,N)*PORD(2,N)+SMALL)
              ENDDO
!
!---        413 Solute (unused)  ---
!
            ELSE
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = 0.D+0
              ENDDO
            ENDIF
!
!---      422 Solute mass  ---
!---      422 Total-species mass, mol  ---
!
          ELSEIF( IPNVX.EQ.22 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              NC = NC + 1
              VARX(NC) = CNVPLOT(NV)*C(N,NSL)*VOL(N)
            ENDDO
!
!---      426 Solute (unused)  ---
!---      426 Total-species solid volumetric conc., mol/m^3  ---
!
          ELSEIF( IPNVX.EQ.26 ) THEN
!
!---        426 Total-species solid volumetric conc., mol/m^3  ---
!
            IF( NSL.GT.NSOLU )THEN
              NEQ = NSL-NSOLU
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = 0.D+0
                DO M = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(M+1,NEQ)
                  IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NC) = VARX(NC) + EQ_C(M,NEQ)*SP_CX
                  ENDIF
                ENDDO
                VARX(NC) = CNVPLOT(NV)*VARX(NC)*1.D-3
              ENDDO
!
!---        426 Solute (unused)  ---
!
            ELSE
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = 0.D+0
              ENDDO
            ENDIF
!
!---      427 Solute (unused)  ---
!---      427 Total-species solid aqueous conc., mol/m^3  ---
!
          ELSEIF( IPNVX.EQ.27 ) THEN
!
!---        427 Total-species solid aqueous conc., mol/m^3  ---
!
            IF( NSL.GT.NSOLU )THEN
              NEQ = NSL-NSOLU
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = 0.D+0
                DO M = 1,IEQ_C(1,NEQ)
                  NSP = IEQ_C(M+1,NEQ)
                  IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
                    IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                      SP_CX = 0.D+0
                    ELSE
                      SP_CX = SP_C(N,NSP)
                    ENDIF
                    VARX(NC) = VARX(NC) + EQ_C(M,NEQ)*SP_CX
                  ENDIF
                ENDDO
                VARX(NC) = CNVPLOT(NV)*VARX(NC)*YL(N,NSL)*1.D-3/
     &            (SL(2,N)*PORD(2,N)+SMALL)
              ENDDO
!
!---        427 Solute (unused)  ---
!
            ELSE
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = 0.D+0
              ENDDO
            ENDIF
          ENDIF
!
!---    Reactive species plot output ---
!
        ELSEIF( IPNV.GT.INDX .AND. IPNV.LE.(INDX+NSPR*33) ) THEN
          NSP = ((IPNV-INDX)/33) + 1
          IPNVX = MOD((IPNV-INDX),33)
          IF( NSP.GT.NSPL+NSPS+NSPE ) THEN
            IDB = INDEX( SPNMG(NSP-NSPL-NSPS-NSPE)(1:),'  ') - 1
            SPNMX = SPNMG(NSP-NSPL-NSPS-NSPE)
          ELSEIF( NSP.GT.NSPL+NSPS ) THEN
            IDB = INDEX( SPNME(NSP-NSPL-NSPS)(1:),'  ') - 1
            SPNMX = SPNME(NSP-NSPL-NSPS)
          ELSEIF( NSP.GT.NSPL ) THEN
            IDB = INDEX( SPNMS(NSP-NSPL)(1:),'  ') - 1
            SPNMX = SPNMS(NSP-NSPL)
          ELSE
            IDB = INDEX( SPNML(NSP)(1:),'  ') - 1
            SPNMX = SPNML(NSP)
          ENDIF
!
!---      401 Reactive species volumetric concentration, mol/m^3  ---
!
          IF( IPNVX.EQ.1 ) THEN
            IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
              NSP_M = NSP-NSPL
              IF( ISP_MN(NSP).EQ.1 ) THEN
                NC = 0
                DO N = 1,NFCGC(ID+1)
!
!---              Skip for ghost cells only  ---
!
                  IF( IGHC(N).EQ.1 ) CYCLE
                  NC = NC + 1
                  VARX(NC) = CNVPLOT(NV)*1.D-3*
     &              (SP_C(N,NSP)+SP_CMN(N,NSP_M))
                ENDDO
              ELSE
                NC = 0
                DO N = 1,NFCGC(ID+1)
!
!---              Skip for ghost cells only  ---
!
                  IF( IGHC(N).EQ.1 ) CYCLE
                  NC = NC + 1
                  VARX(NC) = CNVPLOT(NV)*1.D-3*SP_C(N,NSP)
                ENDDO
              ENDIF
            ELSE
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = CNVPLOT(NV)*1.D-3*SP_C(N,NSP)
              ENDDO
            ENDIF
!
!---      402 Reactive species aqueous concentration, mol/m^3  ---
!
          ELSEIF( IPNVX.EQ.2 ) THEN
            IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
              NSP_M = NSP-NSPL
              IF( ISP_MN(NSP).EQ.1 ) THEN
                NC = 0
                DO N = 1,NFCGC(ID+1)
!
!---              Skip for ghost cells only  ---
!
                  IF( IGHC(N).EQ.1 ) CYCLE
                  NC = NC + 1
                  VARX(NC) = CNVPLOT(NV)*1.D-3*
     &              (SP_C(N,NSP)+SP_CMN(N,NSP_M))/
     &              (SL(2,N)*PORD(2,N)+SMALL)
                ENDDO
              ELSE
                NC = 0
                DO N = 1,NFCGC(ID+1)
!
!---              Skip for ghost cells only  ---
!
                  IF( IGHC(N).EQ.1 ) CYCLE
                  NC = NC + 1
                  VARX(NC) = CNVPLOT(NV)*1.D-3*SP_C(N,NSP)/
     &              (SL(2,N)*PORD(2,N)+SMALL)
                ENDDO
              ENDIF
            ELSE
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = CNVPLOT(NV)*1.D-3*SP_C(N,NSP)/
     &            (SL(2,N)*PORD(2,N)+SMALL)
              ENDDO
            ENDIF
!
!---      411 Reactive species source integral, mol  ---
!
          ELSEIF( IPNVX.EQ.11 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              NC = NC + 1
              VARX(NC) = CNVPLOT(NV)*1.D-3*SRCIC(N,NSL)
            ENDDO
!
!---      424 Reactive species mineral area, m^2  ---
!
          ELSEIF( IPNVX.EQ.24 ) THEN
            IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
              NSP_M = NSP-NSPL
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = CNVPLOT(NV)*SP_AREA(N,NSP_M)
              ENDDO
            ENDIF
!
!---      425 Reactive species mineral rate, mol/s  ---
!
          ELSEIF( IPNVX.EQ.25 ) THEN
            IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
              NSP_M = NSP-NSPL
              NC = 0
              DO N = 1,NFCGC(ID+1)
!
!---            Skip for ghost cells only  ---
!
                IF( IGHC(N).EQ.1 ) CYCLE
                NC = NC + 1
                VARX(NC) = CNVPLOT(NV)*SP_RATE(N,NSP_M)
              ENDDO
            ENDIF
!
!---      426 Reactive species volume fraction  ---
!
          ELSEIF( IPNVX.EQ.26 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              NC = NC + 1
              VARX(NC) = RS_S(3,NSP_M,N)
            ENDDO
!
!---      427 pH  ---
!
          ELSEIF( IPNVX.EQ.26 ) THEN
            NC = 0
            DO N = 1,NFCGC(ID+1)
!
!---          Skip for ghost cells only  ---
!
              IF( IGHC(N).EQ.1 ) CYCLE
              NC = NC + 1
              VARX(NC) = C_PH(N)
            ENDDO
          ENDIF
        ENDIF
!
!---    Plot variable name and unitS  ---
!
        NVAR = 64
        OFFSET = IOFFSET
        IOFFSET = IOFFSET + NVAR
!        PRINT *,'5: OFFSET = ',OFFSET,' IOFFSET = ',IOFFSET,' ID = ',ID
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IPL,OFFSET,PLTVNMX(NV),
     &    NVAR,MPI_CHAR,STATUS,IERR )
!
!---    Plot variable  ---
!
!        PRINT *,'VARX: NV = ',NV,' NFC_L = ',NFC_L,' NFC_G = ',NFC_G,
!     &    ' ID = ',ID
        OFFSET = IOFFSET + NFC_L*NBYTR
        IOFFSET = IOFFSET + NFC_G*NBYTR
!        PRINT *,'VARX: NV = ',NV,'OFFSET = ',OFFSET,
!     &    ' IOFFSET = ',IOFFSET,' ID = ',ID
!        PRINT *,'VARX: NV = ',NV,' NC = ',NC,' ID = ',ID
        CALL MPI_FILE_WRITE_AT( IPL,OFFSET,VARX,NC,MPI_REAL8,
     &    STATUS,IERR )
      ENDDO
!
!---  Deallocate temporary memory for plot file variables  ---
!
      DEALLOCATE( VARX,STAT=ISTAT )
      CHMSG = 'VARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Deallocate temporary memory for plot variable names  ---
!
      DEALLOCATE( PLTVNMX,STAT=ISTAT )
      CHMSG = 'PLTVNMX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Close plot_xxxxxxx.bin file  ---
!
      CALL MPI_FILE_CLOSE( IPL,IERR )
!
!---  Write vertex files once  ---
!
      IF( IPLSKP.EQ.0 ) THEN
        IPLSKP = 1
!
!---    Close plot_xgrid.bin file  ---
!
        CALL MPI_FILE_CLOSE( IPLX,IERR )
!
!---    Close plot_xgrid.bin file  ---
!
        CALL MPI_FILE_CLOSE( IPLY,IERR )
!
!---    Close plot_xgrid.bin file  ---
!
        CALL MPI_FILE_CLOSE( IPLZ,IERR )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRPLOT_W group  ---
!
      RETURN
      END


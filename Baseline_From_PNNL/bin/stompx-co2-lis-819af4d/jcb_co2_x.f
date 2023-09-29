!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCB_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Call the Jacobian matrix loaders for STOMP-CO2
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 4 November 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOURC
      USE SOLTN
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
!----------------------Type Declarations-------------------------------!
!
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCB_CO2'
      ICNV = 3
!
!---  Zero residual array  ---
!
      DO N = 1,NFCGC(ID+1)
        DO M = 1,ISVC
          RSDL(M,N) = 0.D+0
        ENDDO
      ENDDO
!
!---  Load Jacobian matrix for the water equation
!     (zero flux boundary)  ---
!
      CALL JCBW_CO2
!
!---  Load Jacobian matrix for the CO2 equation
!     (zero flux boundary)  ---
!
      CALL JCBA_CO2
!
!---  Load Jacobian matrix for the salt equation
!     (zero flux boundary), isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL JCBS_CO2
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCB_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBP_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Configure the arrays for compressed sparse row matrix storage.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 March 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
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
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBP_CO2'
!
!---  Coupled-well model  ---
!
      IF( N_CW.GT.0  ) THEN
        CALL JCBP_CW_CO2
!
!---  No coupled well model  ---
!
      ELSE
        CALL JCBP_NCW_CO2
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBP_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCB_SV
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Set values of the Jacobian matrix
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 March 2022
!

!----------------------LIS Modules-------------------------------------!
!
      USE LIS_STOMP
!







!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
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
      REAL*8, DIMENSION(NNZFR) :: BUFFER
      INTEGER, DIMENSION(NNZFR) :: ICOL
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCB_SV'

!
!---  Lis solver  ---
!
      CALL lis_matrix_set_csr_f( NNZ,NLU,MLU,DLU,F_MAT,IERR )
      IF( IERR.NE.0 ) THEN
        I_ERR(1) = 0
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
        M_ERR(1) = 'Coupled Flow and Transport Lis Matrix Set CSR'
        RETURN
      ENDIF

!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCB_SV group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBP_CW_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Configure the arrays for compressed sparse row matrix storage
!     for problems with coupled wells.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 March 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE JACOB
      USE GRID
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
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBP_CW_CO2'
!
!---  Determine the number of non-zero elements on each processor,
!     number of equations on each processor
!     
      NNZ = 0
      NROW = 0
      NROW_CW = 0
      NMDMN = NFLD*ISVC
      NCMX = 0
      NCMX_CW = 0
!
!     Loop over local nodes, skipping inactive and ghost cells  ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
        NMD = (IXP(N)-1)*ISVC
        NC = 0
!
!---    Loop over number of equations  ---
!
        DO L = 1,ISVC
!
!---      Node  ---
!
          NNZ = NNZ + ISVC
          NROW = NROW + 1
          NC = NC + 1
          NMDMN = MIN( NMDMN,NMD+L)
!
!---      Bottom ---
!
          NB = ICM(1,N)
          IF( NB.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
!
!---      South ---
!
          NS = ICM(2,N)
          IF( NS.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
!
!---      West ---
!
          NW = ICM(3,N)
          IF( NW.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
!
!---      East ---
!
          NE = ICM(4,N)
          IF( NE.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
!
!---      North ---
!
          NN = ICM(5,N)
          IF( NN.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
!
!---      Top ---
!
          NT = ICM(6,N)
          IF( NT.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
          NCMX = MAX( NCMX,NC )
        ENDDO
      ENDDO
!
!---  Flow and transport equation local offsets ---
!
      IEQ_OFFSET = NMDMN - 1
!
!---  Loop over coupled wells ---
!
      DO NCW = 1,N_CW
        JM_CW(NCW) = (NFLD-NXP)*ISVC + NCW
      ENDDO
!
!---  Loop over coupled wells, adding one non-zero element
!     for the coupled-well unknown, and adding one equation ---
!
      DO NCW = 1,N_CW
        NC = 0
!
!---    Loop over coupled-well well nodes ---
!
        DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(NWN)
          IF( N.EQ.0 ) CYCLE
!
!---      Change in water mass flux into field node with respect
!         to change in coupled-well pressure  ---
!
          NNZ = NNZ + 1
!
!---      Change in CO2 mass flow with respect
!         to change in coupled-well pressure  ---
!
          NNZ = NNZ + 1
!
!---      Skip for isobrine simulations  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Change in salt mass flux into field node with respect
!           to change in coupled-well pressure  ---
!
            NNZ = NNZ + 1
          ENDIF
        ENDDO
!
!---    Coupled well equations located after field node equations
!       on last processor  ---
!
        IF( ID.EQ.(NP-1) ) THEN
          NROW_CW = NROW_CW + 1
!
!---      Change in coupled-well mass balance with respect to
!         change in coupled-well pressure  ---
!
          NC = NC + 1
          NNZ = NNZ + 1
!
!---      Loop over field nodes with coupled-well nodes ---
!
          DO NWF = ID_CW(5,NCW),ID_CW(6,NCW)
!
!---        Change in coupled-well mass balance with respect to
!           change in field node primary variables  ---
!
            DO M = 1,ISVC
              NC = NC + 1
              NNZ = NNZ + 1
            ENDDO
          ENDDO
        ENDIF
        NCMX_CW = MAX( NCMX_CW,NC )
      ENDDO
!
!---  Allocate memory for the CSR pointer (NLU), index (MLU),
!     value (DLU), and equation to index pointer (KLU) ---
!
      ALLOCATE( NLU(1:NROW+NROW_CW+1),STAT=ISTAT )
      CHMSG = 'NLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( MLU(1:NNZ),STAT=ISTAT )
      CHMSG = 'MLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DLU(1:NNZ),STAT=ISTAT )
      CHMSG = 'DLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      LJG = MAX(NROW,1)
      LJH = MAX(NCMX,1)
      ALLOCATE( KLU(1:LJG,1:LJH),STAT=ISTAT )
      CHMSG = 'KLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( KLU1_CW(1:LUK,1:LWN_CW,1:LN_CW),STAT=ISTAT )
      CHMSG = 'KLU1_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( KLU2_CW(1:NCMX_CW,1:LN_CW),STAT=ISTAT )
      CHMSG = 'KLU2_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )

      ALLOCATE( BLU(1:NUKFL(ID+1)),STAT=ISTAT )
      CHMSG = 'BLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )








!      PRINT *,'NROW = ',NROW,'NNZ = ',NNZ,'NCMX = ',NCMX,'LJG = ',LJG,
!     &  'LJH = ',LJH,'LUK = ',LUK,'LWN_CW = ',LWN_CW,'LN_CW',LN_CW,
!     &  'NCMX_CW = ',NCMX_CW,'NROW_CW = ',NROW_CW,'ID = ',ID
!
!---  Initialize memory for the CSR pointer (NLU), index (MLU),
!     value (DLU), and equation to index pointer (KLU) ---
!
      DO L = 1,NROW+NROW_CW+1
        NLU(L) = 0
      ENDDO
      DO L = 1,NNZ
        MLU(L) = 0
        DLU(L) = 0.D+0
      ENDDO
      DO M = 1,LJH
        DO L = 1,LJG
          KLU(L,M) = 0
        ENDDO
      ENDDO
      DO M = 1,LN_CW
        DO L = 1,LWN_CW
          DO K = 1,LUK
            KLU1_CW(K,L,M) = 0
          ENDDO
        ENDDO
        DO L = 1,NCMX_CW
          KLU2_CW(L,M) = 0
        ENDDO
      ENDDO
!
!---  Assign CSR pointer (NLU), index (MLU), and equation to 
!     index pointer (KLU) ---
!
      MC = 0
      NC = 0
      NLU(1) = 0
      NNZFR = 0
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
        MC = MC + 1
        NMD = (IXP(N)-1)*ISVC
!
!---    Loop over number of equations  ---
!
        DO L = 1,ISVC
          MA = 0
!
!---      Node  ---
!
          DO M = 1,ISVC
            NC = NC + 1
            MLU(NC) = NMD + M - 1
            KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
          ENDDO
          MA = MA + ISVC
!
!---      Bottom ---
!
          NB = ICM(1,N)
          IF( NB.NE.0 ) THEN
            NMDB = (IXP(NB)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDB + M - 1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
!
!---      South ---
!
          NS = ICM(2,N)
          IF( NS.NE.0 ) THEN
            NMDS = (IXP(NS)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDS + M - 1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
!
!---      West ---
!
          NW = ICM(3,N)
          IF( NW.NE.0 ) THEN
            NMDW = (IXP(NW)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDW + M - 1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
!
!---      East ---
!
          NE = ICM(4,N)
          IF( NE.NE.0 ) THEN
            NMDE = (IXP(NE)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDE + M - 1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
!
!---      North ---
!
          NN = ICM(5,N)
          IF( NN.NE.0 ) THEN
            NMDN = (IXP(NN)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDN + M - 1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
!
!---      Top ---
!
          NT = ICM(6,N)
          IF( NT.NE.0 ) THEN
            NMDT = (IXP(NT)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDT + M - 1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
!
!---      Node contains a coupled well node,
!         include coupled well equations in nodal equations  ---
!
          DO NCW = 1,N_CW
            DO NWF = ID_CW(5,NCW),ID_CW(6,NCW)
              IF( IWF_CW(NWF).EQ.N ) THEN
                NC = NC+1
                MLU(NC) = JM_CW(NCW)-1
                DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
                  IF( IWN_CW(NWN).EQ.N ) THEN
                    KLU1_CW(L,NWN,NCW) = NC
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDDO
          NLU(NMD+L-IEQ_OFFSET+1) = NC
          NNZFR = MAX( NNZFR,(NC-NLU(NMD+L-IEQ_OFFSET)) )
        ENDDO
      ENDDO
!
!---  Coupled well equations located after field node equations
!     on last processor  ---
!
      IF( ID.EQ.(NP-1) ) THEN
!
!---    Loop over coupled wells ---
!
        DO NCW = 1,N_CW
          NCCW = 0
!
!---      Change in coupled-well mass balance with respect to
!         change in coupled-well pressure  ---
!
          NC = NC + 1
          NCCW = NCCW + 1
          NMD = JM_CW(NCW)
          MLU(NC) = NMD-1
          KLU2_CW(NCCW,NCW) = NC
!
!---      Loop over field nodes with coupled-well nodes ---
!
          DO NWF = ID_CW(5,NCW),ID_CW(6,NCW)
!
!---        Change in coupled-well mass balance with respect to
!           change in field node primary variables  ---
!
            DO M = 1,ISVC
              NC = NC + 1
              NCCW = NCCW + 1
              MLU(NC) = (IXPG_CW(NWF)-1)*ISVC + M - 1
              KLU2_CW(NCCW,NCW) = NC
            ENDDO
          ENDDO
          NLU(NMD-IEQ_OFFSET+1) = NC
          NNZFR = MAX( NNZFR,(NC-NLU(NMD-IEQ_OFFSET)) )
        ENDDO
      ENDIF
!
!---  Solute or reactive species transport solution  ---
!
      IF( IEQC.NE.0 .OR. ISLC(40).NE.0 ) THEN
!
!---    Determine the number of non-zero elements on each processor,
!       number of equations on each processor,
!       
        NNZC = 0
        NROWC = 0
        NMDMNC = NFLD
        NCMXC = 0
!
!       Loop over local nodes, skipping inactive and ghost cells  ---
!
        DO N = 1,NFCGC(ID+1)
!
!---      Skip for inactive nodes or ghost cells  ---
!
          IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
          NMD = IXP(N)
          NC = 0
!
!---      Node  ---
!
          NNZC = NNZC + 1
          NROWC = NROWC + 1
          NMDMNC = MIN( NMDMNC,NMD)
          NC = NC + 1
!
!---      Bottom ---
!
          NB = ICM(1,N)
          IF( NB.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
!
!---      South ---
!
          NS = ICM(2,N)
          IF( NS.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
!
!---      West ---
!
          NW = ICM(3,N)
          IF( NW.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
!
!---      East ---
!
          NE = ICM(4,N)
          IF( NE.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
!
!---      North ---
!
          NN = ICM(5,N)
          IF( NN.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
!
!---      Top ---
!
          NT = ICM(6,N)
          IF( NT.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
          NCMXC = MAX( NCMXC,NC )
        ENDDO
!
!---    Transport equation local offsets ---
!
        IEQC_OFFSET = NMDMNC - 1
!
!---    Allocate memory for the CSR pointer (NLUC), index (MLUC),
!       and equation to index pointer (KLUC) ---
!
        ALLOCATE( NLUC(1:NROWC+1),STAT=ISTAT )
        CHMSG = 'NLUC'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        ALLOCATE( MLUC(1:NNZC),STAT=ISTAT )
        CHMSG = 'MLUC'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        LJK = MAX(NROWC,1)
        LJL = MAX(NCMXC,1)
        ALLOCATE( KLUC(1:LJK,1:LJL),STAT=ISTAT )
        CHMSG = 'KLUC'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---    Initialize memory for the CSR pointer (NLUC), index (MLUC),
!       and equation to index pointer (KLUC) ---
!
        DO L = 1,NROWC+1
          NLUC(L) = 0
        ENDDO
        DO L = 1,NNZC
          MLUC(L) = 0
        ENDDO
        DO M = 1,LJL
          DO L = 1,LJK
            KLUC(L,M) = 0
          ENDDO
        ENDDO
!
!---    Assign CSR pointer (NLUC), index (MLUC), and equation to 
!       index pointer (KLUC) ---
!
        MC = 0
        NC = 0
        NLUC(1) = 0
        NNZTR = 0
        DO N = 1,NFCGC(ID+1)
!
!---      Skip for inactive nodes or ghost cells  ---
!
          IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
          MC = MC + 1
          NMD = IXP(N)
          MA = 1
!
!---      Node  ---
!
          NC = NC + 1
          MLUC(NC) = NMD-1
          KLUC(NMD-IEQC_OFFSET,MA) = NC
          MA = MA + 1
!
!---      Bottom ---
!
          NB = ICM(1,N)
          IF( NB.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NB)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
!
!---      South ---
!
          NS = ICM(2,N)
          IF( NS.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NS)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
!
!---      West ---
!
          NW = ICM(3,N)
          IF( NW.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NW)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
!
!---      East ---
!
          NE = ICM(4,N)
          IF( NE.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NE)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
!
!---      North ---
!
          NN = ICM(5,N)
          IF( NN.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NN)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
!
!---      Top ---
!
          NT = ICM(6,N)
          IF( NT.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NT)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
          NLUC(NMD-IEQC_OFFSET+1) = NC
          NNZTR = MAX( NNZTR,(NC-NLUC(NMD-IEQC_OFFSET)) )
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBP_CW_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBP_NCW_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Configure the arrays for compressed sparse row matrix storage
!     for problems without coupled wells.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 March 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE JACOB
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
!----------------------Type Declarations-------------------------------!
!
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBP_NCW_CO2'
!
!---  Determine the number of non-zero elements on each processor,
!     number of equations on each processor,
!     
      NNZ = 0
      NROW = 0
      NMDMN = NFLD*ISVC
      NCMX = 0
!
!     Loop over local nodes, skipping inactive and ghost cells  ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
        NMD = (IXP(N)-1)*ISVC
        NC = 0
!
!---    Loop over number of equations  ---
!
        DO L = 1,ISVC
!
!---      Node  ---
!
          NNZ = NNZ + ISVC
          NROW = NROW + 1
          NC = NC + 1
          NMDMN = MIN( NMDMN,NMD+L)
!
!---      Bottom ---
!
          NB = ICM(1,N)
          IF( NB.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
!
!---      South ---
!
          NS = ICM(2,N)
          IF( NS.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
!
!---      West ---
!
          NW = ICM(3,N)
          IF( NW.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
!
!---      East ---
!
          NE = ICM(4,N)
          IF( NE.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
!
!---      North ---
!
          NN = ICM(5,N)
          IF( NN.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
!
!---      Top ---
!
          NT = ICM(6,N)
          IF( NT.NE.0 ) THEN
            NNZ = NNZ + ISVC
            NC = NC +1
          ENDIF
          NCMX = MAX( NCMX,NC )
        ENDDO
      ENDDO
!
!---  Flow and transport equation local offsets ---
!
      IEQ_OFFSET = NMDMN - 1
!
!---  Allocate memory for the CSR pointer (NLU), index (MLU),
!     value (DLU), equation to index pointer (KLU),
!     and problem/solution vector (BLU) ---
!
      ALLOCATE( NLU(1:NROW+1),STAT=ISTAT )
      CHMSG = 'NLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( MLU(1:NNZ),STAT=ISTAT )
      CHMSG = 'MLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DLU(1:NNZ),STAT=ISTAT )
      CHMSG = 'DLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )

      ALLOCATE( BLU(1:NUKFL(ID+1)),STAT=ISTAT )
      CHMSG = 'BLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )







      LJG = MAX(NROW,1)
      LJH = MAX(NCMX,1)
      ALLOCATE( KLU(1:LJG,1:LJH),STAT=ISTAT )
      CHMSG = 'KLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Initialize memory for the CSR pointer (NLU), index (MLU),
!     value (DLU), and equation to index pointer (KLU) ---
!
      DO L = 1,NROW+1
        NLU(L) = 0
      ENDDO
      DO L = 1,NNZ
        MLU(L) = 0
        DLU(L) = 0.D+0
      ENDDO
      DO M = 1,LJH
        DO L = 1,LJG
          KLU(L,M) = 0
        ENDDO
      ENDDO
!
!---  Assign CSR pointer (NLU), index (MLU), and equation to 
!     index pointer (KLU) ---
!
      MC = 0
      NC = 0
      NLU(1) = 0
      NNZFR = 0
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
        MC = MC + 1
        NMD = (IXP(N)-1)*ISVC
!
!---    Loop over number of equations  ---
!
        DO L = 1,ISVC
          MA = 0
!
!---      Node  ---
!
          DO M = 1,ISVC
            NC = NC + 1
            MLU(NC) = NMD+M-1
            KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
          ENDDO
          MA = MA + ISVC
!
!---      Bottom ---
!
          NB = ICM(1,N)
          IF( NB.NE.0 ) THEN
            NMDB = (IXP(NB)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDB+M-1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
!
!---      South ---
!
          NS = ICM(2,N)
          IF( NS.NE.0 ) THEN
            NMDS = (IXP(NS)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDS+M-1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
!
!---      West ---
!
          NW = ICM(3,N)
          IF( NW.NE.0 ) THEN
            NMDW = (IXP(NW)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDW+M-1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
!
!---      East ---
!
          NE = ICM(4,N)
          IF( NE.NE.0 ) THEN
            NMDE = (IXP(NE)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDE+M-1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
!
!---      North ---
!
          NN = ICM(5,N)
          IF( NN.NE.0 ) THEN
            NMDN = (IXP(NN)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDN+M-1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
!
!---      Top ---
!
          NT = ICM(6,N)
          IF( NT.NE.0 ) THEN
            NMDT = (IXP(NT)-1)*ISVC
            DO M = 1,ISVC
              NC = NC + 1
              MLU(NC) = NMDT+M-1
              KLU(NMD+L-IEQ_OFFSET,M+MA) = NC
            ENDDO
            MA = MA + ISVC
          ENDIF
          NLU(NMD+L-IEQ_OFFSET+1) = NC
          NNZFR = MAX( NNZFR,(NC-NLU(NMD+L-IEQ_OFFSET)) )
        ENDDO
      ENDDO
!
!---  Solute or reactive species transport solution  ---
!
      IF( IEQC.NE.0 .OR. ISLC(40).NE.0 ) THEN
!
!---    Determine the number of non-zero elements on each processor,
!       number of equations on each processor,
!       
        NNZC = 0
        NROWC = 0
        NMDMNC = NFLD
        NCMXC = 0
!
!       Loop over local nodes, skipping inactive and ghost cells  ---
!
        DO N = 1,NFCGC(ID+1)
!
!---      Skip for inactive nodes or ghost cells  ---
!
          IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
          NMD = IXP(N)
          NC = 0
!
!---      Node  ---
!
          NNZC = NNZC + 1
          NROWC = NROWC + 1
          NMDMNC = MIN( NMDMNC,NMD)
          NC = NC + 1
!
!---      Bottom ---
!
          NB = ICM(1,N)
          IF( NB.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
!
!---      South ---
!
          NS = ICM(2,N)
          IF( NS.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
!
!---      West ---
!
          NW = ICM(3,N)
          IF( NW.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
!
!---      East ---
!
          NE = ICM(4,N)
          IF( NE.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
!
!---      North ---
!
          NN = ICM(5,N)
          IF( NN.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
!
!---      Top ---
!
          NT = ICM(6,N)
          IF( NT.NE.0 ) THEN
            NNZC = NNZC + 1
            NC = NC + 1
          ENDIF
          NCMXC = MAX( NCMXC,NC )
        ENDDO
!
!---    Transport equation local offsets ---
!
        IEQC_OFFSET = NMDMNC - 1
!
!---    Allocate memory for the CSR pointer (NLUC), index (MLUC),
!       and equation to index pointer (KLUC) ---
!
        ALLOCATE( NLUC(1:NROWC+1),STAT=ISTAT )
        CHMSG = 'NLUC'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        ALLOCATE( MLUC(1:NNZC),STAT=ISTAT )
        CHMSG = 'MLUC'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        LJK = MAX(NROWC,1)
        LJL = MAX(NCMXC,1)
        ALLOCATE( KLUC(1:LJK,1:LJL),STAT=ISTAT )
        CHMSG = 'KLUC'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---    Initialize memory for the CSR pointer (NLUC), index (MLUC),
!       and equation to index pointer (KLUC) ---
!
        DO L = 1,NROWC+1
          NLUC(L) = 0
        ENDDO
        DO L = 1,NNZC
          MLUC(L) = 0
        ENDDO
        DO M = 1,LJL
          DO L = 1,LJK
            KLUC(L,M) = 0
          ENDDO
        ENDDO
!
!---    Assign CSR pointer (NLUC), index (MLUC), and equation to 
!       index pointer (KLUC) ---
!
        MC = 0
        NC = 0
        NLUC(1) = 0
        NNZTR = 0
        DO N = 1,NFCGC(ID+1)
!
!---      Skip for inactive nodes or ghost cells  ---
!
          IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
          MC = MC + 1
          NMD = IXP(N)
          MA = 1
!
!---      Node  ---
!
          NC = NC + 1
          MLUC(NC) = NMD-1
          KLUC(NMD-IEQC_OFFSET,MA) = NC
          MA = MA + 1
!
!---      Bottom ---
!
          NB = ICM(1,N)
          IF( NB.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NB)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
!
!---      South ---
!
          NS = ICM(2,N)
          IF( NS.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NS)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
!
!---      West ---
!
          NW = ICM(3,N)
          IF( NW.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NW)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
!
!---      East ---
!
          NE = ICM(4,N)
          IF( NE.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NE)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
!
!---      North ---
!
          NN = ICM(5,N)
          IF( NN.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NN)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
!
!---      Top ---
!
          NT = ICM(6,N)
          IF( NT.NE.0 ) THEN
            NC = NC + 1
            MLUC(NC) = IXP(NT)-1
            KLUC(NMD-IEQC_OFFSET,MA) = NC
            MA = MA + 1
          ENDIF
          NLUC(NMD-IEQC_OFFSET+1) = NC
          NNZTR = MAX( NNZTR,(NC-NLUC(NMD-IEQC_OFFSET)) )
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBP_NCW_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBL_CO2( RSS,RSP,RSA,N,MEQ )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Load the Jacobian matrix.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 4 November 2021.
!

!----------------------LIS Modules-------------------------------------!
!
      USE LIS_STOMP
!







!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
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
      REAL*8 RSP(LUK),RSA(LUK,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBL_CO2'

!
!---  Lis solver  ---
!
!
!---  Node  ---
!
      NMD = (IXP(N)-1)*ISVC
      IROW =  NMD + MEQ
      MA = 0
      DO M = 1,ISVC
        MCOL = KLU(IROW-IEQ_OFFSET,M+MA)
        DLU(MCOL) = DLU(MCOL) + (RSP(M)-RSS)/DNR(M,N)
      ENDDO
      BUFFER = -RSS
      CALL lis_vector_set_value_f( 1,IROW,BUFFER,
     &    F_RHS_VEC,IERR )
      RSDL(MEQ,N) = RSDL(MEQ,N) - RSS
      MA = MA + ISVC
!
!---  Bottom ---
!
      NB = ICM(1,N)
      IF( NB.NE.0 ) THEN
        DO M = 1,ISVC
          MCOL = KLU(IROW-IEQ_OFFSET,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSA(M,1)-RSS)/DNR(M,NB)
        ENDDO
        MA = MA + ISVC
      ENDIF
!
!---  South ---
!
      NS = ICM(2,N)
      IF( NS.NE.0 ) THEN
        DO M = 1,ISVC
          MCOL = KLU(IROW-IEQ_OFFSET,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSA(M,2)-RSS)/DNR(M,NS)
        ENDDO
        MA = MA + ISVC
      ENDIF
!
!---  West ---
!
      NW = ICM(3,N)
      IF( NW.NE.0 ) THEN
        DO M = 1,ISVC
          MCOL = KLU(IROW-IEQ_OFFSET,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSA(M,3)-RSS)/DNR(M,NW)
        ENDDO
        MA = MA + ISVC
      ENDIF
!
!---  East ---
!
      NE = ICM(4,N)
      IF( NE.NE.0 ) THEN
        DO M = 1,ISVC
          MCOL = KLU(IROW-IEQ_OFFSET,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSA(M,4)-RSS)/DNR(M,NE)
        ENDDO
        MA = MA + ISVC
      ENDIF
!
!---  North ---
!
      NN = ICM(5,N)
      IF( NN.NE.0 ) THEN
        DO M = 1,ISVC
          MCOL = KLU(IROW-IEQ_OFFSET,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSA(M,5)-RSS)/DNR(M,NN)
        ENDDO
        MA = MA + ISVC
      ENDIF
!
!---  Top ---
!
      NT = ICM(6,N)
      IF( NT.NE.0 ) THEN
        DO M = 1,ISVC
          MCOL = KLU(IROW-IEQ_OFFSET,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSA(M,6)-RSS)/DNR(M,NT)
        ENDDO
        MA = MA + ISVC
      ENDIF

!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBL_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBA_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Load the Jacobian matrix for the CO2 equation with
!     aqueous and gas contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 4 November 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOURC
      USE SOLTN
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STAX(LUK+1),RAP(LUK),RAA(LUK,6),FA(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBA_CO2'
      IEQAX = IEQA
!
!---  Loop over local nodes  ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---    First-order, forward-difference, time differential  ---
!
        STA1 = PORD(1,N)*(XLA(1,N)*RHOL(1,N)*SL(1,N) +
     &      XGA(1,N)*RHOG(1,N)*SG(1,N))
        DO M = 1,ISVC+1
          MP = M + 1
          STA0 = PORD(MP,N)*(XLA(MP,N)*RHOL(MP,N)*SL(MP,N) +
     &      XGA(MP,N)*RHOG(MP,N)*SG(MP,N))
          STAX(M) = (STA0-STA1)*DTI*VOL(N)
        ENDDO
!
!---  Initialize surface fluxes  ---
!
        DO MD = 1,6
          DO M = 1,ISVF
            FA(M,MD) = 0.D+0
          ENDDO
        ENDDO
!
!---  Compute surface fluxes  ---
!
!---    Bottom ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          DO M = 1,ISVF
            FA(M,1) = -AFZ(1,N)*(WLA(M,1,N)+WGA(M,1,N))
          ENDDO
        ENDIF
!
!---    South ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DO M = 1,ISVF
            FA(M,2) = -AFY(1,N)*(VLA(M,1,N)+VGA(M,1,N))
          ENDDO
        ENDIF
!
!---    West ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          DO M = 1,ISVF
            FA(M,3) = -AFX(1,N)*(ULA(M,1,N)+UGA(M,1,N))
          ENDDO
        ENDIF
!
!---    East ---
!
        NE = ICM(4,N)
        IF( NE.NE.0 ) THEN
          DO M = 1,ISVF
            MF = MFLX(M)
            FA(M,4) = AFX(2,N)*(ULA(MF,2,N)+UGA(MF,2,N))
          ENDDO
        ENDIF
!
!---    North ---
!
        NN = ICM(5,N)
        IF( NN.NE.0 ) THEN
          DO M = 1,ISVF
            MF = MFLX(M)
            FA(M,5) = AFY(2,N)*(VLA(MF,2,N)+VGA(MF,2,N))
          ENDDO
        ENDIF
!
!---    Top ---
!
        NT = ICM(6,N)
        IF( NT.NE.0 ) THEN
          DO M = 1,ISVF
            MF = MFLX(M)
            FA(M,6) = AFZ(2,N)*(WLA(MF,2,N)+WGA(MF,2,N))
          ENDDO
        ENDIF
!
!---    Compute CO2 equation residuals  ---
!
        RAS = STAX(1) - SRCA(2,N)
        DO MD = 1,6
          RAS = RAS + FA(1,MD)
        ENDDO
        DO M = 1,ISVC
          RAP(M) = STAX(M+1) - SRCA(M+2,N)
          MM = 2*M
          DO MD = 1,6
            RAP(M) = RAP(M) + FA(MM,MD)
          ENDDO
        ENDDO
        DO M = 1,ISVC
          MM = 2*M + 1
          DO MD = 1,6
            RAA(M,MD) = RAS - FA(1,MD) + FA(MM,MD)
          ENDDO
        ENDDO
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBL_CO2( RAS,RAP,RAA,N,IEQAX )
!
!---  Continue to Next Node  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBA_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Load the Jacobian matrix for the salt equation with
!     aqueous contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 4 November 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE GRID
      USE FLUX
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STSX(LUK+1),RSP(LUK),RSA(LUK,6),FS(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBS_CO2'
      IEQSX = IEQS
!
!---  Loop over local nodes  ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---    First-order, forward-difference, time differential  ---
!
        DO M = 1,ISVC+1
          MP = M + 1
          STSX(M) = (TMS(MP,N)-TMS(1,N))*VOL(N)*DTI
        ENDDO
!
!---    Initialize surface fluxes  ---
!
        DO MD = 1,6
          DO M = 1,ISVF
            FS(M,MD) = 0.D+0
          ENDDO
        ENDDO
!
!---    Bottom surface flux component  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          DO M = 1,ISVF
            FS(M,1) = -AFZ(1,N)*WS(M,1,N)
          ENDDO
        ENDIF
!
!---    South surface flux component ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DO M = 1,ISVF
            FS(M,2) = -AFY(1,N)*VS(M,1,N)
          ENDDO
        ENDIF
!
!---    West surface flux component ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          DO M = 1,ISVF
            FS(M,3) = -AFX(1,N)*US(M,1,N)
          ENDDO
        ENDIF
!
!---    East surface flux component  ---
!
        NE = ICM(4,N)
        IF( NE.NE.0 ) THEN
          DO M = 1,ISVF
            MF = MFLX(M)
            FS(M,4) = AFX(2,N)*US(MF,2,N)
          ENDDO
        ENDIF
!
!---    North surface flux component  ---
!
        NN = ICM(5,N)
        IF( NN.NE.0 ) THEN
          DO M = 1,ISVF
            MF = MFLX(M)
            FS(M,5) = AFY(2,N)*VS(MF,2,N)
          ENDDO
        ENDIF
!
!---    Top surface flux component  ---
!
        NT = ICM(6,N)
        IF( NT.NE.0 ) THEN
          DO M = 1,ISVF
            MF = MFLX(M)
            FS(M,6) = AFZ(2,N)*WS(MF,2,N)
          ENDDO
        ENDIF
!
!---    Salt equation residuals  ---
!
        RSS = STSX(1) - SRCS(2,N)
        DO MD = 1,6
          RSS = RSS + FS(1,MD)
        ENDDO
        DO M = 1,ISVC
          RSP(M) = STSX(M+1) - SRCS(M+2,N)
          MM = 2*M
          DO MD = 1,6
            RSP(M) = RSP(M) + FS(MM,MD)
          ENDDO
        ENDDO
        DO M = 1,ISVC
          MM = 2*M + 1
          DO MD = 1,6
            RSA(M,MD) = RSS - FS(1,MD) + FS(MM,MD)
          ENDDO
        ENDDO
!
!---    Jacobian matrix loader  --
!
        CALL JCBL_CO2( RSS,RSP,RSA,N,IEQSX )
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Load the Jacobian matrix for the water equation with
!     aqueous and gas contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 4 November 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE GRID
      USE FLUX
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STWX(LUK+1),RWP(LUK),RWA(LUK,6),FW(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBW_CO2'
      IEQWX = IEQW
!
!---  Loop over local nodes  ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---    First-order, forward-difference, time differential  ---
!
        STW1 = PORD(1,N)*(XLW(1,N)*RHOL(1,N)*SL(1,N) +
     &      XGW(1,N)*RHOG(1,N)*SG(1,N))
        DO M = 1,ISVC+1
          MP = M + 1
          STW0 = PORD(MP,N)*(XLW(MP,N)*RHOL(MP,N)*SL(MP,N) +
     &      XGW(MP,N)*RHOG(MP,N)*SG(MP,N))
          STWX(M) = (STW0-STW1)*DTI*VOL(N)
        ENDDO
!
!---  Initialize surface fluxes  ---
!
        DO MD = 1,6
          DO M = 1,ISVF
            FW(M,MD) = 0.D+0
          ENDDO
        ENDDO
!
!---    Compute surface fluxes  ---
!
!---    Bottom ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          DO M = 1,ISVF
            FW(M,1) = -AFZ(1,N)*(WLW(M,1,N)+WGW(M,1,N))
          ENDDO
        ENDIF
!
!---    South ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DO M = 1,ISVF
            FW(M,2) = -AFY(1,N)*(VLW(M,1,N)+VGW(M,1,N))
          ENDDO
        ENDIF
!
!---    West ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          DO M = 1,ISVF
            FW(M,3) = -AFX(1,N)*(ULW(M,1,N)+UGW(M,1,N))
          ENDDO
        ENDIF
!
!---    East ---
!
        NE = ICM(4,N)
        IF( NE.NE.0 ) THEN
          DO M = 1,ISVF
            MF = MFLX(M)
            FW(M,4) = AFX(2,N)*(ULW(MF,2,N)+UGW(MF,2,N))
          ENDDO
        ENDIF
!
!---    North ---
!
        NN = ICM(5,N)
        IF( NN.NE.0 ) THEN
          DO M = 1,ISVF
            MF = MFLX(M)
            FW(M,5) = AFY(2,N)*(VLW(MF,2,N)+VGW(MF,2,N))
          ENDDO
        ENDIF
!
!---    Top ---
!
        NT = ICM(6,N)
        IF( NT.NE.0 ) THEN
          DO M = 1,ISVF
            MF = MFLX(M)
            FW(M,6) = AFZ(2,N)*(WLW(MF,2,N)+WGW(MF,2,N))
          ENDDO
        ENDIF
!
!---    Compute water equation residuals  ---
!
        RWS = STWX(1) - SRCW(2,N)
        DO MD = 1,6
          RWS = RWS + FW(1,MD)
        ENDDO
        DO M = 1,ISVC
          RWP(M) = STWX(M+1) - SRCW(M+2,N)
          MM = 2*M
          DO MD = 1,6
            RWP(M) = RWP(M) + FW(MM,MD)
          ENDDO
        ENDDO
        DO M = 1,ISVC
          MM = 2*M + 1
          DO MD = 1,6
            RWA(M,MD) = RWS - FW(1,MD) + FW(MM,MD)
          ENDDO
        ENDDO
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBL_CO2( RWS,RWP,RWA,N,IEQWX )
!
!---  Continue to Next Node  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBZ_CO2
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Zero the Jacobian matrix.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 March 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE JACOB
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
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBZ_CO2'
!
!---  Zero problem vector array for flow and transport  ---
!

      DO N = 1,NUKFL(ID+1)
        BLU(N) = 0.D+0
      ENDDO








!
!---  Zero Jacobian matrix array for flow and transport  ---
!
      DO N = 1,NNZ
        DLU(N) = 0.D+0
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBZ_CO2 group
!
      RETURN
      END




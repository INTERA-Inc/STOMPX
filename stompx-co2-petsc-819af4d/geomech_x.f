!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCJ_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Modify the geomechanics Jacobian matrix and problem vector
!     for geomechanics boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 November 2016 (Magnus Carlsen
!     birthdate, Norwegian chess grandmaster).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GEO_MECH
      USE CONST







!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP

!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)





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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 BCX(LBCV_GM)
      REAL*8 XEX(5),YEX(5),ZEX(5)
      INTEGER IBX(5),ISX(5),IWX(5),IEX(5),INX(5),ITX(5)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX




      integer(kind=selected_int_kind(5)) :: IERR

!
!----------------------Data Statements---------------------------------!
!
      DATA IBX / 1,3,4,2,1 /
      DATA ISX / 1,2,6,5,1 /
      DATA IWX / 1,5,7,3,1 /
      DATA IEX / 2,4,8,6,2 /
      DATA INX / 3,7,8,4,3 /
      DATA ITX / 5,6,8,7,5 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCJ_GM'
!
!---  For restart simulations, eliminate the reference 
!     boundary condition  ---
!
      IF( IEO.EQ.2 ) THEN
        DO NB = 1,NBC_GM(ID+1)
          MB = IBCIN_GM(NB)
          DO M = 1,IBCM_GM(NB)
            IF( M.LE.IBCR_GM(MB) ) CYCLE
            IF( M.GT.IBCR_GM(MB) ) THEN
              DO N = 1,LBCV_GM
                BC_GM(N,M-1,MB) = BC_GM(N,M,MB)
              ENDDO
            ENDIF
          ENDDO
          IBCM_GM(NB) = IBCM_GM(NB)-1
          IBCR_GM(MB) = 0
        ENDDO
      ENDIF
!
!---  Loop over the geomechanical boundary conditions  ---
!
      DO NB = 1,NBC_GM(ID+1)
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN_GM(NB)
        IF( IBCC_GM(NB).EQ.1 ) TMZ = MOD( TM,BC_GM(1,IBCM_GM(NB),MB) )
        IF( TMZ.LE.BC_GM(1,1,MB) ) CYCLE
!
!---    Reference boundary condition  ---
!
        IF( IREF_GM.EQ.1 .AND. IBCR_GM(MB).GT.0 ) THEN
!
!---      Assign local boundary condition variables  ---
!
          DO N = 1,LBCV_GM
            BCX(N) = BC_GM(N,IBCR_GM(MB),MB)
          ENDDO
!
!---    Non-reference boundary condition  ---
!
        ELSE
!
!---      Single geomechanical boundary condition time, (one time and 
!         no reference value, or two times and reference value )  ---
!
          IF( (IBCM_GM(NB).EQ.1 .AND. IBCR_GM(MB).EQ.0)
     &      .OR. (IBCM_GM(NB).EQ.2 .AND. IBCR_GM(MB).NE.0) ) THEN
            DO M = 1,IBCM_GM(NB)
              IF( M.EQ.IBCR_GM(MB) ) CYCLE
              MA = M
            ENDDO
!
!---        Assign local boundary condition variables  ---
!
            DO N = 1,LBCV_GM
              BCX(N) = BC_GM(N,MA,MB)
            ENDDO
!
!---      Mulitple geomechanical boundary condition times  ---
!
          ELSEIF( IREF_GM.EQ.0 ) THEN
            IFIND = 0
            DO M = 1,IBCM_GM(NB)
              IF( M.EQ.IBCR_GM(MB) ) CYCLE
              IF( TMZ.LE.BC_GM(1,M,MB) ) THEN
                TDBC = (BC_GM(1,M,MB)-BC_GM(1,M-1,MB))
                DTBC = MIN( BC_GM(1,M,MB)-TMZ,DT )
                IF( NSTEP-NRST.EQ.0 .AND. IBCC_GM(NB).EQ.0 ) THEN
                  TMZX = TM
                ELSE
                  TMZX = TMZ
                ENDIF
                TFBC = (TMZX-BC_GM(1,M-1,MB))/TDBC
!
!---            Assign local boundary condition variables  ---
!
                DO N = 1,LBCV_GM
                  BCX(N) = BC_GM(N,M-1,MB) + 
     &              TFBC*(BC_GM(N,M,MB)-BC_GM(N,M-1,MB))
                ENDDO
                IFIND = 1
                EXIT
              ENDIF
            ENDDO
!
!---        Time not within geomechanical boundary condition limits
!           proceed to next geomechanical boundary condition  ---
!
            IF( IFIND.EQ.0 ) CYCLE
          ENDIF
        ENDIF
!
!---    Loop over the x-, y-, and z-directions  ---
!
        DO M0 = 1,3
!
!---      Traction boundary type, compute surface normal and projection
!         of surface normal on the global x-, y-, or z-directions, 
!         reduce boundary surface area by projection and multiple
!         stress by projection surface area, force imposed on each of
!         the FE nodes in element surface is 1/4 of the stress times
!         the projected surface area  ---
!
          IF( IBCT_GM(M0,NB).EQ.1 .OR. IBCT_GM(M0,NB).EQ.2 ) THEN
!
!---        Local field node number  ---
!
            N = IBCN_GM(NB)
!
!---        Skip ghost cells  ---
!
            IF( IGHC(N).EQ.1 ) CYCLE
!
!---        Bottom surface  ---
!
            IF( IBCD_GM(NB).EQ.-3 ) THEN
              DO I = 1,5
                XEX(I) = XE(IBX(I),N)
                YEX(I) = YE(IBX(I),N)
                ZEX(I) = ZE(IBX(I),N)
              ENDDO
              AX = AFZ(1,N)
!
!---        South surface  ---
!
            ELSEIF( IBCD_GM(NB).EQ.-2 ) THEN
              DO I = 1,5
                XEX(I) = XE(ISX(I),N)
                YEX(I) = YE(ISX(I),N)
                ZEX(I) = ZE(ISX(I),N)
              ENDDO
              AX = AFY(1,N)
!
!---        West surface  ---
!
            ELSEIF( IBCD_GM(NB).EQ.-1 ) THEN
              DO I = 1,5
                XEX(I) = XE(IWX(I),N)
                YEX(I) = YE(IWX(I),N)
                ZEX(I) = ZE(IWX(I),N)
              ENDDO
              AX = AFX(1,N)
!
!---        East surface  ---
!
            ELSEIF( IBCD_GM(NB).EQ.1 ) THEN
              DO I = 1,5
                XEX(I) = XE(IEX(I),N)
                YEX(I) = YE(IEX(I),N)
                ZEX(I) = ZE(IEX(I),N)
              ENDDO
              AX = AFX(2,N)
!
!---        North surface  ---
!
            ELSEIF( IBCD_GM(NB).EQ.2 ) THEN
              DO I = 1,5
                XEX(I) = XE(INX(I),N)
                YEX(I) = YE(INX(I),N)
                ZEX(I) = ZE(INX(I),N)
              ENDDO
              AX = AFY(2,N)
!
!---        Top surface  ---
!
            ELSEIF( IBCD_GM(NB).EQ.3 ) THEN
              DO I = 1,5
                XEX(I) = XE(ITX(I),N)
                YEX(I) = YE(ITX(I),N)
                ZEX(I) = ZE(ITX(I),N)
              ENDDO
              AX = AFZ(2,N)
            ENDIF
!
!---        Surface normal x-, y-, and z-components  ---
!
            XNX = 0.D+0
            YNX = 0.D+0
            ZNX = 0.D+0
            DO I = 1,4
              XNX = XNX + (YEX(I)-YEX(I+1))*(ZEX(I)+ZEX(I+1))
              YNX = YNX + (ZEX(I)-ZEX(I+1))*(XEX(I)+XEX(I+1))
              ZNX = ZNX + (XEX(I)-XEX(I+1))*(YEX(I)+YEX(I+1))
            ENDDO
!
!---        Surface x-, y-, and z-projections  ---
!
            XPX = XNX/SQRT((XNX**2) + (YNX**2) + (ZNX**2))
            YPX = YNX/SQRT((XNX**2) + (YNX**2) + (ZNX**2))
            ZPX = ZNX/SQRT((XNX**2) + (YNX**2) + (ZNX**2))
!
!---        FE nodal force, Pa  ---
!
            IF( IBCT_GM(M0,NB).EQ.1 ) THEN
              VARX = BCX((M0-1)*3+2)
!
!---        FE nodal force, Pa w/ z-gradient  ---
!
            ELSEIF( IBCT_GM(M0,NB).EQ.2 ) THEN
              ZERX = 0.D+0
              DO I = 1,4
                ZERX = ZERX + ZEX(I)
              ENDDO
              ZERX = 2.5D-1*ZERX
              VARX = BCX((M0-1)*3+2) + (ZERX-BCX((M0-1)*3+3))*
     &          BCX((M0-1)*3+4)    
            ENDIF
            IF( M0.EQ.1 ) THEN
              FX = AX*XPX*VARX
            ELSEIF( M0.EQ.2 ) THEN
              FX = AX*YPX*VARX
            ELSEIF( M0.EQ.3 ) THEN
              FX = AX*ZPX*VARX
            ENDIF
!
!---        Loop over the local FE nodes on the element surface  ---
!
            DO I = 1,4
              IF( IBCD_GM(NB).EQ.-3 ) THEN
                NFEN = ND_GM(IBX(I),N)
              ELSEIF( IBCD_GM(NB).EQ.-2 ) THEN
                NFEN = ND_GM(ISX(I),N)
              ELSEIF( IBCD_GM(NB).EQ.-1 ) THEN
                NFEN = ND_GM(IWX(I),N)
              ELSEIF( IBCD_GM(NB).EQ.1 ) THEN
                NFEN = ND_GM(IEX(I),N)
              ELSEIF( IBCD_GM(NB).EQ.2 ) THEN
                NFEN = ND_GM(INX(I),N)
              ELSEIF( IBCD_GM(NB).EQ.3 ) THEN
                NFEN = ND_GM(ITX(I),N)
              ENDIF
!
!---          Equation number  ---
!
              NMD = (IM_GM(NFEN)-1)*3 + M0
              BUFFER = -2.5D-1*FX





              IROW_P = NMD - 1
              CALL VecSetValues( G_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,
     &          IERR )

            ENDDO
!
!---      Displacement boundary type, assign boundary condition
!         value to problem vector, set diagonal term to 1.0, and set
!         off diagonal terms to 0.0  ---
!
          ELSEIF( IBCT_GM(M0,NB).EQ.3 ) THEN
!
!---        Local FE node number  ---
!
            NFEN = IBCN_GM(NB)
            NMD = (IM_GM(NFEN)-1)*3 + M0
!
!---        Skip ghost finite-element nodes  ---
!
            IF( IGHN(NFEN).EQ.1 ) CYCLE
!
!---        Equation number  ---
!
            NMD = (IM_GM(NFEN)-1)*3 + M0
            NMDX = NMD-IEQ_OFFSET_GM
            ISKPX = NLU_GM(NMDX+1) -  NLU_GM(NMDX)
!
!---        Loop over the elements connected to the FE node  ---
!
            DO M1 = 1,8
              NX = NE_GM(M1,NFEN)
              IF( NX.EQ.0 ) CYCLE
!
!---          Loop over the FE nodes of the element  ---
!
              DO M2 = 1,8
                MFEN = ND_GM(M2,NX)
!
!---            Loop over the displacement equations  ---
!
                DO M3 = 1,3
!
!---              Global column number  ---
!
                  MCOL = (IM_GM(MFEN)-1)*3 + M3
!
!---              Sparse row data structure pointer  ---
!
                  ISTCX = NK_GM(M2,M1)
                  MROW = (M0-1)*ISKPX + (KLU_GM(ISTCX,NFEN)-1) + M3
!
!---              Diagonal term  ---
!
                  IF( NMD.EQ.MCOL ) THEN
                    DIAGX = DLU(MROW)
!                    DLU(MROW) = 1.D+0
!
!---              Off diagonal terms  ---
!
                  ELSE
                    DLU(MROW) = 0.D+0
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            BUFFER = BCX((M0-1)*3+2)*DIAGX





            IROW_P = NMD - 1
            CALL VecSetValues( G_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,
     &        IERR )

!
!---      Force boundary type, assign boundary condition
!         value to problem vector  ---
!
          ELSEIF( IBCT_GM(M0,NB).EQ.5 ) THEN
!
!---        FE node  ---
!
            NFEN = IBCN_GM(NB)
!
!---        Equation number  ---
!
            NMD = (IM_GM(NFEN)-1)*3 + M0
            BUFFER = BCX((M0-1)*3+2)





            IROW_P = NMD - 1
            CALL VecSetValues( G_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,
     &        IERR )

          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCJ_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DISP_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Update displacements at FE ghost nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 31 January 2017 (3M begins marketing
!     Scotch Tape in 1930).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GEO_MECH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DISP_GM'
!
!---  Loop over the finite-element nodes ---
!
      M0 = 0
      DO NFEN = 1,NFNGN(ID+1)






        U_GM(2,NFEN) = BLU(M0+1)
        V_GM(2,NFEN) = BLU(M0+2)
        W_GM(2,NFEN) = BLU(M0+3)
        M0 = M0 + 3
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DISP_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DISP_GN_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Update the displacements on finite-element ghost nodes
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 November 2022.
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
      USE GEO_MECH
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
      SUB_LOG(ISUB_LOG) = '/DISP_GN_GM'
      NPVX = 3
!
!---  Load sending buffer for bottom ghost FE nodes for processors
!     with bottom ghost FE nodes, and send buffer to receiving
!     processors  ---
!
      IF( NCGN(1,ID+1).GT.0 ) THEN
        MCS = 0
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGN(1,ID+1) = ',NCGN(1,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGN(1,ID+1)
          SBFB(NCS+1) = U_GM(2,NLSGN(M+MCS))
          SBFB(NCS+2) = V_GM(2,NLSGN(M+MCS))
          SBFB(NCS+3) = W_GM(2,NLSGN(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGN(1,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGN(1,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Bottom Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFB,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post Bottom Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending bottom ghost FE nodes,
!     and unload receiving buffer  ---
!
      IF( NCGN(6,ID+1).GT.0 ) THEN
        NRCVX = NCGN(6,ID+1)*NPVX
        IDSNDX = NPGN(6,ID+1) - 1
        IDRCVX = NPGN(1,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Bottom Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFB,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post Bottom Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,5
          MCR = MCR + NCGN(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGN(6,ID+1) = ',NCGN(6,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGN(6,ID+1)
          U_GM(2,NLRGN(M+MCR)) = RBFB(NCR+1)
          V_GM(2,NLRGN(M+MCR)) = RBFB(NCR+2)
          W_GM(2,NLRGN(M+MCR)) = RBFB(NCR+3)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for south ghost FE nodes for processors
!     with south ghost FE nodes, and send buffer to receiving
!     processors  ---
!
      IF( NCGN(2,ID+1).GT.0 ) THEN
        MCS = NCGN(1,ID+1)
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGN(2,ID+1) = ',NCGN(2,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGN(2,ID+1)
          SBFS(NCS+1) = U_GM(2,NLSGN(M+MCS))
          SBFS(NCS+2) = V_GM(2,NLSGN(M+MCS))
          SBFS(NCS+3) = W_GM(2,NLSGN(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGN(2,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGN(2,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'South Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFS,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post South Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending south ghost FE nodes,
!     and unload receiving buffer  ---
!
      IF( NCGN(5,ID+1).GT.0 ) THEN
        NRCVX = NCGN(5,ID+1)*NPVX
        IDSNDX = NPGN(5,ID+1) - 1
        IDRCVX = NPGN(2,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'South Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFS,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post South Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,4
          MCR = MCR + NCGN(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGN(5,ID+1) = ',NCGN(5,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGN(5,ID+1)
          U_GM(2,NLRGN(M+MCR)) = RBFS(NCR+1)
          V_GM(2,NLRGN(M+MCR)) = RBFS(NCR+2)
          W_GM(2,NLRGN(M+MCR)) = RBFS(NCR+3)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for west ghost FE nodes for processors
!     with west ghost FE nodes, and send buffer to receiving
!     processors  ---
!
      IF( NCGN(3,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,2
          MCS = MCS + NCGN(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGN(3,ID+1) = ',NCGN(3,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGN(3,ID+1)
!          IF( ID.EQ.2 )
!     &      PRINT *,'S2: ND(',NLSGN(M+MCS),') = ',ND(NLSGN(M+MCS))
          SBFW(NCS+1) = U_GM(2,NLSGN(M+MCS))
          SBFW(NCS+2) = V_GM(2,NLSGN(M+MCS))
          SBFW(NCS+3) = W_GM(2,NLSGN(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGN(3,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGN(3,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'West Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFW,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post West Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending west ghost FE nodes,
!     and unload receiving buffer  ---
!
      IF( NCGN(4,ID+1).GT.0 ) THEN
        NRCVX = NCGN(4,ID+1)*NPVX
        IDSNDX = NPGN(4,ID+1) - 1
        IDRCVX = NPGN(3,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'West Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFW,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post West Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,3
          MCR = MCR + NCGN(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGN(4,ID+1) = ',NCGN(4,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGN(4,ID+1)
!          IF( ID.EQ.1 )
!     &      PRINT *,'R1: ND(',NLSGN(M+MCR),') = ',ND(NLSGN(M+MCR))
          U_GM(2,NLRGN(M+MCR)) = RBFW(NCR+1)
          V_GM(2,NLRGN(M+MCR)) = RBFW(NCR+2)
          W_GM(2,NLRGN(M+MCR)) = RBFW(NCR+3)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for east ghost FE nodes for processors
!     with east ghost FE nodes, and send buffer to receiving
!     processors  ---
!
      IF( NCGN(4,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,3
          MCS = MCS + NCGN(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGN(4,ID+1) = ',NCGN(4,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGN(4,ID+1)
          SBFE(NCS+1) = U_GM(2,NLSGN(M+MCS))
          SBFE(NCS+2) = V_GM(2,NLSGN(M+MCS))
          SBFE(NCS+3) = W_GM(2,NLSGN(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGN(4,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGN(4,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'East Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFE,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post East Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending east ghost FE nodes,
!     and unload receiving buffer  ---
!
      IF( NCGN(3,ID+1).GT.0 ) THEN
        NRCVX = NCGN(3,ID+1)*NPVX
        IDSNDX = NPGN(3,ID+1) - 1
        IDRCVX = NPGN(4,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'East Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFE,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post East Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,2
          MCR = MCR + NCGN(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGN(3,ID+1) = ',NCGN(3,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGN(3,ID+1)
          U_GM(2,NLRGN(M+MCR)) = RBFE(NCR+1)
          V_GM(2,NLRGN(M+MCR)) = RBFE(NCR+2)
          W_GM(2,NLRGN(M+MCR)) = RBFE(NCR+3)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for north ghost FE nodes for processors
!     with north ghost FE nodes, and send buffer to receiving
!     processors  ---
!
      IF( NCGN(5,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,4
          MCS = MCS + NCGN(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGN(5,ID+1) = ',NCGN(5,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGN(5,ID+1)
          SBFN(NCS+1) = U_GM(2,NLSGN(M+MCS))
          SBFN(NCS+2) = V_GM(2,NLSGN(M+MCS))
          SBFN(NCS+3) = W_GM(2,NLSGN(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGN(5,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGN(5,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'North Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFN,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post North Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending north ghost FE nodes,
!     and unload receiving buffer  ---
!
      IF( NCGN(2,ID+1).GT.0 ) THEN
        NRCVX = NCGN(2,ID+1)*NPVX
        IDSNDX = NPGN(2,ID+1) - 1
        IDRCVX = NPGN(5,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'North Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFN,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post North Receive: IERR = ',IERR,' ID = ',ID
        MCR = NCGN(1,ID+1)
!        PRINT *,' MCR = ',MCR,' NCGN(2,ID+1) = ',NCGN(2,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGN(2,ID+1)
          U_GM(2,NLRGN(M+MCR)) = RBFN(NCR+1)
          V_GM(2,NLRGN(M+MCR)) = RBFN(NCR+2)
          W_GM(2,NLRGN(M+MCR)) = RBFN(NCR+3)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for top ghost FE nodes for processors
!     with top ghost FE nodes, and send buffer to receiving
!     processors  ---
!
      IF( NCGN(6,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,5
          MCS = MCS + NCGN(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGN(6,ID+1) = ',NCGN(6,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGN(6,ID+1)
          SBFT(NCS+1) = U_GM(2,NLSGN(M+MCS))
          SBFT(NCS+2) = V_GM(2,NLSGN(M+MCS))
          SBFT(NCS+3) = W_GM(2,NLSGN(M+MCS))
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGN(6,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGN(6,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Top Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFT,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post Top Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending top ghost FE nodes,
!     and unload receiving buffer  ---
!
      IF( NCGN(1,ID+1).GT.0 ) THEN
        NRCVX = NCGN(1,ID+1)*NPVX
        IDSNDX = NPGN(1,ID+1) - 1
        IDRCVX = NPGN(6,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Top Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFT,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post Top Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
!        PRINT *,' MCR = ',MCR,' NCGN(1,ID+1) = ',NCGN(1,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGN(1,ID+1)
          U_GM(2,NLRGN(M+MCR)) = RBFT(NCR+1)
          V_GM(2,NLRGN(M+MCR)) = RBFT(NCR+2)
          W_GM(2,NLRGN(M+MCR)) = RBFT(NCR+3)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DISP_GN_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ESDM( PXIX,PETAX,PMUX,XEX,YEX,ZEX,BX,DAX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Compute the 6 x 24 components of the strain displacement matrix
!     for a hexahedron FE at the location in the natural
!     coordinate system Xi, Eta, and Mu
!
!     PXIX input value of natural coordinate Xi -1 to +1
!     PETAX input value of natural coordinate Eta -1 to +1
!     PMUX input value of natural coordinate Mu -1 to +1
!     XE(1:8) input value of x global Cartesian hexahedron nodes
!     YE(1:8) input value of y global Cartesian hexahedron nodes
!     ZE(1:8) input value of z global Cartesian hexahedron nodes
!
!     where nodes are numbered in the order around the hexahedron
!     grid cells as
!
!     1 i,j,k
!     2 i+1,j,k
!     3 i,j+1,k
!     4 i+1,j+1,k
!     5 i,j,k+1
!     6 i+1,j,k+1
!     7 i,j+1,k+1
!     8 i+1,j+1,k+1  ---
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2016 (Maiden flight of
!     Howard Hughes Spruce Goose).
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
      REAL*8 XIX(8),ETAX(8),MUX(8)
      REAL*8 XEX(8),YEX(8),ZEX(8)
      REAL*8 SAX(3,3),AIX(3,3),BX(6,24)
!
!----------------------Data Statements---------------------------------!
!
      DATA XIX / -1.D+0,1.D+0,-1.D+0,1.D+0,-1.D+0,1.D+0,-1.D+0,1.D+0 /
      DATA ETAX / -1.D+0,-1.D+0,1.D+0,1.D+0,-1.D+0,-1.D+0,1.D+0,1.D+0 /
      DATA MUX / -1.D+0,-1.D+0,-1.D+0,-1.D+0,1.D+0,1.D+0,1.D+0,1.D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ESDM'
!
!---  Initialize a  ---
!
      DO J = 1,3
        DO I = 1,3
          SAX(I,J) = 0.D+0
        ENDDO
      ENDDO
!
!---  Loop over hexahedron nodes  ---
!
      DO M = 1,8
!
!---    Xi multiplier  ---
!
        VXIX = XIX(M)*(1.D+0+PETAX*ETAX(M))*(1.D+0+PMUX*MUX(M))
!
!---    Eta multiplier  ---
!
        VETAX = ETAX(M)*(1.D+0+PXIX*XIX(M))*(1.D+0+PMUX*MUX(M))
!
!---    Mu multiplier  ---
!
        VMUX = MUX(M)*(1.D+0+PETAX*ETAX(M))*(1.D+0+PXIX*XIX(M))
!
!---    dx/dxi at Xi, Eta, and Mu  ---
!
        SAX(1,1) = SAX(1,1) + XEX(M)*VXIX
!
!---    dy/dxi at Xi, Eta, and Mu  ---
!
        SAX(1,2) = SAX(1,2) + YEX(M)*VXIX
!
!---    dz/dxi at Xi, Eta, and Mu  ---
!
        SAX(1,3) = SAX(1,3) + ZEX(M)*VXIX
!
!---    dx/deta at Xi, Eta, and Mu  ---
!
        SAX(2,1) = SAX(2,1) + XEX(M)*VETAX
!
!---    dy/deta at Xi, Eta, and Mu  ---
!
        SAX(2,2) = SAX(2,2) + YEX(M)*VETAX
!
!---    dz/deta at Xi, Eta, and Mu  ---
!
        SAX(2,3) = SAX(2,3) + ZEX(M)*VETAX
!
!---    dx/dmu at Xi, Eta, and Mu  ---
!
        SAX(3,1) = SAX(3,1) + XEX(M)*VMUX
!
!---    dy/dmu at Xi, Eta, and Mu  ---
!
        SAX(3,2) = SAX(3,2) + YEX(M)*VMUX
!
!---    dz/dmu at Xi, Eta, and Mu  ---
!
        SAX(3,3) = SAX(3,3) + ZEX(M)*VMUX
      ENDDO
      DO J = 1,3
        DO I = 1,3
          SAX(I,J) = 1.25D-1*SAX(I,J)
        ENDDO
      ENDDO
!
!---  Inverse Jacobian matrix  ---
!
      AIX(1,1) = SAX(2,2)*SAX(3,3) - SAX(2,3)*SAX(3,2)
      AIX(1,2) = SAX(1,3)*SAX(3,2) - SAX(1,2)*SAX(3,3)
      AIX(1,3) = SAX(1,2)*SAX(2,3) - SAX(1,3)*SAX(2,2)
      AIX(2,1) = SAX(3,1)*SAX(2,3) - SAX(2,1)*SAX(3,3)
      AIX(2,2) = SAX(1,1)*SAX(3,3) - SAX(1,3)*SAX(3,1)
      AIX(2,3) = SAX(1,3)*SAX(2,1) - SAX(1,1)*SAX(2,3)
      AIX(3,1) = SAX(2,1)*SAX(3,2) - SAX(3,1)*SAX(2,2)
      AIX(3,2) = SAX(3,1)*SAX(1,2) - SAX(1,1)*SAX(3,2)
      AIX(3,3) = SAX(1,1)*SAX(2,2) - SAX(1,2)*SAX(2,1)
!
!---  Determinant of the Jacobian matrix  ---
!
      DAX = SAX(1,1)*AIX(1,1) + SAX(2,1)*AIX(1,2) + SAX(3,1)*AIX(1,3)
!
!---  Divide by determinant to complete the inverse  ---
!
      DO J = 1,3
        DO I = 1,3
          AIX(I,J) = AIX(I,J)/DAX
        ENDDO
      ENDDO
!
!---  Zero element strain displacement matrix  ---
!
      DO J = 1,24
        DO I = 1,6
          BX(I,J) = 0.D+0
        ENDDO
      ENDDO
!
!---  Element strain displacement matrix at Xi, Eta, and Mu  ---
!
      DO M = 1,8
        MX = (M-1)*3 + 1
        MY = (M-1)*3 + 2
        MZ = (M-1)*3 + 3
!
!---    Partial of shape factor with respect to hexahedron 
!       coordinates (Xi, Eta, and Mu)  ---
!
        DXIX = 1.25D-1*XIX(M)*(1.D+0 + PETAX*ETAX(M))*
     &    (1.D+0 + PMUX*MUX(M))
        DETAX = 1.25D-1*ETAX(M)*(1.D+0 + PMUX*MUX(M))*
     &    (1.D+0 + PXIX*XIX(M))
        DMUX = 1.25D-1*MUX(M)*(1.D+0 + PXIX*XIX(M))*
     &    (1.D+0 + PETAX*ETAX(M))
!
!---    Partial of shape factor with respect to global coordinates
!       (x, y, and z)  ---
!
        DXX = (AIX(1,1)*DXIX + AIX(1,2)*DETAX + AIX(1,3)*DMUX)
        DYX = (AIX(2,1)*DXIX + AIX(2,2)*DETAX + AIX(2,3)*DMUX)
        DZX = (AIX(3,1)*DXIX + AIX(3,2)*DETAX + AIX(3,3)*DMUX)
!
!---    du/dx for strain epsilon xx  ---
!
        BX(1,MX) = BX(1,MX) + DXX
!
!---    dv/dy for strain epsilon yy  ---
!
        BX(2,MY) = BX(2,MY) + DYX
!
!---    dw/dz for strain epsilon zz  ---
!
        BX(3,MZ) = BX(3,MZ) + DZX
!
!---    (dv/dz + dw/dy) for strain 2*epsilon yz  ---
!
        BX(4,MY) = BX(4,MY) + DZX
        BX(4,MZ) = BX(4,MZ) + DYX
!
!---    (du/dz + dw/dx) for strain 2*epsilon xz  ---
!
        BX(5,MX) = BX(5,MX) + DZX
        BX(5,MZ) = BX(5,MZ) + DXX
!
!---    (du/dy + dv/dx) for strain 2*epsilon xy  ---
!
        BX(6,MX) = BX(6,MX) + DYX
        BX(6,MY) = BX(6,MY) + DXX
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ESDM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE GIESM( EX,ESMX,PTOX,PTVX,XEX,YEX,ZEX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Gauss integration of the element stiffness matrix.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 November 2016 (Ensisheim meteorite).
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
      REAL*8 EX(6,6),ESMX(24,24)
      REAL*8 XEX(8),YEX(8),ZEX(8)
      REAL*8 BTEBX(24,24),BTEX(24,6)
      REAL*8 BX(6,24),BTX(24,6)
      REAL*8 WX(5),PX(5)
      REAL*8 BTPTX(24),PTOX(6),PTVX(24)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/GIESM'
!
!---  Initialization  ---
!
      DO M = 1,5
        WX(M) = 0.D+0
        PX(M) = 0.D+0
      ENDDO
      DO L = 1,6
        DO M = 1,24
          BX(L,M) = 0.D+0
        ENDDO
      ENDDO
      DO M = 1,24
        PTVX(M) = 0.D+0
      ENDDO
!
!---  Number of sample points  ---
!
      IF( ISLC(72).EQ.2 ) THEN
        WX(1) = 1.D+0
        WX(2) = 1.D+0
        PX(1) = -1.D+0/SQRT(3.D+0)
        PX(2) = 1.D+0/SQRT(3.D+0)
      ELSEIF( ISLC(72).EQ.3 ) THEN
        WX(1) = 5.D+0/9.D+0
        WX(2) = 8.D+0/9.D+0
        WX(3) = 5.D+0/9.D+0
        PX(1) = -SQRT(3.D+0/5.D+0)
        PX(2) = 0.D+0
        PX(3) = SQRT(3.D+0/5.D+0)
      ELSEIF( ISLC(72).EQ.4 ) THEN
        WX(1) = 5.D-1 - (1.D+0/6.D+0)*SQRT(5.D+0/6.D+0)
        WX(2) = 5.D-1 + (1.D+0/6.D+0)*SQRT(5.D+0/6.D+0)
        WX(3) = 5.D-1 + (1.D+0/6.D+0)*SQRT(5.D+0/6.D+0)
        WX(4) = 5.D-1 - (1.D+0/6.D+0)*SQRT(5.D+0/6.D+0)
        PX(1) = -SQRT((3.D+0 + 2.D+0*SQRT(6.D+0/5.D+0))/7.D+0)
        PX(2) = -SQRT((3.D+0 - 2.D+0*SQRT(6.D+0/5.D+0))/7.D+0)
        PX(3) = SQRT((3.D+0 - 2.D+0*SQRT(6.D+0/5.D+0))/7.D+0)
        PX(4) = SQRT((3.D+0 + 2.D+0*SQRT(6.D+0/5.D+0))/7.D+0)
      ELSEIF( ISLC(72).EQ.5 ) THEN
        WX(1) = (3.22D+2 - 13.D+0*SQRT(7.D+1))/9.D+2
        WX(2) = (3.22D+2 + 13.D+0*SQRT(7.D+1))/9.D+2
        WX(3) = 5.12D+2/9.D+2
        WX(4) = (3.22D+2 + 13.D+0*SQRT(7.D+1))/9.D+2
        WX(5) = (3.22D+2 - 13.D+0*SQRT(7.D+1))/9.D+2
        PX(1) = -(1.D+0/3.D+0)*SQRT((5.D+0 + 2.D+0*SQRT(1.D+1/7.D+0)))
        PX(2) = -(1.D+0/3.D+0)*SQRT((5.D+0 - 2.D+0*SQRT(1.D+1/7.D+0)))
        PX(3) = 0.D+0
        PX(4) = (1.D+0/3.D+0)*SQRT((5.D+0 - 2.D+0*SQRT(1.D+1/7.D+0)))
        PX(5) = (1.D+0/3.D+0)*SQRT((5.D+0 + 2.D+0*SQRT(1.D+1/7.D+0)))
      ENDIF
!
!---  Numerical integration over hexahedra, loop over the three 
!     conical coordinate directions  ---
!
      DO K = 1,ISLC(72)
        DO J = 1,ISLC(72)
          DO I = 1,ISLC(72)
            CALL ESDM( PX(I),PX(J),PX(K),XEX,YEX,ZEX,BX,DETJX )
            WTX = WX(K)*WX(J)*WX(I)*DETJX
            CALL MATTRP( BX,BTX,6,24 )
            CALL MATMUL( BTX,EX,BTEX,24,6,6 )
            CALL MATMUL( BTX,PTOX,BTPTX,24,6,1 )
            CALL MATMUL( BTEX,BX,BTEBX,24,6,24 )
            DO M = 1,24
              DO L = 1,24
                ESMX(L,M) = ESMX(L,M) + WTX*BTEBX(L,M)  
              ENDDO
              PTVX(M) = PTVX(M) + WTX*BTPTX(M)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of GIESM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE GSM( ESMX,GBDFX,PTVX,N )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Global stiffness matrix
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 November 2016 (Ensisheim meteorite).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE GEO_MECH
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
      REAL*8 EX(6,6),ESMX(24,24)
      REAL*8 PTOX(6),PTVX(24)
      REAL*8 PROP_GMX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/GSM'
!
!---  Hydrate C-Factor composite model  ---
!
      IF( IHCM_GM(N).EQ.1 ) THEN
        PROP_GMX(1) = PROP_GM(1,N) + 
     &    PROP_GM(5,N)*SH(2,N)*PROP_GM(6,N)
      ELSE
        PROP_GMX(1) = PROP_GM(1,N)
      ENDIF
      PROP_GMX(2) = PROP_GM(2,N)
!
!---  Stress-strain matrix  ---
!
      CALL SSM( EX,PROP_GMX )
!
!---  Gauss integration of the element stiffness matrix  ---
!
      DO J = 1,24
        DO I = 1,24
          ESMX(I,J) = 0.D+0
        ENDDO
        PTVX(J) = 0.D+0
      ENDDO
!
!---  Combined poroelasticity + thermoelasticity vector  ---
!
      DO M = 1,6
        PTOX(M) = 0.D+0
      ENDDO
!
!---  Bypass restart check for geomechanics options  ---
!
      ISLC50X = ABS(ISLC(50))
!
!---  Bulk modulus times 3  ---
!
      BLK3X = PROP_GMX(1)/(1.0D+0-2.0D+0*PROP_GMX(2))
!
!---  Poroelasticity vector  ---
!
      IF( MOD(ISLC50X,10).NE.2 .AND. MOD(ISLC50X,10).NE.4 ) THEN
        PX = MAX( PG(2,N),PL(2,N),PN(2,N) ) + PATM
        PTOX(1) = PTOX(1) - PROP_GM(3,N)*(PX-PCMP(N))
        PTOX(2) = PTOX(2) - PROP_GM(3,N)*(PX-PCMP(N))
        PTOX(3) = PTOX(3) - PROP_GM(3,N)*(PX-PCMP(N))
      ENDIF
!
!---  Thermoelasticity vector  ---
!
      IF( MOD(ISLC50X,10).NE.3 .AND. MOD(ISLC50X,10).NE.4 ) THEN
        PTOX(1) = PTOX(1) - PROP_GM(4,N)*(T(2,N)-TCMP(N))/BLK3X
        PTOX(2) = PTOX(2) - PROP_GM(4,N)*(T(2,N)-TCMP(N))/BLK3X
        PTOX(3) = PTOX(3) - PROP_GM(4,N)*(T(2,N)-TCMP(N))/BLK3X
      ENDIF
!
!---  Finite-element node gravity body force  ---
!
      IF( ISLC50X.LT.10 ) THEN
        GBDFX = 1.25D-1*VOL(N)*(1.D+0-POR(2,N))*RHOS(N)*GRAV_GM
      ELSE
        GBDFX = 0.D+0
      ENDIF
!
!---  Gauss integration of the element stiffness matrix.
!
      CALL GIESM( EX,ESMX,PTOX,PTVX,XE(1,N),YE(1,N),ZE(1,N) )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of GSM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCB_SV_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Set values of the geomechanics Jacobian matrix
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2022
!







!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
!

!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GEO_MECH
      USE FDVP
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
!  Include file for Fortran use of the type(tVec) package in PETSc
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
!----------------------Type Declarations-------------------------------!
!




      integer(kind=selected_int_kind(5)) :: IERR

!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(NNZGR) :: BUFFER
      INTEGER, DIMENSION(NNZGR) :: ICOL
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCB_SV_GM'

!
!---  PETSc solver  ---
!
      NC = 0
      DO NEQ = 1,NUKGL(ID+1)
        MC = 0
        DO NCOL = NLU_GM(NEQ),(NLU_GM(NEQ+1)-1)
          NC = NC + 1
          MC = MC + 1
          BUFFER(MC) = DLU(NC)
          ICOL(MC) = MLU_GM(NC)
        ENDDO
        IROW = NEQ + NUKGO(ID+1) - 1
        CALL MatSetValues( G_MAT,1,IROW,MC,ICOL,BUFFER,
     &    INSERT_VALUES,IERR )
      ENDDO

!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCB_SV_GM group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBL_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Load the Jacobian matrix
!
!     Each finite-element node point has three unknowns, the global
!     cartesian component of x-, y-, and z-displacments.
!
!     Equations are ordered by their finite-element node numbers.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 November 2016.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GEO_MECH







!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP

!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)





!
!
!  Include file for Fortran use of the type(tVec) package in PETSc
!


!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 ESMX(24,24),PTVX(24)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX




      integer(kind=selected_int_kind(5)) :: IERR

!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBL_GM'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    Compute element contribution to Jacobian matrix ---
!
        CALL GSM( ESMX,GBDFX,PTVX,N )
!
!---    Loop over the FE Nodes of the element  ---
!
        DO M1 = 1,8
!
!---      Global FE node number  ---
!
          NFEN = ND_GM( M1,N )
          IF( NFEN.EQ.0 ) CYCLE
          IF( IGHN(NFEN).EQ.1 ) CYCLE
          MX = (IM_GM(NFEN)-1)*3 + 1
          NEPE = NLU_GM(MX+1-IEQ_OFFSET_GM) - NLU_GM(MX-IEQ_OFFSET_GM)
!
!---      Loop over the displacement equations  ---
!
          DO M2 = 1,3
!
!---        Equation number  ---
!
            KROW = (M1-1)*3 + M2
            NMD = (IM_GM(NFEN)-1)*3 + M2
!
!---        Loop over the FE Nodes of the element  ---
!
            DO M3 = 1,8
!
!---          Loop over the displacement equations  ---
!
              DO M4 = 1,3
!
!---            Column number  ---
!
                KCOL = (M3-1)*3 + M4
!
!---            Row number  ---
!
                M5 = NK_GM(M3,M1)
                MCOL = (M2-1)*NEPE + (KLU_GM(M5,NFEN)-1) + M4
                DLU(MCOL) = DLU(MCOL) + ESMX(KROW,KCOL)
              ENDDO
            ENDDO
!
!---        Poroelasticity + thermoelasticity  ---
!
            BUFFER = -PTVX(KROW)





            IROW_P = NMD - 1
            CALL VecSetValues( G_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,
     &        IERR )

!
!---        Gravitational body force  ---
!
            IF( M2.EQ.3 ) THEN
              BUFFER = -GBDFX





              IROW_P = NMD - 1
              CALL VecSetValues( G_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,
     &          IERR )

            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBL_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBP_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Configure the Jacobian matrix pointer arrays for geomechanics.
!
!     Each finite-element node point has three unknowns, the global
!     cartesian component of x-, y-, and z-displacments.
!
!     Equations are ordered by their finite-element node numbers.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 5 October 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GEO_MECH
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER ISTCX(27)
      INTEGER NKX(8,8)
      INTEGER NUKGLX(NP)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
      CHARACTER(32) :: CHMSG
!
!----------------------Data Statements---------------------------------!
!
      DATA NKX /14,15,17,18,23,24,26,27,
     &          13,14,16,17,22,23,25,26,
     &          11,12,14,15,20,21,23,24,
     &          10,11,13,14,19,20,22,23,
     &          5,6,8,9,14,15,17,18,
     &          4,5,7,8,13,14,16,17,
     &          2,3,5,6,11,12,14,15,
     &          1,2,4,5,10,11,13,14 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBP_GM'
!
!---  Initializing global counters for the number of geomechanics
!     equations  ---
!
      DO N = 1,NP
        NUKGLX(N) = 0
        NUKGL(N) = 0
      ENDDO
!!
!!---  Loop over number of local finite-element nodes to load
!!     the Jacobian matrix pointer array  ---
!!
!      DO NFEN = 1,NFNGN(ID+1)
!        IM_GM(NFEN) = 0
!!
!!---    Finite-element node is inactive, only if all elements
!!       that contain the FE node are inactive grid cells  ---
!!
!        IFX = 0
!        DO M = 1,8
!          N = NE_GM(M,NFEN)
!          IF( N.EQ.0 ) CYCLE
!          IF( IXP(N).NE.0 ) THEN
!            IFX = 1
!            EXIT
!          ENDIF
!        ENDDO
!        IF( IFX.EQ.1 ) THEN
!          IM_GM(NFEN) = NF_GM(NFEN)
!        ENDIF
!      ENDDO
!
!---  Loop over number of local finite-element nodes to load
!     the Jacobian matrix pointer array  ---
!
      NROW = 0
      NMDMN = NFNGN_G*3
      DO NFEN = 1,NFNGN(ID+1)
!
!---    Skip ghost finite-element nodes  ---
!
        IF( IGHN(NFEN).EQ.1 ) CYCLE
!
!---    Finite-element node is inactive, only if all elements
!       that contain the FE node are inactive grid cells  ---
!
        IFX = 0
        DO M = 1,8
          N = NE_GM(M,NFEN)
          IF( N.EQ.0 ) CYCLE
          IF( IXP(N).NE.0 ) THEN
            IFX = 1
            EXIT
          ENDIF
        ENDDO
        IF( IFX.EQ.1 ) THEN
          NROW = NROW + 3
          NMD = (IM_GM(NFEN)-1)*3
          DO L = 1,3
            NMDMN = MIN( NMDMN,NMD+L)
          ENDDO
        ENDIF
      ENDDO
      NUKGLX(ID+1) = NROW
      CALL MPI_ALLREDUCE( NUKGLX,NUKGL,NP,MPI_INTEGER,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
      NUKGG = 0
      DO N = 1,NP
        NUKGG = NUKGG + NUKGL(N)
      ENDDO
      NUKGO(1) = 0
      DO N = 2,NP
        NUKGO(N) = 0
        DO M = 2,N
          NUKGO(N) = NUKGO(N) + NUKGL(M-1)
        ENDDO
      ENDDO
!      DO N = 1,NP
!        PRINT *,'NUKGL(',N,') = ',NUKGL(N),'NUKGO(',N,') = ',NUKGO(N),
!     &    'NUKGG = ',NUKGG,'ID = ',ID
!      ENDDO
!
!---  Geomechanics equation local offsets ---
!
      IEQ_OFFSET_GM = NMDMN - 1
!
!---  SPLIB .or. Lis .or. PETSc Solver  ---
!
      NC = 0
!
!---  Loop over number of finite-element nodes  ---
!
      DO NFEN = 1,NFNGN(ID+1)
!
!---    Skip ghost finite-element nodes  ---
!
        IF( IGHN(NFEN).EQ.1 ) CYCLE
!
!---    Determine the connection stencil for the finite-element
!       node  ---
!
        DO M1 = 1,27
          ISTCX(M1) = 0
        ENDDO
        LJO_GM = 1
!
!---    Loop over the connections to adjacent finite elements  ---
!
        DO M2 = 1,8
          N = NE_GM(M2,NFEN)
          IF( N.EQ.0 ) CYCLE
!
!---      Loop over the FE nodes in the adjacent finite element  ---
!
          DO M3 = 1,8
            MFEN = ND_GM(M3,N)
            IF( MFEN.EQ.0 ) CYCLE
            ISTCX(NKX(M2,M3)) = MFEN
          ENDDO
        ENDDO
!
!---    Loop over displacement equations  ---
!
        DO M4 = 1,3
          NMD = (IM_GM(NFEN)-1)*3 + M4
!
!---      Loop over the active connections  ---
!
          DO M5 = 1,27
            IF( ISTCX(M5).EQ.0 ) CYCLE
            MFEN = ISTCX(M5)
!
!---        Loop over displacement equations  ---
!
            DO M6 = 1,3
              NC = NC + 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      NNZ_GM = MAX( NC,LJO_GM )
!
!---  Reallocate memory for the Jacobian matrix array (DLU) if
!     more memory is needed for geomechanics ---
!
!      PRINT *,'NNZ = ',NNZ,'NNZ_GM = ',NNZ_GM,'NUKGL(ID+1) = ',
!     &  NUKGL(ID+1),'NUKFL(ID+1) = ',NUKFL(ID+1),'ID = ',ID
      IF( NNZ_GM.GT.NNZ ) THEN
        IF( ALLOCATED(DLU) ) THEN
          DEALLOCATE( DLU,STAT=ISTAT )
          CHMSG = 'DLU'
          CALL DEALLOC_ERROR( CHMSG,ISTAT )
        ENDIF
        ALLOCATE( DLU(1:NNZ_GM),STAT=ISTAT )
        CHMSG = 'DLU'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
      ENDIF

!
!---  Reallocate memory for the problem/solution vector array (BLU) if
!     more memory is needed for geomechanics ---
!
      NUKFX = NFCGC(ID+1)*ISVC
      IF( ID+1.EQ.NP ) NUKFX = NUKFX + N_CW
      NUKGX = NFNGN(ID+1)
      IF( NUKGX.GT.NUKFX ) THEN     
        IF( ALLOCATED(BLU) ) THEN
          DEALLOCATE( BLU,STAT=ISTAT )
          CHMSG = 'BLU'
          CALL DEALLOC_ERROR( CHMSG,ISTAT )
        ENDIF
        ALLOCATE( BLU(1:NUKGX),STAT=ISTAT )
        CHMSG = 'BLU'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
      ENDIF

!
!---  Allocate memory for the CSR pointer (NLU_GM), index (MLU_GM),
!     and equation to index pointer (KLU_GM) ---
!
      ALLOCATE( NLU_GM(1:NROW+1),STAT=ISTAT )
      CHMSG = 'NLU_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( MLU_GM(1:NNZ_GM),STAT=ISTAT )
      CHMSG = 'MLU_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( KLU_GM(1:27,1:NFNGN(ID+1)),STAT=ISTAT )
      CHMSG = 'KLU_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      NC = 0
      NLU_GM(1) = 0
      NNZGR = 0
!
!---  Loop over number of local finite-element nodes to load
!     the Jacobian matrix pointer array  ---
!
      DO NFEN = 1,NFNGN(ID+1)
!
!---    Skip ghost finite-element nodes  ---
!
        IF( IGHN(NFEN).EQ.1 ) CYCLE
!
!---    Determine the connection stencil for the finite-element
!       node  ---
!
        DO M1 = 1,27
          ISTCX(M1) = 0
        ENDDO
!
!---    Loop over the connections to adjacent finite elements  ---
!
        DO M2 = 1,8
          N = NE_GM(M2,NFEN)
          IF( N.EQ.0 ) CYCLE
!
!---      Loop over the FE nodes in the adjacent finite element  ---
!
          DO M3 = 1,8
            MFEN = ND_GM(M3,N)
            IF( MFEN.EQ.0 ) CYCLE
            ISTCX(NK_GM(M3,M2)) = MFEN
          ENDDO
        ENDDO
!
!---    Loop over displacement equations  ---
!
        DO M4 = 1,3
          NMD = (IM_GM(NFEN)-1)*3 + M4
!
!---      Loop over the active connections  ---
!
          DO M5 = 1,27
            IF( ISTCX(M5).EQ.0 ) CYCLE
            MFEN = ISTCX(M5)
!
!---        Loop over displacement equations  ---
!
            DO M6 = 1,3
              NC = NC + 1
              IF( M4*M6.EQ.1 ) KLU_GM(M5,NFEN) = NC
              MLU_GM(NC) = (IM_GM(MFEN)-1)*3 + M6 - 1
            ENDDO
          ENDDO
          NLU_GM(NMD+1-IEQ_OFFSET_GM) = NC
          NNZGR = MAX( NNZGR,(NC-NLU_GM(NMD-IEQ_OFFSET_GM)) )
        ENDDO
      ENDDO
      MK_GM = NC
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBP_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBZ_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Zero the geomechanics Jacobian matrix and solution/problem
!     vector.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 October 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GEO_MECH
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
      SUB_LOG(ISUB_LOG) = '/JCBTZ_GM'

!
!---  Zero problem vector array for solute transport  ---
!
      NUKGX = NFNGN(ID+1)
      DO N = 1,NUKGX
        BLU(N) = 0.D+0
      ENDDO

!
!---  Zero Jacobian matrix array for solute transport  ---
!
      DO N = 1,NNZ_GM
        DLU(N) = 0.D+0
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTZ_GM group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LD_GM( INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     INDX = 1
!     Load old time step values of volumetric strain and pore pressure.
!     INDX = 2
!     Load k values of volumetric strain and pore pressure.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 1 March 2017 (Yellowstone National 
!     Park is established as the world's first national park. in 1872).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE GEO_MECH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LD_GM'
!
!---  Loop over the finite elements (i.e., grid cells) ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
        EPSV_GM(INDX,N) = EPSV_GM(INDX+1,N)
        SIGV_GM(INDX,N) = SIGV_GM(INDX+1,N)
        P_GM(INDX,N) = P_GM(INDX+1,N)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LD_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDDISP_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Load reference displacements at finite element nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 February 2017 (Wake Forest 
!     University is established in 1834).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE GEO_MECH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LDDISP_GM'
!
!---  Loop over the finite-element nodes, including FE ghost nodes ---
!
      DO NFEN = 1,NFNGN(ID+1)
        U_GM(1,NFEN) = U_GM(2,NFEN)
        V_GM(1,NFEN) = V_GM(2,NFEN)
        W_GM(1,NFEN) = W_GM(2,NFEN)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDDISP_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE MATTRP( AX,BX,IX,JX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Return the transpose of matrix AX (IX x JX) in BX (JX x IX).
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 November 2016 (Ensisheim meteorite).
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
      REAL*8 AX(IX,JX),BX(JX,IX)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/MATTRP'
!
!---  Matrix transpose  ---
!
      DO J = 1,JX
        DO I = 1,IX
          BX(J,I) = AX(I,J)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of MATTRP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE MATMUL( AX,BX,CX,IX,JX,KX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Multiple the matrix AX (IX X JX) and BX (JX X KX), returning the
!     result in matrix CX (IX X KX).
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 November 2016 (Ensisheim meteorite).
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AX(IX,JX),BX(JX,KX),CX(IX,KX)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/MATMUL'
!
!---  Initialize CX  ---
!
      DO K = 1,KX
        DO I = 1,IX
          CX(I,K) = 0.D+0
        ENDDO
      ENDDO
!
!---  Matrix multiply  ---
!
      DO I = 1,IX
        DO K = 1,KX
          DO J = 1,JX
            CX(I,K) = CX(I,K) + AX(I,J)*BX(J,K)
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of MATMUL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PERMRF_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Permeability as a function of mean effective stress via porosity.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 27 March 2017 (Construction of the 
!     Trans-Alaska Pipeline System begins in 1975).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE JACOB
      USE GRID
      USE GEO_MECH
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PERMRF_GM'
!
!---  Loop over the finite elements (i.e., grid cells) ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    Davis-Davis intrinsic permeability-porosity function  ---
!
        IF( IPROP_GM(2,N).EQ.1 ) THEN
!
!---      Loop over increment indices  ---
!
          DO M = 2,ISVC+2
            PERMRF(M,N) = EXP(PROP_GM(7,N)*
     &       MIN(((PORD(M,N)/POR(2,N))-1.D+0),0.D+0))
          ENDDO
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PERMRF_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PRESS_GM( INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Set k iterate value of pore pressure.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 4 October 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE GEO_MECH
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
      SUB_LOG(ISUB_LOG) = '/PRESS_GM'
!
!---  Loop over the finite elements ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        P_GM(INDX,N) = MAX( PG(2,N),PL(2,N),PN(2,N) ) + PATM
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PRESS_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_BIN_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Read binary geomechanics files.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 6 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE PROP
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
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_BIN_GM'
!
!---  Read gmgd.bin for geomechanical grid data  ---
!
      IF( ISLC(50).NE.0 ) CALL READ_GMGD
!      PRINT *,'Post READ_GMBC: ID = ',ID
!
!---  Read gmec.bin for geomechanical property data onto 
!     field nodes  ---
!
      IF( ISLC(50).NE.0 ) CALL READ_GMEC
!      PRINT *,'Post READ_GMEC: ID = ',ID
!
!---  Read gmbc.bin for geomechanical boundary condition data  ---
!
      IF( ISLC(50).NE.0 ) CALL READ_GMBC
!      PRINT *,'Post READ_GMBC: ID = ',ID
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_BIN_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_GMBC
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Read binary gmbc.bin file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE PROP
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
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_GMBC'
!
!---  Set local starting point for local copies of nodal variables  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NFCGC(I)
      ENDDO
!
!---  Open gmbc.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'gmbc.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Allocate memory for NBC_GM  ---
!
      ALLOCATE( NBC_GM(1:NP),STAT=ISTAT )
      CHMSG = 'NBC_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read number of geomechanical boundary condition nodes on each
!     processor  ---
!
      NVAR = NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NBC_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Allocate memory for geomechanics boundary condition arrays  ---
!
      CALL ALLOC_GMBC
!      PRINT *,'Post ALLOC_GMBC: ID = ',ID
!
!---  Initialize memory for geomechanics boundary condition arrays  ---
!
      CALL INTLZ_GMBC
!      PRINT *,'Post INTLZ_GMBC: ID = ',ID
!
!---  Read geomechanical boundary condition variables 
!     (duplicated across processors)  ---
!
      NVAR = LBCV_GM*LBTM_GM*LBCIN_GM
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,BC_GM,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read reference geomechanical boundary condition index 
!     (duplicated across processors)  ---
!
      NVAR = LBCIN_GM
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCR_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Assign offsets for reading geomechanical boundary condition
!      indices on processors  ---
!
      NC_GM = 0
      DO I = 1,ID
        NC_GM = NC_GM + NBC_GM(I)
      ENDDO
      NCP_GM = 0
      DO N = 1,NP
        NCP_GM = NCP_GM + NBC_GM(N)
      ENDDO
!      PRINT *,'NC_GM = ',NC_GM,'NCP_GM = ',NCP_GM,'ID = ',ID
!
!---  Index array of geomechanical boundary condition finite-element
!     nodes  ---
!
      NVAR = NBC_GM(ID+1)
      OFFSET = IOFFSET + NC_GM*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP_GM*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCN_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'Post IBCN_GM: ID = ',ID
!      PRINT *,'IBCN_GM(1) = ',IBCN_GM(1),
!     &  'IBCN_GM(',NBC_GM(ID+1),') = ',IBCN_GM(NBC_GM(ID+1)),'ID = ',ID
!
!---  Index array of geomechanical boundary condition directions  ---
!
      NVAR = NBC_GM(ID+1)
      OFFSET = IOFFSET + NC_GM*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP_GM*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCD_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'Post IBCD_GM: ID = ',ID
!      PRINT *,'IBCD_GM(1) = ',IBCD_GM(1),
!     &  'IBCD_GM(',NBC_GM(ID+1),') = ',IBCD_GM(NBC_GM(ID+1)),'ID = ',ID
!
!---  Index array of geomechanical boundary condition number of 
!     time points  ---
!
      NVAR = NBC_GM(ID+1)
      OFFSET = IOFFSET + NC_GM*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP_GM*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCM_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'Post IBCM_GM: ID = ',ID
!      PRINT *,'IBCM_GM(1) = ',IBCM_GM(1),
!     &  'IBCM_GM(',NBC_GM(ID+1),') = ',IBCM_GM(NBC_GM(ID+1)),'ID = ',ID
!
!---  Index array of geomechanical boundary condition input links  ---
!
      NVAR = NBC_GM(ID+1)
      OFFSET = IOFFSET + NC_GM*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP_GM*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCIN_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'Post IBCIN_GM: ID = ',ID
!      PRINT *,'IBCIN_GM(1) = ',IBCIN_GM(1),
!     &  'IBCIN_GM(',NBC_GM(ID+1),') = ',IBCIN_GM(NBC_GM(ID+1)),
!     &  'ID = ',ID
!
!---  Index array of geomechancial boundary condition cycling 
!     options  ---
!
      NVAR = NBC_GM(ID+1)
      OFFSET = IOFFSET + NC_GM*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP_GM*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCC_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'Post IBCC_GM: ID = ',ID
!      PRINT *,'IBCC_GM(1) = ',IBCC_GM(1),
!     &  'IBCC_GM(',NBC_GM(ID+1),') = ',IBCC_GM(NBC_GM(ID+1)),
!     &  'ID = ',ID
!
!---  Index array of geomechanical boundary condition types  ---
!
      NVAR = NBC_GM(ID+1)*3
      OFFSET = IOFFSET + NC_GM*3*NBYTI + NBYTB
      IOFFSET = IOFFSET + NCP_GM*3*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IBCT_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'IBCT_GM(1,1) = ',IBCT_GM(1,1),
!     &  'IBCT_GM(2,1) = ',IBCT_GM(2,1),
!     &  'IBCT_GM(3,1) = ',IBCT_GM(3,1),
!     &  'IBCT_GM(1,',NBC_GM(ID+1),') = ',IBCT_GM(1,NBC_GM(ID+1)),
!     &  'IBCT_GM(2,',NBC_GM(ID+1),') = ',IBCT_GM(2,NBC_GM(ID+1)),
!     &  'IBCT_GM(3,',NBC_GM(ID+1),') = ',IBCT_GM(3,NBC_GM(ID+1)),
!     &  'ID = ',ID
!      PRINT *,'Post IBCT_GM: ID = ',ID
!
!---  Close the gmbc.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_GMBC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_GMEC
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Read binary gmec.bin file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 September 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE PROP
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
      SUB_LOG(ISUB_LOG) = '/READ_GMEC'
!
!---  Set local starting point for local copies of nodal variables  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NFCGC(I)
      ENDDO
!
!---  Open gmec.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'gmec.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read RSDM_GM array (duplicated across processors)  ---
!
      NVAR = LEPD
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,RSDM_GM,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!
!---  Read local copies of PROP_GM array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*7
      OFFSET = IOFFSET + NBYTB + NC*7*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*7*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,PROP_GM,NVAR,MPI_REAL8,
     &  STATUS,IERR )
!      DO M = 1,7
!        PRINT *,'PROP_GM(',M,',1) = ',PROP_GM(M,1),'ID = ',ID
!      ENDDO
!
!---  Read local copies of IPROP_GM array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)*3
      OFFSET = IOFFSET + NBYTB + NC*3*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*3*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IPROP_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      DO M = 1,3
!        PRINT *,'IPROP_GM(',M,',1) = ',IPROP_GM(M,1),'ID = ',ID
!      ENDDO
!
!---  Read local copies of IHCM_GM array (including ghost cells)  ---
!
      NVAR = NFCGC(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IHCM_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      PRINT *,'IHCM_GM(1) = ',IHCM_GM(1),'ID = ',ID
!
!---  Close the gmec.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_GMEC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_GMGD
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Read binary gmgd.bin file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 September 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE PROP
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
      INTEGER, DIMENSION(:), ALLOCATABLE :: ND_GMX,NE_GMX,NGHNX
      INTEGER NKX(8,8)
      CHARACTER(32) :: CHMSG
!
!----------------------Data Statements---------------------------------!
!
      DATA NKX /14,15,17,18,23,24,26,27,
     &          13,14,16,17,22,23,25,26,
     &          11,12,14,15,20,21,23,24,
     &          10,11,13,14,19,20,22,23,
     &          5,6,8,9,14,15,17,18,
     &          4,5,7,8,13,14,16,17,
     &          2,3,5,6,11,12,14,15,
     &          1,2,4,5,10,11,13,14 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_GMGD'
!
!---  Allocate memory for NUKGL ---
!
      ALLOCATE( NUKGL(1:NP),STAT=ISTAT )
      CHMSG = 'NUKGL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NUKGO ---
!
      ALLOCATE( NUKGO(1:NP),STAT=ISTAT )
      CHMSG = 'NUKGO'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NUKGP ---
!
      ALLOCATE( NUKGP(1:NP),STAT=ISTAT )
      CHMSG = 'NUKGP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for IDP_GM ---
!
      ALLOCATE( IDP_GM(1:2,1:NP),STAT=ISTAT )
      CHMSG = 'IDP_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for JDP_GM ---
!
      ALLOCATE( JDP_GM(1:2,1:NP),STAT=ISTAT )
      CHMSG = 'JDP_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for KDP_GM ---
!
      ALLOCATE( KDP_GM(1:2,1:NP),STAT=ISTAT )
      CHMSG = 'KDP_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Set local starting point for local copies of nodal variables  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NFCGC(I)
      ENDDO
!
!---  Open gmgd.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'gmgd.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRD,IERR )
      IOFFSET = 0
!
!---  Read IDP_GM array (duplicated across processors)  ---
!
      NVAR = 2*NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IDP_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read JDP_GM array (duplicated across processors)  ---
!
      NVAR = 2*NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,JDP_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read KDP_GM array (duplicated across processors)  ---
!
      NVAR = 2*NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,KDP_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Allocate memory for ND_GM  ---
!
      ALLOCATE( ND_GM(1:8,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ND_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for ND_GMX  ---
!
      ALLOCATE( ND_GMX(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ND_GMX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read local copies of ND_GM array (including ghost cells)
!     with ND_GM pointing to global finite-element node numbers  ---
!
      NVAR = NFCGC(ID+1)*8
      OFFSET = IOFFSET + NBYTB + NC*8*NBYTI
      IOFFSET = IOFFSET + NFCGC_G*8*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,ND_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Allocate memory for NFNGN, number of finite-element nodes and
!     ghost finite-element nodes on each processor  ---
!
      ALLOCATE( NFNGN(1:NP),STAT=ISTAT )
      CHMSG = 'NFNGN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read NFNGN array (duplicated across processors)  ---
!
      NVAR = NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET =  IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NFNGN,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
      NFNGN_G = 0
      DO N = 1,NP
        NFNGN_G = NFNGN_G + NFNGN(N)
      ENDDO
!
!---  Allocate memory for geomechanics variables  ---
!
      CALL ALLOC_GMEC
!      PRINT *,'Post ALLOC_GMEC: ID = ',ID
!
!---  Initialize memory for geomechanics variables  ---
!
      CALL INTLZ_GMEC
!      PRINT *,'Post INTLZ_GMEC: ID = ',ID
!
!---  Set local starting point for local copies of finite-element
!     nodal variables  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NFNGN(I)
      ENDDO
!
!---  Read local copies of NE_GM array (including ghost cells)
!     with NE_GM pointing to global field node numbers  ---
!
      NVAR = NFNGN(ID+1)*8
      OFFSET = IOFFSET + NBYTB + NC*8*NBYTI
      IOFFSET = IOFFSET + NFNGN_G*8*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NE_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read local copies of NF_GM array (including ghost cells)  ---
!
      NVAR = NFNGN(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFNGN_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NF_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read local copies of IM_GM array (including ghost cells)  ---
!
      NVAR = NFNGN(ID+1)
      OFFSET = IOFFSET + NBYTB + NC*NBYTI
      IOFFSET = IOFFSET + NFNGN_G*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,IM_GM,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Convert ND_GM to local finite-element node numbers  ---
!
      DO M = 1,8
        DO N = 1,NFCGC(ID+1)
          ND_GMX(N) = 0
          IF( ND_GM(M,N).EQ.0 ) CYCLE
          DO NX = 1,NFNGN(ID+1)
            IF( NF_GM(NX).EQ.ND_GM(M,N) ) THEN
              ND_GMX(N) = NX
              EXIT
            ENDIF
          ENDDO
        ENDDO
        DO N = 1,NFCGC(ID+1)
          ND_GM(M,N) = ND_GMX(N)
        ENDDO
      ENDDO
!      IF( ID.EQ.3 ) THEN
!        DO N = 1,NFCGC(ID+1)
!          DO M = 1,8
!            PRINT *,'ND_GM(',M,',',ND(N),') = ',ND_GM(M,N),
!     &        'NF_GM(ND_GM(M,N)) = ',NF_GM(ND_GM(M,N))
!          ENDDO
!        ENDDO
!      ENDIF
!
!---  Dellocate memory for ND_GMX  ---
!
      DEALLOCATE( ND_GMX,STAT=ISTAT )
      CHMSG = 'ND_GMX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NE_GMX  ---
!
      ALLOCATE( NE_GMX(1:NFNGN(ID+1)),STAT=ISTAT )
      CHMSG = 'NE_GMX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Convert NE_GM to local field node numbers  ---
!
      DO M = 1,8
        DO N = 1,NFNGN(ID+1)
          NE_GMX(N) = 0
          IF( NE_GM(M,N).EQ.0 ) CYCLE
          DO NX = 1,NFCGC(ID+1)
            IF( ND(NX).EQ.NE_GM(M,N) ) THEN
              NE_GMX(N) = NX
              EXIT
            ENDIF
          ENDDO
        ENDDO
        DO N = 1,NFNGN(ID+1)
          NE_GM(M,N) = NE_GMX(N)
        ENDDO
      ENDDO
!
!---  Dellocate memory for NE_GMX  ---
!
      DEALLOCATE( NE_GMX,STAT=ISTAT )
      CHMSG = 'NE_GMX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NGHN, number of ghost finite-element nodes
!     on each processor  ---
!
      ALLOCATE( NGHNX(1:NP),STAT=ISTAT )
      CHMSG = 'NGHNX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NGHN(1:NP),STAT=ISTAT )
      CHMSG = 'NGHN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NCGN ---
!
      ALLOCATE( NCGN(1:6,1:NP),STAT=ISTAT )
      CHMSG = 'NCGN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NPGN ---
!
      ALLOCATE( NPGN(1:6,1:NP),STAT=ISTAT )
      CHMSG = 'NPGN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Initialize finite-element ghost node memory ---
!
      DO N = 1,NP
        NGHNX(N) = 0
        NGHN(N) = 0
        DO M = 1,6
          NCGN(M,N) = 0
          NPGN(M,N) = 0
        ENDDO
      ENDDO
!
!---  Initialize ghost node map  ---
!
      DO L = 1,NFNGN(ID+1)
        IGHN(L) = 0
      ENDDO
!
!---  Identify local ghost finite-element nodes  ---
!
      KFLD_L = KDP_GM(2,ID+1)-KDP_GM(1,ID+1)+1
      JFLD_L = JDP_GM(2,ID+1)-JDP_GM(1,ID+1)+1
      IFLD_L = IDP_GM(2,ID+1)-IDP_GM(1,ID+1)+1
      DO K = KDP_GM(1,ID+1),KDP_GM(2,ID+1)
        DO J = JDP_GM(1,ID+1),JDP_GM(2,ID+1)
          DO I = IDP_GM(1,ID+1),IDP_GM(2,ID+1)
            KX = K - KDP_GM(1,ID+1) + 1
            JX = J - JDP_GM(1,ID+1) + 1
            IX = I - IDP_GM(1,ID+1) + 1
            NFEN = (KX-1)*IFLD_L*JFLD_L + (JX-1)*IFLD_L + IX
!
!---        Bottom ghost nodes  ---
!
            IF( KX.EQ.1 .AND. K.GT.1 ) THEN
              IGHN(NFEN) = 1
            ENDIF
!
!---        South ghost nodes  ---
!
            IF( JX.EQ.1 .AND. J.GT.1 ) THEN
              IGHN(NFEN) = 1
            ENDIF
!
!---        West ghost nodes  ---
!
            IF( IX.EQ.1 .AND. I.GT.1 ) THEN
              IGHN(NFEN) = 1
            ENDIF
!
!---        East ghost nodes  ---
!
            IF( IX.EQ.IFLD_L .AND. I.LT.IDP_GM(2,NP) ) THEN
              IGHN(NFEN) = 1
            ENDIF
!
!---        North ghost nodes  ---
!
            IF( JX.EQ.JFLD_L .AND. J.LT.JDP_GM(2,NP) ) THEN
              IGHN(NFEN) = 1
            ENDIF
!
!---        Top ghost nodes  ---
!
            IF( KX.EQ.KFLD_L .AND. K.LT.KDP_GM(2,NP) ) THEN
              IGHN(NFEN) = 1
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Count ghost finite-element nodes  ---
!
      DO L = 1,NFNGN(ID+1)
        NGHNX(ID+1) = NGHNX(ID+1) + IGHN(L)
      ENDDO
      CALL MPI_ALLREDUCE( NGHNX,NGHN,NP,MPI_INTEGER,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
!      DO L = 1,NFNGN(ID+1)
!        IF( IGHN(L).EQ.1 ) PRINT *,'IGHN(',L,') = ',IGHN(L),'ID = ',ID
!      ENDDO
!      DO M = 1,NP
!        PRINT *,'NGHN(',M,') = ',NGHN(M),'ID = ',ID
!      ENDDO
!
!---  Deallocate temporary local memory  ---
!
      DEALLOCATE( NGHNX,STAT=ISTAT )
      CHMSG = 'NGHNX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!      DO N = 1,NP
!        PRINT *,'NGHN(',N,') = ',NGHN(N),'ID = ',ID
!      ENDDO
!
!---  Set NK_GM on all processors  ---
!
      DO K = 1,8
        DO L = 1,8
          NK_GM(K,L) = NKX(K,L)
        ENDDO
      ENDDO
!
!---  Read NCGN array, number of finite-element ghost nodes to be 
!     sent in the stencil directions, (duplicated across processors)  ---
!
      NVAR = 6*NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NCGN,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read NPGN array, receiving processors in the
!     stencil directions, (duplicated across processors)  ---
!
      NVAR = 6*NP
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NPGN,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Allocate memory for NLSGN ---
!
      ALLOCATE( NLSGN(1:NGHN(ID+1)),STAT=ISTAT )
      CHMSG = 'NLSGN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate memory for NLRGN ---
!
      ALLOCATE( NLRGN(1:NGHN(ID+1)),STAT=ISTAT )
      CHMSG = 'NLRGN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Initialize NLSGN and NLRGN ---
!
      DO M = 1,NGHN(ID+1)
        NLSGN(M) = 0
        NLRGN(M) = 0
      ENDDO
!
!---  Set offsets for FE ghost node sending and receiving
!     indices  ---
!
      NCX = 0
      DO I = 1,ID
        NCX = NCX + NGHN(I)
      ENDDO
      NCSX = 0
      DO I = 1,NP
        NCSX = NCSX + NGHN(I)
      ENDDO
!
!---  Read local copies of NLSGN array  ---
!
      NVAR = NGHN(ID+1)
      OFFSET = IOFFSET + NBYTB + NCX*NBYTI
      IOFFSET = IOFFSET + NCSX*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NLSGN,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!
!---  Read local copies of NLRGN array  ---
!
      NVAR = NGHN(ID+1)
      OFFSET = IOFFSET + NBYTB + NCX*NBYTI
      IOFFSET = IOFFSET + NCSX*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IRD,OFFSET,NLRGN,NVAR,MPI_INTEGER,
     &  STATUS,IERR )
!      DO N = 1,NGHN(ID+1)
!        PRINT *,'NLSGN(',N,') = ',NLSGN(N),'ID = ',ID
!        PRINT *,'NLRGN(',N,') = ',NLRGN(N),'ID = ',ID
!      ENDDO
!
!---  Close the gmgd.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRD,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_GMGD group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Convergence check for sequential coupled flow and transport and
!     geomechanics  
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 February 2017 (Robert II becomes 
!     King of Scotland, beginning the Stuart dynasty in 1371).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE GEO_MECH
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
      REAL*8 PROP_GMX(2)
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL_GM'
!
!---  Initialize residual  ---
!
      RSD_GM = 0.D+0
      RSD_GML = 0.D+0
      NSD_GM = 0
      NSD_GML = 0
!
!---  No poroelasticity option ---
!
      ISLC50X = ABS(ISLC(50))
      IF( MOD(ISLC50X,10).EQ.2 .OR. MOD(ISLC50X,10).EQ.4 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Loop over the finite elements (i.e., grid cells) ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---    Flow porosity at iterate level k+1  --
!
        PORDFWK1X = PORD(2,N)
!
!---    Davis-Davis porosity versus mean effective stress model  ---
!
        IF( IPROP_GM(1,N).EQ.1 ) THEN
!
!---      Volumetric stress (effective), at iterate level k+1  ---
!
          SIGV_GM(3,N) = (SIG_GM(1,N)+SIG_GM(2,N)+SIG_GM(3,N))/3.D+0
!
!---      Mechanical porosity at iterate level k+1  --
!
          PORDMCHK1X = (POR(1,N)-PROP_GM(5,N))*
     &      EXP(-PROP_GM(6,N)*SIGV_GM(3,N)) +  PROP_GM(5,N)
!
!---    Classical model  ---
!
        ELSE
!
!---      Hydrate C-Factor composite model  ---
!
          IF( IHCM_GM(N).EQ.1 ) THEN
            PROP_GMX(1) = PROP_GM(1,N) + 
     &        PROP_GM(5,N)*SH(2,N)*PROP_GM(6,N)
          ELSE
            PROP_GMX(1) = PROP_GM(1,N)
          ENDIF
          PROP_GMX(2) = PROP_GM(2,N)
!
!---      Drained bulk modulus  --
!
          BLKDRNX = PROP_GMX(1)/(3.0E+0*(1.0D+0-2.0D+0*PROP_GMX(2)))
!
!---      1/N  ---
!
          OONMODX = (PROP_GM(3,N)-POR(1,N))*
     &      (1.0D+0-PROP_GM(3,N))/BLKDRNX
!
!---      Volumetric strain differential, at iterate level k  ---
!
          DEVKX = EPSV_GM(2,N)-EPSV_CMP(N)
!
!---      Pressure differential, at iterate level k  ---
!
          DPKX = P_GM(2,N)-PCMP(N)
!
!---      Mechanical porosity at iterate level k  --
!
          PORDMCHKX = POR(1,N) - PROP_GM(3,N)*DEVKX + OONMODX*DPKX
!
!---      Volumetric strain, at iterate level k+1  ---
!
          EPSV_GM(3,N) = EPS_GM(1,N) + EPS_GM(2,N) + EPS_GM(3,N)
!
!---      Strain differential, at iterate level k+1  ---
!
          DEVK1X = EPSV_GM(3,N)-EPSV_GM(2,N)
!
!---      Pressure differential, at iterate level k+1  ---
!
          DPK1X = P_GM(3,N)-P_GM(2,N)
!
!---      Mechanical porosity at iterate level k+1  --
!
          PORDMCHK1X = PORDMCHKX + OONMODX*DPK1X + PROP_GM(3,N)*DEVK1X
        ENDIF
!
!---    Find maximum local residual  --
!
        RSDX = ABS((PORDFWK1X-PORDMCHK1X)/PORDMCHK1X)
!        IF( N.EQ.1 ) PRINT *,'ND(N) = ',ND(N),'PORDFWK1X = ',PORDFWK1X,
!     &    'PORDMCHK1X = ',PORDMCHK1X,'RSDX = ',RSDX,
!     &    'RSD_GML = ',RSD_GML,'ID = ',ID
        IF( RSDX.GT.RSD_GML ) THEN
          NSD_GML = ND(N)
          RSD_GML = RSDX
        ENDIF
      ENDDO
!
!---  Maximum global residuals  ---
!
      CALL MPI_ALLREDUCE( RSD_GML,RSD_GM,1,MPI_REAL8,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
!
!---  Identify processor with maximum residual  ---
!
      IDLX = -1
      IF( ABS((RSD_GML-RSD_GM)/EPSL).LT.EPSL ) IDLX = ID
      IF( ID.EQ.IDLX ) THEN
        NSD_GM = NSD_GML
      ELSE
        NSD_GM = 0
      ENDIF
      CALL MPI_ALLREDUCE( IDLX,IDX,1,MPI_INTEGER,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
      CALL MPI_BCAST( NSD_GM,1,MPI_INTEGER,IDX,MPI_COMM_WORLD,
     &   IERR )
!      PRINT *,'RSD_GM = ',RSD_GM,'NSD_GM = ',NSD_GM,'ID = ',ID
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDL_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SSM( EX,PX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Compute the 6 x 6 components of the stress-strain matrix
!
!     EX(6,6) - stress-strain matrix
!     PX(1) - Young's modulus, Pa
!     PX(2) - Poisson's ratio
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 November 2016 (Ensisheim meteorite).
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
      REAL*8 PX(2),EX(6,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SSM'
      VARX = PX(1)/((1.D+0+PX(2))*(1.D+0-2.D+0*PX(2)))
      EX(1,1) = VARX*(1.D+0-PX(2))
      EX(1,2) = VARX*PX(2)
      EX(1,3) = VARX*PX(2)
      EX(1,4) = 0.D+0
      EX(1,5) = 0.D+0
      EX(1,6) = 0.D+0
      EX(2,1) = VARX*PX(2)
      EX(2,2) = VARX*(1.D+0-PX(2))
      EX(2,3) = VARX*PX(2)
      EX(2,4) = 0.D+0
      EX(2,5) = 0.D+0
      EX(2,6) = 0.D+0
      EX(3,1) = VARX*PX(2)
      EX(3,2) = VARX*PX(2)
      EX(3,3) = VARX*(1.D+0-PX(2))
      EX(3,4) = 0.D+0
      EX(3,5) = 0.D+0
      EX(3,6) = 0.D+0
      EX(4,1) = 0.D+0
      EX(4,2) = 0.D+0
      EX(4,3) = 0.D+0
      EX(4,4) = VARX*(1.D+0-2.D+0*PX(2))/2.D+0
      EX(4,5) = 0.D+0
      EX(4,6) = 0.D+0
      EX(5,1) = 0.D+0
      EX(5,2) = 0.D+0
      EX(5,3) = 0.D+0
      EX(5,4) = 0.D+0
      EX(5,5) = VARX*(1.D+0-2.D+0*PX(2))/2.D+0
      EX(5,6) = 0.D+0
      EX(6,1) = 0.D+0
      EX(6,2) = 0.D+0
      EX(6,3) = 0.D+0
      EX(6,4) = 0.D+0
      EX(6,5) = 0.D+0
      EX(6,6) = VARX*(1.D+0-2.D+0*PX(2))/2.D+0
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SSM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STATIC_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Static porothermoelastic geomechanics.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 February 2017 (The erroneous word
!     "dord" is discovered in the Webster's New International 
!     Dictionary, Second Edition, prompting an investigation in 1939).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE GEO_MECH







!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP

!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)





!
!
!  Include file for Fortran use of the type(tMat) package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
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




      integer(kind=selected_int_kind(5)) :: IERR

!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STATIC_GM'
!
!---  Zero geomechanical Jacobian matrix  ---
!
      CALL JCBZ_GM
!
!---  Load geomechanical Jacobian matrix  ---
!
      CALL JCBL_GM
!
!---  Modify the geomechanical Jacobian matrix and problem vector
!     for boundary conditions  ---
!
      CALL BCJ_GM
!
!---  Set values of the Jacobian matrix  ---
!
      CALL JCB_SV_GM
!
!---  Linear equation solver  ---
!






      INDX = 2
      CALL STOMP_PETSC_SOLVE( G_KSP,G_MAT,G_RHS_VEC,G_SOL_VEC,
     &  G_SOL_VEC_S,G_SCATTER,INDX )

!
!---  Displacements at FE nodes  ---
!
      CALL DISP_GM






!
!---  Strains at FE centroids  ---
!
      CALL STRAIN_GM
!
!---  Stresses at FE centroids  ---
!
      CALL STRESS_GM
!
!---  Update stress and strain tensor on ghost cells  ---
!
      CALL UPDT_GC_GM
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of STATIC_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STRAIN_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Strain tensor at finite-element centroid.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 31 January 2017 (3M begins marketing
!     Scotch Tape in 1930).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE GEO_MECH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 BX(6,24),UX(24)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STRAIN_GM'
!
!---  Finite-element centroid in the natural
!     coordinate system Xi, Eta, and Mu ---
!
      PXIX = 0.D+0
      PETAX = 0.D+0
      PMUX = 0.D+0
!
!---  Loop over the finite elements (i.e., grid cells) ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---    Loop over the FE nodes ---
!
        DO M0 = 1,8
          M1 = (M0-1)*3
          NFEN = ND_GM(M0,N)
          UX(M1+1) = U_GM(2,NFEN)
          UX(M1+2) = V_GM(2,NFEN)
          UX(M1+3) = W_GM(2,NFEN)
        ENDDO
        CALL ESDM( PXIX,PETAX,PMUX,XE(1,N),YE(1,N),ZE(1,N),BX,DETJX )
        DO M2 = 1,6
          EPS_GM(M2,N) = 0.D+0
        ENDDO
        DO M3 = 1,24
          EPS_GM(1,N) = EPS_GM(1,N) + BX(1,M3)*UX(M3)
          EPS_GM(2,N) = EPS_GM(2,N) + BX(2,M3)*UX(M3)
          EPS_GM(3,N) = EPS_GM(3,N) + BX(3,M3)*UX(M3)
          EPS_GM(4,N) = EPS_GM(4,N) + BX(4,M3)*UX(M3)
          EPS_GM(5,N) = EPS_GM(5,N) + BX(5,M3)*UX(M3)
          EPS_GM(6,N) = EPS_GM(6,N) + BX(6,M3)*UX(M3)
        ENDDO
!
!---    Change sign to make stress/strain to be positive in
!       compression  ---
!
        EPS_GM(1,N) = -EPS_GM(1,N)
        EPS_GM(2,N) = -EPS_GM(2,N)
        EPS_GM(3,N) = -EPS_GM(3,N)
        EPS_GM(4,N) = -5.D-1*EPS_GM(4,N)
        EPS_GM(5,N) = -5.D-1*EPS_GM(5,N)
        EPS_GM(6,N) = -5.D-1*EPS_GM(6,N)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of STRAIN_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STRESS_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Stress tensor at finite-element centroid.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 31 January 2017 (3M begins marketing
!     Scotch Tape in 1930).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE GEO_MECH
      USE FDVH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 EX(6,6),FX(6)
      REAL*8 PROP_GMX(2)
!
!----------------------Data Statements---------------------------------!
!
      DATA FX / 1.D+0,1.D+0,1.D+0,2.D+0,2.D+0,2.D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STRESS_GM'
!
!---  Loop over the finite elements (i.e., grid cells) ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---    Hydrate C-Factor composite model  ---
!
        IF( IHCM_GM(N).EQ.1 ) THEN
          PROP_GMX(1) = PROP_GM(1,N) + 
     &      PROP_GM(5,N)*SH(2,N)*PROP_GM(6,N)
        ELSE
          PROP_GMX(1) = PROP_GM(1,N)
        ENDIF
        PROP_GMX(2) = PROP_GM(2,N)
!
!---    Stress-strain matrix  ---
!
        CALL SSM( EX,PROP_GMX )
        DO M0 = 1,6
          SIG_GM(M0,N) = 0.D+0
        ENDDO
        DO M1 = 1,6
          DO M2 = 1,6
            SIG_GM(M1,N) = SIG_GM(M1,N) + EX(M1,M2)*FX(M2)*EPS_GM(M2,N)
!            IF( ID.EQ.0 .AND. N.EQ.1 ) PRINT *,'M1 = ',M1,'M2 = ',M2,
!     &        'EX(M1,M2) = ',EX(M1,M2),'FX(M2) = ',FX(M2),
!     &        'EPS_GM(M2,N) = ',EPS_GM(M2,N)
          ENDDO
        ENDDO
!!
!!---    Bypass restart check for geomechanics options  ---
!!
!        ISLC50X = ABS(ISLC(50))
!!
!!---    Bulk modulus times 3  ---
!!
!        BLK3X = PROP_GM(1,N)/(1.0D+0-2.0D+0*PROP_GM(2,N))
!!
!!---    Poroelasticity ---
!!
!        IF( MOD(ISLC50X,10).NE.2 .AND. MOD(ISLC50X,10).NE.4 ) THEN
!          PX = MAX( PG(2,N),PL(2,N),PN(2,N) ) + PATM
!          SIG_GM(1,N) = SIG_GM(1,N) + PROP_GM(3,N)*(PX-PCMP(N))
!          SIG_GM(2,N) = SIG_GM(2,N) + PROP_GM(3,N)*(PX-PCMP(N))
!          SIG_GM(3,N) = SIG_GM(3,N) + PROP_GM(3,N)*(PX-PCMP(N))
!        ENDIF
!!
!!---  Thermoelasticity  ---
!!
!        IF( MOD(ISLC50X,10).NE.3 .AND. MOD(ISLC50X,10).NE.4 ) THEN
!          SIG_GM(1,N) = SIG_GM(1,N) + PROP_GM(4,N)*BLK3X*
!     &      (T(2,N)-TCMP(N))
!          SIG_GM(2,N) = SIG_GM(2,N) + PROP_GM(4,N)*BLK3X*
!     &      (T(2,N)-TCMP(N))
!          SIG_GM(3,N) = SIG_GM(3,N) + PROP_GM(4,N)*BLK3X*
!     &      (T(2,N)-TCMP(N))
!        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of STRESS_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_GC_GM
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Update the stress and strain tensor on ghost cells
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 November 2022.
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
      USE GEO_MECH
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_GC_GM'
      NPVX = 12
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
          SBFB(NCS+1) = EPS_GM(1,NLSGC(M+MCS))
          SBFB(NCS+2) = EPS_GM(2,NLSGC(M+MCS))
          SBFB(NCS+3) = EPS_GM(3,NLSGC(M+MCS))
          SBFB(NCS+4) = EPS_GM(4,NLSGC(M+MCS))
          SBFB(NCS+5) = EPS_GM(5,NLSGC(M+MCS))
          SBFB(NCS+6) = EPS_GM(6,NLSGC(M+MCS))
          SBFB(NCS+7) = SIG_GM(1,NLSGC(M+MCS))
          SBFB(NCS+8) = SIG_GM(2,NLSGC(M+MCS))
          SBFB(NCS+9) = SIG_GM(3,NLSGC(M+MCS))
          SBFB(NCS+10) = SIG_GM(4,NLSGC(M+MCS))
          SBFB(NCS+11) = SIG_GM(5,NLSGC(M+MCS))
          SBFB(NCS+12) = SIG_GM(6,NLSGC(M+MCS))
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
        NRCVX = NCGC(6,ID+1)*NPVX
        IDSNDX = NPGC(6,ID+1) - 1
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
          EPS_GM(1,NLRGC(M+MCR)) = RBFB(NCR+1)
          EPS_GM(2,NLRGC(M+MCR)) = RBFB(NCR+2)
          EPS_GM(3,NLRGC(M+MCR)) = RBFB(NCR+3)
          EPS_GM(4,NLRGC(M+MCR)) = RBFB(NCR+4)
          EPS_GM(5,NLRGC(M+MCR)) = RBFB(NCR+5)
          EPS_GM(6,NLRGC(M+MCR)) = RBFB(NCR+6)
          SIG_GM(1,NLRGC(M+MCR)) = RBFB(NCR+7)
          SIG_GM(2,NLRGC(M+MCR)) = RBFB(NCR+8)
          SIG_GM(3,NLRGC(M+MCR)) = RBFB(NCR+9)
          SIG_GM(4,NLRGC(M+MCR)) = RBFB(NCR+10)
          SIG_GM(5,NLRGC(M+MCR)) = RBFB(NCR+11)
          SIG_GM(6,NLRGC(M+MCR)) = RBFB(NCR+12)
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
          SBFS(NCS+1) = EPS_GM(1,NLSGC(M+MCS))
          SBFS(NCS+2) = EPS_GM(2,NLSGC(M+MCS))
          SBFS(NCS+3) = EPS_GM(3,NLSGC(M+MCS))
          SBFS(NCS+4) = EPS_GM(4,NLSGC(M+MCS))
          SBFS(NCS+5) = EPS_GM(5,NLSGC(M+MCS))
          SBFS(NCS+6) = EPS_GM(6,NLSGC(M+MCS))
          SBFS(NCS+7) = SIG_GM(1,NLSGC(M+MCS))
          SBFS(NCS+8) = SIG_GM(2,NLSGC(M+MCS))
          SBFS(NCS+9) = SIG_GM(3,NLSGC(M+MCS))
          SBFS(NCS+10) = SIG_GM(4,NLSGC(M+MCS))
          SBFS(NCS+11) = SIG_GM(5,NLSGC(M+MCS))
          SBFS(NCS+12) = SIG_GM(6,NLSGC(M+MCS))
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
        NRCVX = NCGC(5,ID+1)*NPVX
        IDSNDX = NPGC(5,ID+1) - 1
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
          EPS_GM(1,NLRGC(M+MCR)) = RBFS(NCR+1)
          EPS_GM(2,NLRGC(M+MCR)) = RBFS(NCR+2)
          EPS_GM(3,NLRGC(M+MCR)) = RBFS(NCR+3)
          EPS_GM(4,NLRGC(M+MCR)) = RBFS(NCR+4)
          EPS_GM(5,NLRGC(M+MCR)) = RBFS(NCR+5)
          EPS_GM(6,NLRGC(M+MCR)) = RBFS(NCR+6)
          SIG_GM(1,NLRGC(M+MCR)) = RBFS(NCR+7)
          SIG_GM(2,NLRGC(M+MCR)) = RBFS(NCR+8)
          SIG_GM(3,NLRGC(M+MCR)) = RBFS(NCR+9)
          SIG_GM(4,NLRGC(M+MCR)) = RBFS(NCR+10)
          SIG_GM(5,NLRGC(M+MCR)) = RBFS(NCR+11)
          SIG_GM(6,NLRGC(M+MCR)) = RBFS(NCR+12)
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
          SBFW(NCS+1) = EPS_GM(1,NLSGC(M+MCS))
          SBFW(NCS+2) = EPS_GM(2,NLSGC(M+MCS))
          SBFW(NCS+3) = EPS_GM(3,NLSGC(M+MCS))
          SBFW(NCS+4) = EPS_GM(4,NLSGC(M+MCS))
          SBFW(NCS+5) = EPS_GM(5,NLSGC(M+MCS))
          SBFW(NCS+6) = EPS_GM(6,NLSGC(M+MCS))
          SBFW(NCS+7) = SIG_GM(1,NLSGC(M+MCS))
          SBFW(NCS+8) = SIG_GM(2,NLSGC(M+MCS))
          SBFW(NCS+9) = SIG_GM(3,NLSGC(M+MCS))
          SBFW(NCS+10) = SIG_GM(4,NLSGC(M+MCS))
          SBFW(NCS+11) = SIG_GM(5,NLSGC(M+MCS))
          SBFW(NCS+12) = SIG_GM(6,NLSGC(M+MCS))
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
        NRCVX = NCGC(4,ID+1)*NPVX
        IDSNDX = NPGC(4,ID+1) - 1
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
          EPS_GM(1,NLRGC(M+MCR)) = RBFW(NCR+1)
          EPS_GM(2,NLRGC(M+MCR)) = RBFW(NCR+2)
          EPS_GM(3,NLRGC(M+MCR)) = RBFW(NCR+3)
          EPS_GM(4,NLRGC(M+MCR)) = RBFW(NCR+4)
          EPS_GM(5,NLRGC(M+MCR)) = RBFW(NCR+5)
          EPS_GM(6,NLRGC(M+MCR)) = RBFW(NCR+6)
          SIG_GM(1,NLRGC(M+MCR)) = RBFW(NCR+7)
          SIG_GM(2,NLRGC(M+MCR)) = RBFW(NCR+8)
          SIG_GM(3,NLRGC(M+MCR)) = RBFW(NCR+9)
          SIG_GM(4,NLRGC(M+MCR)) = RBFW(NCR+10)
          SIG_GM(5,NLRGC(M+MCR)) = RBFW(NCR+11)
          SIG_GM(6,NLRGC(M+MCR)) = RBFW(NCR+12)
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
          SBFE(NCS+1) = EPS_GM(1,NLSGC(M+MCS))
          SBFE(NCS+2) = EPS_GM(2,NLSGC(M+MCS))
          SBFE(NCS+3) = EPS_GM(3,NLSGC(M+MCS))
          SBFE(NCS+4) = EPS_GM(4,NLSGC(M+MCS))
          SBFE(NCS+5) = EPS_GM(5,NLSGC(M+MCS))
          SBFE(NCS+6) = EPS_GM(6,NLSGC(M+MCS))
          SBFE(NCS+7) = SIG_GM(1,NLSGC(M+MCS))
          SBFE(NCS+8) = SIG_GM(2,NLSGC(M+MCS))
          SBFE(NCS+9) = SIG_GM(3,NLSGC(M+MCS))
          SBFE(NCS+10) = SIG_GM(4,NLSGC(M+MCS))
          SBFE(NCS+11) = SIG_GM(5,NLSGC(M+MCS))
          SBFE(NCS+12) = SIG_GM(6,NLSGC(M+MCS))
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
        NRCVX = NCGC(3,ID+1)*NPVX
        IDSNDX = NPGC(3,ID+1) - 1
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
          EPS_GM(1,NLRGC(M+MCR)) = RBFE(NCR+1)
          EPS_GM(2,NLRGC(M+MCR)) = RBFE(NCR+2)
          EPS_GM(3,NLRGC(M+MCR)) = RBFE(NCR+3)
          EPS_GM(4,NLRGC(M+MCR)) = RBFE(NCR+4)
          EPS_GM(5,NLRGC(M+MCR)) = RBFE(NCR+5)
          EPS_GM(6,NLRGC(M+MCR)) = RBFE(NCR+6)
          SIG_GM(1,NLRGC(M+MCR)) = RBFE(NCR+7)
          SIG_GM(2,NLRGC(M+MCR)) = RBFE(NCR+8)
          SIG_GM(3,NLRGC(M+MCR)) = RBFE(NCR+9)
          SIG_GM(4,NLRGC(M+MCR)) = RBFE(NCR+10)
          SIG_GM(5,NLRGC(M+MCR)) = RBFE(NCR+11)
          SIG_GM(6,NLRGC(M+MCR)) = RBFE(NCR+12)
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
          SBFN(NCS+1) = EPS_GM(1,NLSGC(M+MCS))
          SBFN(NCS+2) = EPS_GM(2,NLSGC(M+MCS))
          SBFN(NCS+3) = EPS_GM(3,NLSGC(M+MCS))
          SBFN(NCS+4) = EPS_GM(4,NLSGC(M+MCS))
          SBFN(NCS+5) = EPS_GM(5,NLSGC(M+MCS))
          SBFN(NCS+6) = EPS_GM(6,NLSGC(M+MCS))
          SBFN(NCS+7) = SIG_GM(1,NLSGC(M+MCS))
          SBFN(NCS+8) = SIG_GM(2,NLSGC(M+MCS))
          SBFN(NCS+9) = SIG_GM(3,NLSGC(M+MCS))
          SBFN(NCS+10) = SIG_GM(4,NLSGC(M+MCS))
          SBFN(NCS+11) = SIG_GM(5,NLSGC(M+MCS))
          SBFN(NCS+12) = SIG_GM(6,NLSGC(M+MCS))
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
        NRCVX = NCGC(2,ID+1)*NPVX
        IDSNDX = NPGC(2,ID+1) - 1
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
          EPS_GM(1,NLRGC(M+MCR)) = RBFN(NCR+1)
          EPS_GM(2,NLRGC(M+MCR)) = RBFN(NCR+2)
          EPS_GM(3,NLRGC(M+MCR)) = RBFN(NCR+3)
          EPS_GM(4,NLRGC(M+MCR)) = RBFN(NCR+4)
          EPS_GM(5,NLRGC(M+MCR)) = RBFN(NCR+5)
          EPS_GM(6,NLRGC(M+MCR)) = RBFN(NCR+6)
          SIG_GM(1,NLRGC(M+MCR)) = RBFN(NCR+7)
          SIG_GM(2,NLRGC(M+MCR)) = RBFN(NCR+8)
          SIG_GM(3,NLRGC(M+MCR)) = RBFN(NCR+9)
          SIG_GM(4,NLRGC(M+MCR)) = RBFN(NCR+10)
          SIG_GM(5,NLRGC(M+MCR)) = RBFN(NCR+11)
          SIG_GM(6,NLRGC(M+MCR)) = RBFN(NCR+12)
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
          SBFT(NCS+1) = EPS_GM(1,NLSGC(M+MCS))
          SBFT(NCS+2) = EPS_GM(2,NLSGC(M+MCS))
          SBFT(NCS+3) = EPS_GM(3,NLSGC(M+MCS))
          SBFT(NCS+4) = EPS_GM(4,NLSGC(M+MCS))
          SBFT(NCS+5) = EPS_GM(5,NLSGC(M+MCS))
          SBFT(NCS+6) = EPS_GM(6,NLSGC(M+MCS))
          SBFT(NCS+7) = SIG_GM(1,NLSGC(M+MCS))
          SBFT(NCS+8) = SIG_GM(2,NLSGC(M+MCS))
          SBFT(NCS+9) = SIG_GM(3,NLSGC(M+MCS))
          SBFT(NCS+10) = SIG_GM(4,NLSGC(M+MCS))
          SBFT(NCS+11) = SIG_GM(5,NLSGC(M+MCS))
          SBFT(NCS+12) = SIG_GM(6,NLSGC(M+MCS))
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
        NRCVX = NCGC(1,ID+1)*NPVX
        IDSNDX = NPGC(1,ID+1) - 1
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
          EPS_GM(1,NLRGC(M+MCR)) = RBFT(NCR+1)
          EPS_GM(2,NLRGC(M+MCR)) = RBFT(NCR+2)
          EPS_GM(3,NLRGC(M+MCR)) = RBFT(NCR+3)
          EPS_GM(4,NLRGC(M+MCR)) = RBFT(NCR+4)
          EPS_GM(5,NLRGC(M+MCR)) = RBFT(NCR+5)
          EPS_GM(6,NLRGC(M+MCR)) = RBFT(NCR+6)
          SIG_GM(1,NLRGC(M+MCR)) = RBFT(NCR+7)
          SIG_GM(2,NLRGC(M+MCR)) = RBFT(NCR+8)
          SIG_GM(3,NLRGC(M+MCR)) = RBFT(NCR+9)
          SIG_GM(4,NLRGC(M+MCR)) = RBFT(NCR+10)
          SIG_GM(5,NLRGC(M+MCR)) = RBFT(NCR+11)
          SIG_GM(6,NLRGC(M+MCR)) = RBFT(NCR+12)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_GC_GM group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VOLSS_GM( INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Reference volumetric strain at finite-element centroid.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 21 February 2017 (Initial issue of 
!     the Cherokee Phoenix in 1828).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE GEO_MECH
      USE FDVH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 PROP_GMX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VOLSS_GM'
!
!---  Loop over the finite elements (i.e., grid cells) ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    Reference volumetric strains  ---
!
        IF( INDX.EQ.0 ) THEN
          EPSV_CMP(N) = EPS_GM(1,N) + EPS_GM(2,N) + EPS_GM(3,N)
          SIGV_CMP(N) = (SIG_GM(1,N) + SIG_GM(2,N) + SIG_GM(3,N))/3.D+0
!
!---    k or k+1 level volumetric strains  ---
!
        ELSE
          EPSV_GM(INDX,N) = EPS_GM(1,N) + EPS_GM(2,N) + EPS_GM(3,N)
          SIGV_GM(INDX,N) = (SIG_GM(1,N)+SIG_GM(2,N)+SIG_GM(3,N))/3.D+0
!
!---      Hydrate C-Factor composite model  ---
!
          IF( IHCM_GM(N).EQ.1 ) THEN
            PROP_GMX(1) = PROP_GM(1,N) + 
     &        PROP_GM(5,N)*SH(2,N)*PROP_GM(6,N)
          ELSE
            PROP_GMX(1) = PROP_GM(1,N)
          ENDIF
          PROP_GMX(2) = PROP_GM(2,N)
!
!---      Bulk modulus  ---
!
          BLKX = PROP_GMX(1)/(3.D+0*(1.0D+0-2.0D+0*PROP_GMX(2)))
          EPSV_CMP(N) = SIGV_CMP(N)/BLKX
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VOLSS_GM group  ---
!
      RETURN
      END




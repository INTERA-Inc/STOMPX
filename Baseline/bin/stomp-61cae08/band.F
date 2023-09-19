!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BAND( ISLV,MU,ML,INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Banded matrix linear system solver.
!
!     LINPACK User's Guide. J.J. Dongarra, C.B. Moler, J.R. Bunch, and
!     G.W. Stewart. SIAM, Philadelphia, 1979.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, March, 1993.
!     Last Modified by MD White, Battelle, PNL, April 14, 1994.
!     band.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_CL
      USE SOLTN
      USE LEAK_WELL
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE GEO_MECH
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
      CHARACTER*64 FDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BAND'
      IERR = 0
!
!---  Number of unknowns equals number of active nodes + 
!     number of fracture triangles + number of borehole nodes for
!     transport equations ---
!
      IF( INDX.EQ.1 ) THEN
        N = NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC + NBN_BH - 
     &    NXP_BH
!
!---  Geomechanics - number of unknowns equals number of active
!     finite element nodes x 3
!
      ELSEIF( INDX.EQ.2 ) THEN
        N = 3*NFEN_GM
!
!---  Number of unknowns equals number of number of coupled equations x
!     (active nodes + fracture triangles) + number of coupled wells  ---
!
      ELSEIF( N_CW.GT.0 ) THEN
        N = ISLV*(NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC) + N_CW
!
!---  Number of unknowns equals number of number of coupled equations
!     x (number of active nodes + number of wells +
!     number of surface nodes)
!
      ELSEIF( LSPILL.GT.0 ) THEN
        N = ISLV*(NFBN - NRFN - NXP + NWLN + NWN_LW + IJFLD)
!
!---  Number of unknowns equals number of number of coupled equations
!     x (number of active nodes + number of wells + fracture triangles
!     + number of borehole nodes)
!
      ELSE
        N = ISLV*(NFBN - NRFN - NXP + NWLN + NWN_LW + NT_FRC - NXP_FRC
     &    + NBN_BH - NXP_BH)
      ENDIF
!
!---  Dual-porosity option  ---
!
      IF( ISLC(11).EQ.1 .AND. INDX.EQ.1 ) 
     &  N = N + ISLV*(NFBN-NRFN-NXP+NWN_LW)
      MDX = ML + MU + 1
!
!---  Solute transport ---
!
      IF( INDX.EQ.1 ) THEN
!
!---    Output transport problem vector  ---
!
        IF( ISLC(34).EQ.11 ) THEN
          FDUM = 'rhs_vec_c.dat'
          CALL BAND_OUTPUT_VECTOR( N,FDUM )
          STOP
        ENDIF
!
!---    Output transport matrix  ---
!
        IF( ISLC(34).EQ.12 ) THEN
          FDUM = 'mat_c.dat'
          CALL BAND_OUTPUT_MATRIX( MDT,MLT,MUT,N,FDUM )
          STOP
        ENDIF
!
!---  Geomechancis  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
!
!---    Output geomechanics problem vector  ---
!
        IF( ISLC(34).EQ.21 ) THEN
          FDUM = 'rhs_vec_g.dat'
          CALL BAND_OUTPUT_VECTOR( N,FDUM )
          STOP
        ENDIF
!
!---    Output geomechanics matrix  ---
!
        IF( ISLC(34).EQ.22 ) THEN
          FDUM = 'mat_g.dat'
          CALL BAND_OUTPUT_MATRIX( MDX,ML,MU,N,FDUM )
          STOP
        ENDIF
!
!---  Coupled flow  ---
!
      ELSE
!
!---    Output coupled-flow problem vector  ---
!
        IF( ISLC(34).EQ.1 ) THEN
          FDUM = 'rhs_vec.dat'
          CALL BAND_OUTPUT_VECTOR( N,FDUM )
          STOP
        ENDIF
!
!---    Output coupled-flow matrix  ---
!
        IF( ISLC(34).EQ.2 ) THEN
          FDUM = 'mat.dat'
          CALL BAND_OUTPUT_MATRIX( MDC,MLC,MUC,N,FDUM )
          STOP
        ENDIF
      ENDIF
!
!---  Identify zero valued solution vector elements  ---
!
      DO 100 M = 1,N
        JLU(M) = 0
        IF( ABS(BLU(M)).LE.SMALL ) JLU(M) = 1
        IF( ABS(ALU(MDX,M))/EPSL.LE.EPSL ) THEN
          IERR = M
        ENDIF
  100 CONTINUE
!
!---  Factorization  ---
!
      IF( IERR.EQ.0 ) THEN
        CALL DBGBFA( N,MU,ML,IERR )
      ENDIF
!
!---  Singular matrix error message  ---
!
      IF( IERR.NE.0 ) THEN
        NODE = 0
!
!---   Check coupled well equations  ---
!
       DO 230 NCW = 1,N_CW
         IF( JM_CW(NCW).EQ.IERR ) THEN
           WRITE(ISC,'(2X,A)') '---  Singular Matrix  ---'
           WRITE(ISC,'(2X,A,I6)') 'Jacobian Matrix Index = ',IERR
           WRITE(ISC,'(2X,A,I6)') 'Coupled Well Number   = ',NCW
           WRITE(IWR,'(2X,A)') '---  Singular Matrix  ---'
           WRITE(IWR,'(2X,A,I6)') 'Jacobian Matrix Index = ',IERR
           WRITE(IWR,'(2X,A,I6)') 'Coupled Well Number   = ',NCW
           ICNV = 1
           ISUB_LOG = ISUB_LOG-1
           RETURN
         ENDIF
  230  CONTINUE
!
!---   Check field node equations  ---
!
        DO 250 N = 1,NFBN
          IF( IXP(N).EQ.0 ) GOTO 250
          NMD = IXP(N)
          DO 240 M = 1,ISLV
            IF( (NMD-1)*ISLV+M.EQ.IERR ) THEN
              NODE = N
              GOTO 260
            ENDIF
  240   CONTINUE
  250   CONTINUE
  260   CONTINUE
        WRITE(ISC,'(2X,A)') '---  Singular Matrix  ---'
        WRITE(ISC,'(2X,A,I6)') 'Jacobian Matrix Index = ',IERR
        WRITE(ISC,'(2X,A,I6)') 'Node Number           = ',NODE
        WRITE(IWR,'(2X,A)') '---  Singular Matrix  ---'
        WRITE(IWR,'(2X,A,I6)') 'Jacobian Matrix Index = ',IERR
        WRITE(IWR,'(2X,A,I6)') 'Node Number           = ',NODE
        ICNV = 1
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Solve linear equations  ---
!
      CALL DBGBSL( N,MU,ML )
!
!---  Solute transport ---
!
      IF( INDX.EQ.1 ) THEN
!
!---    Output transport solution vector  ---
!
        IF( ISLC(34).EQ.13 ) THEN
          FDUM = 'fsol_vec_c.dat'
          CALL BAND_OUTPUT_VECTOR( N,FDUM )
          STOP
        ENDIF
!
!---  Geomechancis  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
!
!---    Output geomechanics solution vector  ---
!
        IF( ISLC(34).EQ.23 ) THEN
          FDUM = 'fsol_vec_g.dat'
          CALL BAND_OUTPUT_VECTOR( N,FDUM )
          STOP
        ENDIF
!
!---  Fracture flow  ---
!
      ELSEIF( INDX.EQ.3 ) THEN
!
!---    Output geomechanics solution vector  ---
!
        IF( ISLC(34).EQ.33 ) THEN
          FDUM = 'fsol_vec_frc.dat'
          CALL BAND_OUTPUT_VECTOR( N,FDUM )
          STOP
        ENDIF
!
!---  Coupled flow  ---
!
      ELSE
!
!---    Output coupled-flow solution vector  ---
!
        IF( ISLC(34).EQ.3 ) THEN
          FDUM = 'fsol_vec.dat'
          CALL BAND_OUTPUT_VECTOR( N,FDUM )
          STOP
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BAND group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BAND_OUTPUT_VECTOR( NUK,FDUM )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!    Output the problem or solution vector
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 03 May 2016
!
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
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,FDUM
      CHARACTER*6 FORM1
      CHARACTER*16 FORM2
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1)'/
      DATA FORM2 /'(I1,2X,1pe27.20)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BAND_OUTPUT_VECTOR'
      NCH = LEN_TRIM(FDUM)
      OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='UNKNOWN',FORM='FORMATTED')
      CLOSE(UNIT=26,STATUS='DELETE')
      OPEN(UNIT=26, FILE=FDUM(1:NCH), STATUS='NEW', FORM='FORMATTED')
      WRITE(26,'(A)') '%%MatrixMarket vector coordinate real general'
      ICX = ICOUNT(NUK)
      WRITE(FORM1(3:3),'(I1)') ICX
      WRITE(26,FORM1) NUK
      DO N = 1,NUK
        ICX = ICOUNT(N)
        WRITE(FORM2(3:3),'(I1)') ICX
        WRITE(ADUM,FORM2) N,BLU(N)
        CALL LCASE( ADUM )
        NCH = LEN_TRIM(ADUM)
        WRITE(26,'(A)') ADUM(1:NCH)
      ENDDO
      CLOSE(UNIT=26)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BAND_OUTPUT_VECTOR group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BAND_OUTPUT_MATRIX( MDX,MLX,MUX,NUK,FDUM )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!    Output the Jacobian or transport matrix
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 December 2018
!
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
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,FDUM
      CHARACTER*16 FORM1
      CHARACTER*22 FORM2
      CHARACTER*21 FORM3
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1,1X,I1,1X,I1)'/
      DATA FORM2 /'(I1,1X,I1,2X,1PE27.20)'/
      DATA FORM3 /'(I1,1X,I1,2X,1PE12.5)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BAND_OUTPUT_MATRIX'
      NCH = LEN_TRIM(FDUM)
      OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='UNKNOWN',FORM='FORMATTED')
      CLOSE(UNIT=26,STATUS='DELETE')
      OPEN(UNIT=26, FILE=FDUM(1:NCH), STATUS='NEW', FORM='FORMATTED')
      WRITE(26,'(A)') '%%MatrixMarket matrix coordinate real general'
      ICX = ICOUNT(NUK)
      WRITE(FORM1(3:3),'(I1)') ICX
      WRITE(FORM1(9:9),'(I1)') ICX
      NNZ = NUK*(MUX+MLX+1)
      ICX = ICOUNT(NNZ)
      WRITE(FORM1(15:15),'(I1)') ICX
      WRITE(26,FORM1) NUK,NUK,NNZ
!
!---  Loop over equations  ---
!
      DO N = 1,NUK
!
!---    Loop over entries in equations, between upper and lower
!       half band widths  ---
!
        M1 = MAX( 1,N-MUX )
        M2 = MIN( NUK,N+MLX )
        DO M = M1,M2       
          ICX = ICOUNT(N)
          WRITE(FORM3(3:3),'(I1)') ICX
          ICX = ICOUNT(M)
          WRITE(FORM3(9:9),'(I1)') ICX
          IROW = MDX - M + N
          ICOL = M
          WRITE(ADUM,FORM3) N,M,ALU(IROW,ICOL)
          CALL LCASE( ADUM )
          NCH = LEN_TRIM(ADUM)
          WRITE(26,'(A)') ADUM(1:NCH)
        ENDDO
      ENDDO
      CLOSE(UNIT=26)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BAND_OUTPUT_MATRIX group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DIAG_SC( IEQX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
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
!     Diagonal scaling of the Jacobian matrix and solution vector
!     for the banded matrix storage.
!
!----------------------Authors-----------------------------------------!
!
!     Written MD White, PNNL, 12 July 2002.
!     Last Modified by MD White, PNNL, 12 July 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
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
      SUB_LOG(ISUB_LOG) = '/DIAG_SC'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO K = 1,KFLD
      DO J = 1,JFLD
      DO I = 1,IFLD
        N = ND(I,J,K)
        IF( IXP(N).EQ.0 ) CYCLE
        IF( ILES.EQ.1 ) THEN
!
!---  Node  ---
!
          NMD = IXP(N)
          MP = IM(IEQX,NMD)
          DIAGX = ALU(MDC,MP)
          DO 100 M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAGX
  100     CONTINUE
          BLU(MP) = BLU(MP)/DIAGX
!
!---  Bottom ---
!
          IF( K.NE.1 ) THEN
            NB = N-IJFLD
            IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 210
            NMD = IXP(NB)
            DO 200 M = 1,ISVC
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAGX
  200       CONTINUE
          ENDIF
  210   CONTINUE
!
!---  South ---
!
          IF( J.NE.1 ) THEN
            NS = N-IFLD
            IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 310
            NMD = IXP(NS)
            DO 300 M = 1,ISVC
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAGX
  300       CONTINUE
          ENDIF
  310     CONTINUE
!
!---  West ---
!
          IF( I.NE.1 ) THEN
            NW = N-1
            IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 410
            NMD = IXP(NW)
            DO 400 M = 1,ISVC
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAGX
  400       CONTINUE
          ENDIF
  410     CONTINUE
!
!---  East ---
!
          IF( I.NE.IFLD ) THEN
            NE = N+1
            IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 510
            NMD = IXP(NE)
            DO 500 M = 1,ISVC
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAGX
  500       CONTINUE
          ENDIF
  510     CONTINUE
!
!---  North ---
!
          IF( J.NE.JFLD ) THEN
            NN = N+IFLD
            IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 610
            NMD = IXP(NN)
            DO 600 M = 1,ISVC
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAGX
  600       CONTINUE
          ENDIF
  610     CONTINUE
!
!---  Top ---
!
          IF( K.NE.KFLD ) THEN
            NT = N+IJFLD
            IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 710
            NMD = IXP(NT)
            DO 700 M = 1,ISVC
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAGX
  700       CONTINUE
          ENDIF
  710     CONTINUE
        ENDIF
      ENDDO
      ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DIAG_SC group  ---
!
      RETURN
      END


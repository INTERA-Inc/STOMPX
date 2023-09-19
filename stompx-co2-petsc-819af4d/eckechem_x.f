!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_PTZRCOEF
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
!     Allocate global array memory for Pitzer coefficients
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 February 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE PTZRCOEF
      USE GRID
      USE GLB_PAR
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
      SUB_LOG(ISUB_LOG) = '/ALLOC_PTZRCOEF'
      ALLOCATE( JPA(1:LANI),STAT=ISTAT )
      CHMSG = 'JPA'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( JPC(1:LCAT),STAT=ISTAT )
      CHMSG = 'JPC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( JPN(1:LNEU),STAT=ISTAT )
      CHMSG = 'JPN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( B0(1:LCAT,1:LANI),STAT=ISTAT )
      CHMSG = 'B0'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( B1(1:LCAT,1:LANI),STAT=ISTAT )
      CHMSG = 'B1'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( B2(1:LCAT,1:LANI),STAT=ISTAT )
      CHMSG = 'B2'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CMXX(1:LCAT,1:LANI),STAT=ISTAT )
      CHMSG = 'CMXX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TCC(1:LNCF),STAT=ISTAT )
      CHMSG = 'TCC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TAA(1:LNAF),STAT=ISTAT )
      CHMSG = 'TAA'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PSIC(1:LNCF,1:LANI),STAT=ISTAT )
      CHMSG = 'PSIC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PSIA(1:LNAF,1:LCAT),STAT=ISTAT )
      CHMSG = 'PSIA'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ALAMB(1:LNEU,1:LANI),STAT=ISTAT )
      CHMSG = 'ALAMB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CLAMB(1:LNEU,1:LCAT),STAT=ISTAT )
      CHMSG = 'CLAMB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ELAMB(1:LNEU,1:LNEU),STAT=ISTAT )
      CHMSG = 'ELAMB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( HOLAMB(1:LNNF,1:LANI),STAT=ISTAT )
      CHMSG = 'HOLAMB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( BPPR(1:LCAT,1:LANI),STAT=ISTAT )
      CHMSG = 'BPPR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( BPHI(1:LCAT,1:LANI),STAT=ISTAT )
      CHMSG = 'BPHI'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( BPR(1:LCAT,1:LANI),STAT=ISTAT )
      CHMSG = 'BPR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( BMMX(1:LCAT,1:LANI),STAT=ISTAT )
      CHMSG = 'BMMX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATB0(1:LCAT,1:LANI,1:8),STAT=ISTAT )
      CHMSG = 'ATB0'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATB1(1:LCAT,1:LANI,1:8),STAT=ISTAT )
      CHMSG = 'ATB1'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATB2(1:LCAT,1:LANI,1:8),STAT=ISTAT )
      CHMSG = 'ATB2'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATCMX(1:LCAT,1:LANI,1:8),STAT=ISTAT )
      CHMSG = 'ATCMX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATNLAM(1:LNEU,1:LNEU,1:6),STAT=ISTAT )
      CHMSG = 'ATNLAM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATCLAM(1:LNEU,1:LCAT,1:6),STAT=ISTAT )
      CHMSG = 'ATCLAM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATALAM(1:LNEU,1:LANI,1:6),STAT=ISTAT )
      CHMSG = 'ATALAM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATHLAM(1:LNNF,1:LANI,1:6),STAT=ISTAT )
      CHMSG = 'ATHLAM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATTC(1:LNCF,1:6),STAT=ISTAT )
      CHMSG = 'ATTC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATPC(1:LNCF,1:LANI,1:6),STAT=ISTAT )
      CHMSG = 'ATPC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATTA(1:LNAF,1:6),STAT=ISTAT )
      CHMSG = 'ATTA'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ATPA(1:LNAF,1:LCAT,1:6),STAT=ISTAT )
      CHMSG = 'ATPA'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CTCPH(1:LNCF),STAT=ISTAT )
      CHMSG = 'CTCPH'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CTC(1:LNCF),STAT=ISTAT )
      CHMSG = 'CTC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CTCPR(1:LNCF),STAT=ISTAT )
      CHMSG = 'CTCPR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CTCPPR(1:LNCF),STAT=ISTAT )
      CHMSG = 'CTCPPR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CTAPH(1:LNAF),STAT=ISTAT )
      CHMSG = 'CTAPH'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CTA(1:LNAF),STAT=ISTAT )
      CHMSG = 'CTA'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CTAPR(1:LNAF),STAT=ISTAT )
      CHMSG = 'CTAPR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CTAPPR(1:LNAF),STAT=ISTAT )
      CHMSG = 'CTAPPR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ETH(1:LMCG*LMCG),STAT=ISTAT )
      CHMSG = 'ETH'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ETHP(1:LMCG*LMCG),STAT=ISTAT )
      CHMSG = 'ETHP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ETHP2(1:LMCG*LMCG),STAT=ISTAT )
      CHMSG = 'ETHP2'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IDD_PZ(1:LSPL),STAT=ISTAT )
      CHMSG = 'IDD_PZ'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_PTZRCOEF group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCK_STP
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
!     ECKEChemX
!
!     Back step.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE HYST
      USE GRID
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
      SUB_LOG(ISUB_LOG) = '/BCK_STP'
      NSTEP = NSTEP-1
!
!---  STOMPX-CO2 Operational Mode  ---
!
      IF( IOM.EQ.32 ) THEN
!
!---    Loop over all nodes, skipping inactive nodes  ---
!
        DO N = 1,NFCGC(ID+1)
          IF( IXP(N).EQ.0 ) CYCLE
          T(2,N) = T(1,N)
          PL(2,N) = PL(1,N)
          PG(2,N) = PG(1,N)
          SL(2,N) = SL(1,N)
          SGT(2,N) = SGT(1,N)
          ASLMIN(2,N) = ASLMIN(1,N)
          SG(2,N) = SG(1,N)
          XLA(2,N) = XLA(1,N)
          YLS(2,N) = YLS(1,N)
          NPHAZ(2,N) = NPHAZ(1,N)
          DO NSL = 1,NSOLU
            C(N,NSL) = CO(N,NSL)
          ENDDO
          DO NSP = 1,NSPR
            SP_C(N,NSP) = SP_CO(N,NSP)
          ENDDO
        ENDDO

      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCK_STP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BDOT( ACTVX,SP_CX,DSP_CX,SLX,PORDX,RHOLX,TX,XLWX )
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
!     ECKEChemX
! 
!     B-Dot or extended Debye-Huckel model for activity coefficient.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 ACTVX(LSPL+1,LSPL),SP_CX(LSPL),DSP_CX(LSPL)
      REAL*8 CPIX(LSPL+1)
      REAL*8 ACOF(5),BCOF(5),BDCOF(5)
      REAL*8 SACOF(4),SBCOF(4),SCCOF(4)
!
!----------------------Data Statements---------------------------------!
!
      DATA ACOF / 0.49463D+0,0.000398438D+0,9.9092D-6,-4.36244D-8,
     &  1.09811D-10 /
      DATA BCOF / 0.325325D+0,0.000122253D+0,6.68944D-7,-2.57037D-9,
     &  4.89847D-12 /
      DATA BDCOF / 0.0373904D+0,0.000191209D+0,-2.05222D-6,
     &  1.31092D-8,-3.26031D-11 /
      DATA SACOF / 0.131678D+0,-0.000836829D+0,3.07179D-06,
     &  1.46701D-09 /
      DATA SBCOF / -0.0186731D+0,0.00039022D+0,-2.62611D-06,
     &  4.40918D-09 /
      DATA SCCOF / 0.00288841D+0,-6.70405D-05,5.65666D-07,
     &  -1.34012D-09 /
!      DATA TSMI / -1.D+3 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BDOT'
!
!---  Coefficients as a function of temperature  ---
!
      ACPX = ACOF(5)
      BCPX = BCOF(5)
      BDCPX = BDCOF(5)
      DO I = 4,1,-1
        ACPX = ACPX*TX + ACOF(I)
        BCPX = BCPX*TX + BCOF(I)
        BDCPX = BDCPX*TX + BDCOF(I)
      ENDDO
      ASMX = SACOF(4)
      BSMX = SBCOF(4)
      CSMX = SCCOF(4)
      DO I = 3,1,-1
        ASMX = ASMX*TX + SACOF(I)
        BSMX = BSMX*TX + SBCOF(I)
        CSMX = CSMX*TX + SCCOF(I)
      ENDDO
!
!---  Ionic strength of the aqueous solution  ---
!
      DO M = 1,NSPL+1
        CPIX(M) = 0.D+0
        SUM_MZ = 0.D+0
        SUM_M = 0.D+0
        DO NSP = 1,NSPL
!
!---      Skip neutral species  ---
!
          IF( ABS(SP_L(1,NSP)).LT.EPSL ) CYCLE
          IF( NSP.EQ.(M-1) ) THEN
            CLX = SP_CX(NSP) + DSP_CX(NSP)
          ELSE
            CLX = SP_CX(NSP)
          ENDIF
!
!---      Molarity in mol solute/m^3 aqueous
!         or mol solute/l aqueous  ---
!
          CMX = CLX/(SLX*PORDX)
!
!---      Molality in mol solute/kg water  ---
!
          CMX = CMX/(RHOLX*XLWX)
          SUM_M = SUM_M + CMX
          SUM_MZ = SUM_MZ + CMX*SP_L(1,NSP)
          CPIX(M) = CPIX(M) + CMX*(SP_L(1,NSP)**2)
        ENDDO
!
!---    Correct for electrical neutrality  ---
!
        IF( ABS(SUM_MZ).GT.EPSL )
     &    CPIX(M) = CPIX(M) + (SUM_MZ**2)/SUM_M
        CPIX(M) = CPIX(M)*5.D-1
      ENDDO
!
!---  Activity coefficients for aqueous species  ---
!
      DO M = 1,NSPL+1
        DO NSP = 1,NSPL
          IF( ABS(SP_L(1,NSP)).GT.EPSL ) THEN
!
!---        Activity coefficients for charged species  ---
!
            ACTVX(M,NSP) = -ACPX*(SP_L(1,NSP)**2)*SQRT(CPIX(M))
     &      /(1.D+0 + SP_L(2,NSP)*1.D+10*BCPX*SQRT(CPIX(M)))
     &      + BDCPX*CPIX(M)
            ACTVX(M,NSP) = EXP(TOLN*ACTVX(M,NSP))
!
!---        Activity coefficients for electrically neutral species  ---
!
          ELSE
            ACTVX(M,NSP) = EXP(TOLN*(ASMX*CPIX(M) + BSMX*(CPIX(M)**2) +
     &        CSMX*(CPIX(M)**3)))
          ENDIF
          ACTVX(M,NSP) = MIN( ACTVX(M,NSP),1.D+1 )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BDOT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CECHEM( ACTVX,AJM,BJM,CX,SP_CX,N,NEQ,INDX )
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
!     ECKEChemX
! 
!     Conservation Equation CHEMistry
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AJM(LSPR,LSPR),BJM(LSPR)
      REAL*8 SP_CX(LSPR)
      REAL*8 CX(LEQC+LEQK)
      REAL*8 ACTVX(LSPL+1,LSPL)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CECHEM'
!
!---  Volumetric concentration to aqueous concentration,
!     mol/m^3 -> mol/m^3 aqu  ---
!
      VTOLX = 1.D+0/(SL(2,N)*PORD(2,N))
!
!---  Loop over conservation equation species  ---
!
      NEQX = NEQ - NEQE
      BJM(NEQ) = CX(NEQX)*VTOLX
!
!---  Skip for initial pH  ---
!
      IF( ISPLK(1).GT.1000 ) THEN
        IF( NEQ.EQ.IEQ_S(MOD(ISPLK(1),1000)) .AND.
     &    NSTEP.EQ.0 ) THEN
          AJM(NEQ,NEQ) = 1.D+0
          BJM(NEQ) = 0.D+0
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
      ENDIF
!
!---  Fixed species concentration or activity  ---
!
      DO NSPKX = 1,NSPLK
        IF( ISPLK(14+NSPKX).LT.0 ) THEN
          NSPX = ABS(ISPLK(14+NSPKX))
          IF( NSPX.GT.1000 ) NSPX = NSPX - 1000
          IF( NEQ.EQ.IEQ_S(NSPX) ) THEN
            AJM(NEQ,NEQ) = 1.D+0
            BJM(NEQ) = 0.D+0
            ISUB_LOG = ISUB_LOG-1
            RETURN
          ENDIF
        ENDIF
      ENDDO
      L1: DO M = 1,IEQ_C(1,NEQX)
        NSP = IEQ_C(M+1,NEQX)
        BJM(NEQ) = BJM(NEQ) - EQ_C(M,NEQX)*SP_CX(NSP)*VTOLX
!
!---    Skip for initial pH  ---
!
        IF( ISPLK(1).GT.1000 ) THEN
          IF( NSP.EQ.MOD(ISPLK(1),1000) .AND.
     &      NSTEP.EQ.0 ) CYCLE L1
        ENDIF
!
!---    Skip fixed species concentration ---
!
        DO NSPKX = 1,NSPLK
          IF( ISPLK(14+NSPKX).LT.0 ) THEN
            NSPX = ABS(ISPLK(14+NSPKX))
            IF( NSP.EQ.NSPX) CYCLE L1
          ENDIF
        ENDDO
        AJM(NEQ,IEQ_S(NSP)) = - EQ_C(M,NEQX)*VTOLX
        DO NSPKX = 1,NSPLK
          IF( ISPLK(14+NSPKX).LT.0 ) THEN
            NSPX = ABS(ISPLK(14+NSPKX))
            IF( NSPX .GT. 1000 ) NSPX = NSPX - 1000
            IF( NSP.EQ.NSPX ) AJM(NEQ,IEQ_S(NSP)) =
     &        -1.D+0/ACTVX(1,NSP)*EQ_C(M,NEQX)*VTOLX
          ENDIF
        ENDDO
      ENDDO L1
!
!---  Return residual vector and Jacobian matrix  ---
!
      IF( ABS(INDX).EQ.1 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      BJM(NEQ) = -BJM(NEQ)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CECHEM group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CECHEM_R( ACTVX,AJM,BJM,CX,SP_CX,N,NEQ,INDX )
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
!     ECKEChemX
! 
!     Conservation Equation CHEMistry
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AJM(LSPR,LSPR),BJM(LSPR)
      REAL*8 SP_CX(LSPR)
      REAL*8 CX(LEQC+LEQK)
      REAL*8 ACTVX(LSPL+1,LSPL)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CECHEM_R'
!
!---  Volumetric concentration to aqueous concentration,
!     mol/m^3 -> mol/m^3 aqu  ---
!
      VTOLX = 1.D+0/(SL(2,N)*PORD(2,N))
!
!---  Loop over conservation equation species  ---
!
      NEQX = NEQ - NEQE
      NROW = NEQX
      BJM(NROW) = CX(NEQX)*VTOLX
!
!---  Skip for initial pH  ---
!
      IF( ISPLK(1).GT.1000 ) THEN
        IF( NEQ.EQ.IEQ_S(MOD(ISPLK(1),1000)) .AND.
     &    NSTEP.EQ.0 ) THEN
          AJM(NROW,NROW) = 1.D+0
          BJM(NROW) = 0.D+0
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
      ENDIF
!
!--- fixed species concentration or activity
!
      DO NSPKX = 1,NSPLK
        IF( ISPLK(14+NSPKX).LT.0 ) THEN
          NSPX = ABS(ISPLK(14+NSPKX))
          IF( NSPX.GT.1000 ) NSPX = NSPX - 1000
          IF( NEQ.EQ.IEQ_S(NSPX) ) THEN
            AJM(NROW,NROW) = 1.D+0
            BJM(NROW) = 0.D+0
            ISUB_LOG = ISUB_LOG-1
            RETURN
          ENDIF
        ENDIF
      ENDDO
!      
      L1: DO M = 1,IEQ_C(1,NEQX)
        NSP = IEQ_C(M+1,NEQX)
        NEQXX = IEQ_S(NSP)
        BJM(NROW) = BJM(NROW) - EQ_C(M,NEQX)*SP_CX(NSP)*VTOLX
        IF( NEQXX.GT.NEQE ) THEN
!
!---      Skip for initial pH  ---
!
          IF( ISPLK(1).GT.1000 ) THEN
            IF( NSP.EQ.MOD(ISPLK(1),1000) .AND. NSTEP.EQ.0 ) CYCLE L1
          ENDIF
!
!---      Skip fixed species concentration
!
          DO NSPKX = 1,NSPLK
            IF( ISPLK(14+NSPKX).LT.0 ) THEN
              NSPX = ABS(ISPLK(14+NSPKX))
              IF( NSP.EQ.NSPX) CYCLE L1
            ENDIF
         ENDDO
         NCOL = NEQXX-NEQE
         AJM(NROW,NCOL) = AJM(NROW,NCOL)- EQ_C(M,NEQX)*VTOLX
         DO NSPKX = 1,NSPLK
           IF( ISPLK(14+NSPKX).LT.0 ) THEN
             NSPX = ABS(ISPLK(14+NSPKX))
             IF( NSPX .GT. 1000 ) NSPX = NSPX - 1000
             IF( NSP.EQ.NSPX ) AJM(NROW,NCOL) = AJM(NROW,NCOL)-
     &         1.D0/ACTVX(1,NSP)*EQ_C(M,NEQX)*VTOLX
            ENDIF
         ENDDO
        ELSE
!
!---      Equilibrium mass action
!
          NSE = IEQ_E(1,NEQXX)
!
!---      Loop over 1 species  ---
!
          DO MSX = 2,NSE
            NCM_SP = IEQ_E(MSX+1,NEQXX)
            NCOL = IEQ_S(NCM_SP)-NEQE
            AJM(NROW,NCOL) = AJM(NROW,NCOL)-EQ_C(M,NEQX)*VTOLX*
     &        EQ_E(MSX-1,NEQXX)*SP_CX(NSP)/SP_CX(NCM_SP)
          ENDDO
        ENDIF
      ENDDO L1
!
!---  Return residual vector and Jacobian matrix  ---
!
      IF( ABS(INDX).EQ.1 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      BJM(NROW) = -BJM(NROW)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CECHEM_R group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DAVIES( ACTVX,SP_CX,DSP_CX,SLX,PORDX,RHOLX,TX,XLWX )
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
!     ECKEChemX
! 
!     B-Dot or extended Debye-Huckel model for activity coefficient.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 ACTVX(LSPL+1,LSPL),SP_CX(LSPL),DSP_CX(LSPL)
      REAL*8 CPIX(LSPL+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DAVIES'
!
!---  Ionic strength of the aqueous solution  ---
!
      DO M = 1,NSPL+1
        CPIX(M) = 0.D+0
        SUM_MZ = 0.D+0
        SUM_M = 0.D+0
        DO NSP = 1,NSPL
!
!---      Skip neutral species  ---
!
          IF( ABS(SP_L(1,NSP)).LT.EPSL ) CYCLE
          IF( NSP.EQ.(M-1) ) THEN
            CLX = SP_CX(NSP) + DSP_CX(NSP)
          ELSE
            CLX = SP_CX(NSP)
          ENDIF
!
!---      Molarity in mol solute/m^3 aqueous
!         or mol solute/l aqueous  ---
!
          CMX = CLX/(SLX*PORDX)
!
!---      Molality in mol solute/kg water  ---
!
          CMX = CMX/(RHOLX*XLWX)
          SUM_M = SUM_M + CMX
          SUM_MZ = SUM_MZ + CMX*SP_L(1,NSP)
          CPIX(M) = CPIX(M) + CMX*(SP_L(1,NSP)**2)
        ENDDO
!
!---    Correct for electrical neutrality  ---
!
        CPIX(M) = CPIX(M)+DABS(SUM_MZ)
        CPIX(M) = CPIX(M)*5.D-1
      ENDDO
!
!---  Activity coefficients for aqueous species  ---
!
      DO M = 1,NSPL+1
        DO NSP = 1,NSPL
!
!---      Activity coefficients for charged species  ---
!
          IF( ABS(SP_L(1,NSP)).GT.EPSL ) THEN
            ACTVX(M,NSP) = 5.D-1*(SP_L(1,NSP)**2)*(SQRT(CPIX(M))
     &        /(1.D+0 + SQRT(CPIX(M)))- 2.4D-1*CPIX(M))
            IF(ACTVX(M,NSP) > 1.D+1) ACTVX(M,NSP) = 0.D+0
            ACTVX(M,NSP) = 1.D+1**(-ACTVX(M,NSP))
!
!---      Activity coefficients for electrically neutral species  ---
!
          ELSE
            ACTVX(M,NSP) = 1.D+0
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DAVIES group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ECKECHEM
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
!     ECKEChemX
! 
!     Equilibrium, Conservation, and Kinetic Equation CHEMistry.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNNL, 2 February 2022
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 SP_CX(LSPR),SP_CXX(LSPR)
      REAL*8 CX(LEQC+LEQK),COX(LEQC+LEQK)
      REAL*8 ACTCX(LSPL)
      REAL*8 ACTVX(LSPL+1,LSPL)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ECKECHEM'
      ECKE_ER = .FALSE.
!
!---  Total number of equations  ---
!
      NEQR = NEQE + NEQK + NEQC
!
!---  Total number species  ---
!
      NSPR = NSPG + NSPL + NSPN + NSPS + NSPE
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      ND_ERL = NFLD+1
      L1: DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        IF( SL(2,N).LT.1.D-4 ) CYCLE
        N_DB = ND(N)
!
!---    Assign local values of old-time-step species concentration  ---
!
        DO NSP = 1,NSPR
          IF( ISLC(57).EQ.0 ) THEN
            SP_CX(NSP) = SP_C(N,NSP)
          ELSE
            SP_CX(NSP) = MAX(SP_C(N,NSP),CMIN)
          ENDIF
        ENDDO
!
!---    Assign local values of component species concentration  ---
!
        DO NSP = 1,NEQC
          NSL = NSP + NSOLU
          CX(NSP) = C(N,NSL)
          COX(NSP) = C(N,NSL)
        ENDDO
!
!---    Assign local values of kinetic species concentration  ---
!
        DO NSP = 1,NEQK
          NSPX = NSP + NEQC
          NSL = NSPX + NSOLU
          CX(NSPX) = C(N,NSL)
          COX(NSPX) = C(N,NSL)
        ENDDO
!
!---    Guess specie concentrations, assuming activity
!       coefficients of 1.0  ---
!
        IF( ISLC(42).EQ.1 .AND. NSTEP.EQ.0 ) THEN
!
!---      Check if neighboring cell has identical
!         conservation component species concentrations  ---
!
          ISKIP = 1
          IF( N.GT.1 .AND. NFCGC(ID+1).GT.1 ) THEN
            DO NSL = NSOLU+1,NSOLU+NEQC
              IF( ABS(CO(N-1,NSL)-C(N,NSL))/EPSL.GT.EPSL ) THEN
                ISKIP = 0
                EXIT
              ENDIF
            ENDDO
!
!---        Skip guess calculation, use neighbor-cell species
!           concentrations  ---
!
            IF( ISKIP.EQ.1 ) THEN
              DO NSP = 1,NSPR
                SP_CX(NSP) = SP_C(N-1,NSP)
              ENDDO
            ELSE
              CALL GUSPCN( CX,SP_CX,N )
            ENDIF
          ELSE
            CALL GUSPCN( CX,SP_CX,N )
          ENDIF
        ENDIF
!
!---    Conventional Newton method (i.e., global solve)  ---
!
        IF (ISLC(57).EQ.0) THEN
          CALL ECKECN( ACTCX,CX,COX,SP_CX,N )
!
!---    Reduced Equations  ---
!
        ELSE
          CALL ECKECN_R( ACTCX,CX,COX,SP_CX,N )
        ENDIF
!
!---    Convergence failure, store failed node number exit loop  ---
!
        IF( ECKE_ER ) THEN
          ND_ERL = ND(N)
          EXIT L1
        ENDIF
!
!---   Fixed species activity
!
        DO NSP = 1,NSPR
          SP_CXX(NSP) = SP_CX(NSP)
          DO NSPKX = 1,NSPLK
           IF( ISPLK(14+NSPKX).LT.-1000 ) THEN
             NSPX = ABS(ISPLK(14+NSPKX))-1000
             IF(NSP.EQ.NSPX) SP_CX(NSP) = SP_CX(NSP)/ACTVX(1,NSP)
            ENDIF
          ENDDO
        ENDDO
!
!---    Conservation-component species, mobile conservation-component
!       species, and old time step mobile conservation-component
!       species  ---
!
        DO NEQ = 1,NEQC
          NSL = NSOLU + NEQ
          C(N,NSL) = 0.D+0
          CMBX = 0.D+0
!
!---      Loop over conservation-component species  ---
!
          DO M = 1,IEQ_C(1,NEQ)
            NSP = IEQ_C(M+1,NEQ)
            IF( ABS(SP_CX(NSP)).LT.CMIN ) THEN
              SP_CX(NSP) = 0.D+0
            ENDIF
            C(N,NSL) = C(N,NSL) + EQ_C(M,NEQ)*SP_CX(NSP)
!
!---        Mobile species ---
!
            IF( NSP.LE.NSPL .OR. NSP.GT.(NSPL+NSPS+NSPE) )
     &        CMBX = CMBX + EQ_C(M,NEQ)*SP_CX(NSP)
!
!---        Immobile species ---
!
            IF( NSP.GT.NSPL .AND. NSP.LE.(NSPL+NSPS+NSPE) ) THEN
              IF( ABS(SP_CO(N,NSP)).LT.CMIN ) THEN
                SP_COX = 0.D+0
              ELSE
                SP_COX = SP_CO(N,NSP)
              ENDIF
              SPI_CX = EQ_C(M,NEQ)*SP_CX(NSP)
              SPI_COX = EQ_C(M,NEQ)*SP_COX
!
!---          Load CO2 sources for linked aqueous CO2,
!             skip for initial conditions  ---
!
              IF( ISPLK(6).EQ.NSL .AND. (NSTEP-NRST).GT.0 ) THEN
                SRCAX = -1.D-3*(SPI_CX-SPI_COX)*VOL(N)*WTMA
                IF( ABS(SRCAX).LT.1.D-16 ) SRCAX = 0.D+0
                SRCA(1,N) = SRCA(1,N) + SRCAX
              ENDIF

!
!---          Load air sources for linked aqueous air,
!             skip for initial conditions  ---
!
              IF( ISPLK(4).EQ.NSL .AND. (NSTEP-NRST).GT.0 ) THEN
                SRCAX = -1.D-3*(SPI_CX-SPI_COX)*VOL(N)*WTMA
                IF( ABS(SRCAX).LT.1.D-16 ) SRCAX = 0.D+0
                SRCA(1,N) = SRCA(1,N) + SRCAX
              ENDIF
            ENDIF
          ENDDO
        ENDDO
!
!---    Kinetic-component species, mobile kinetic-component species
!       and old time step mobile kinetic-component species  ---
!
        DO NEQ = 1,NEQK
          NSL = NSOLU + NEQC + NEQ
          C(N,NSL) = 0.D+0
          CMBX = 0.D+0
!
!---      Loop over kinetic-component species  ---
!
          DO M = 1,IEQ_K(1,NEQ)
            NSP = IEQ_K(M+1,NEQ)
            C(N,NSL) = C(N,NSL) + EQ_K(M,NEQ)*SP_CX(NSP)
!
!---        Mobile species ---
!
            IF( NSP.LE.NSPL .OR. NSP.GT.(NSPL+NSPS+NSPE) )
     &        CMBX = CMBX + EQ_K(M,NEQ)*SP_CX(NSP)
!
!---        Immobile species ---
!
            IF( NSP.GT.NSPL .AND. NSP.LE.(NSPL+NSPS+NSPE) ) THEN
              IF( ABS(SP_CO(N,NSP)).LT.CMIN ) THEN
                SP_COX = 0.D+0
              ELSE
                SP_COX = SP_CO(N,NSP)
              ENDIF
              SPI_CX = EQ_K(M,NEQ)*SP_CX(NSP)
              SPI_COX = EQ_K(M,NEQ)*SP_COX
!
!---          Load CO2 sources for linked aqueous CO2,
!             skip for initial conditions  ---
!
              IF( ISPLK(6).EQ.NSL .AND. (NSTEP-NRST).GT.0 ) THEN
                SRCAX = -1.D-3*(SPI_CX-SPI_COX)*VOL(N)*WTMA
                IF( ABS(SRCAX).LT.1.D-16 ) SRCAX = 0.D+0
                SRCA(1,N) = SRCA(1,N) + SRCAX
              ENDIF

!
!---          Load air sources for linked aqueous air,
!             skip for initial conditions  ---
!
              IF( ISPLK(4).EQ.NSL .AND. (NSTEP-NRST).GT.0 ) THEN
                SRCAX = -1.D-3*(SPI_CX-SPI_COX)*VOL(N)*WTMA
                IF( ABS(SRCAX).LT.1.D-16 ) SRCAX = 0.D+0
                SRCA(1,N) = SRCA(1,N) + SRCAX
              ENDIF
            ENDIF
          ENDDO
        ENDDO
!
!---    Assign global species concentrations  ---
!
        DO NSP = 1,NSPR
          IF( ABS(SP_CX(NSP)).LT.CMIN ) THEN
            SP_C(N,NSP) = 0.D+0
          ELSE
            SP_C(N,NSP) = SP_CX(NSP)
          ENDIF
!
!---   Fixed species activity
!
          DO NSPKX = 1,NSPLK
            IF( ISPLK(14+NSPKX).LT.-1000 ) THEN
              NSPX = ABS(ISPLK(14+NSPKX))-1000
              SP_C(N,NSP) = SP_CXX(NSP)
            ENDIF
          ENDDO
        ENDDO
!
!---    Loop over solid species to determine current mineral
!       volume fractions and porosity
!       POR0(1,N) total porosity
!       POR0(2,N) diffusive porosity
!       POR0(3,N) total unreactive mineral volume   ---
!
        POR0(1,N) = 1.D+0 - POR0(3,N)
        POR0(2,N) = 1.D+0 - POR0(3,N) - POR(2,N) + POR(1,N)
        DO NSP = NSPL+1,NSPL+NSPS
          NSP_M = NSP-NSPL
!
!---      Minerals ---
!
          IF( ISP_MN(NSP).EQ.1 ) THEN
            IF( SP_CX(NSP)+SP_CMN(N,NSP_M).LT.CMIN ) THEN
              RS_S(3,NSP_M,N) = 0.D+0
            ELSE
              RS_S(3,NSP_M,N) = 1.D-3*(SP_CX(NSP)+SP_CMN(N,NSP_M))*
     &          SP_S(2,NSP_M)/SP_S(1,NSP_M)
            ENDIF
!
!---      Non minerals ---
!
          ELSE
            IF( SP_CX(NSP).GT.CMIN ) THEN
              RS_S(3,NSP_M,N) = 1.D-3*SP_CX(NSP)*SP_S(2,NSP_M)/
     &          SP_S(1,NSP_M)
            ELSE
              RS_S(3,NSP_M,N) = 0.D+0
            ENDIF
          ENDIF
          RS_S(3,NSP_M,N) = MIN(RS_S(3,NSP_M,N),1.D+0)
          POR0(1,N) = POR0(1,N) - RS_S(3,NSP_M,N)
          POR0(2,N) = POR0(2,N) - RS_S(3,NSP_M,N)
          POR0(1,N) = MAX(POR0(1,N),1.D-12)
          POR0(2,N) = MAX(POR0(2,N),1.D-12)
        ENDDO
!
!---    Load pH variable for output  ---
!
        IF( ISPLK(1).GT.1 ) THEN
          NSPH = MOD(ISPLK(1),1000)
          IF( SL(2,N)*PORD(2,N).GT.SMALL ) THEN
            VARX = ACTCX(NSPH)*SP_C(N,NSPH)/(SL(2,N)*PORD(2,N))
          ELSE
            VARX = 0.D+0
          ENDIF
          VARX = VARX*1.D-3
          C_PH(N) = -LOG10( VARX )
        ENDIF
      ENDDO L1
      CALL MPI_ALLREDUCE( ND_ERL,ND_ERG,1,MPI_INTEGER,MPI_MIN,
     &  MPI_COMM_WORLD,IERR )
!
!---  Convergence failure, double the number of
!     sub-time steps  ---
!
      IF( ND_ERG.LT.(NFLD+1) ) THEN
        CALL RESET_SP
        N_RST = 2*N_RST
        REALX = REAL(N_RST)
        DT = DT_RST/REALX
        DTI = 1.D+0/(DT+EPSL)
        TM = TM_RST - DT_RST
        IF( ID.EQ.0 ) THEN
          PRINT *,'          ---  ECKEChem Convergence Failure, Node = '
     &      ,ND_ERG,'  ---'
          PRINT *,'          ---  ECKEChem Sub-Time Stepping Factor = '
     &      ,N_RST,'  ---'
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ECKECHEM group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ECKECN( ACTCX,CX,COX,SP_CX,N )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
! 
!     Conventional Newton scheme for solving the nonlinear
!     1, conservation, and kinetic chemistry equations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 2 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
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
      REAL*8 AJM(LSPR,LSPR),BJM(LSPR)
      REAL*8 SP_CX(LSPR)
      REAL*8 CX(LEQC+LEQK),COX(LEQC+LEQK)
      REAL*8 DCLX(LSPL),CLX(LSPL),ACTVX(LSPL+1,LSPL)
      REAL*8 ACTCX(LSPL)
      REAL*8 RSDX(LSPR)
      REAL*8 RSDMAX(128)
      INTEGER IJM(LSPR)
      CHARACTER(64) GETSPNM
      CHARACTER(32) :: CHMSG
      EXTERNAL GETSPNM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ECKECN'
!
!---  Total number of equations  ---
!
      NEQR = NSPR
      RLXF = 1.D+0
      NCX = 128
      IF( NSTEP-NRST.EQ.0 ) NCX = 128
      DO NC = 1,NCX
        RSDMAX(NC) = 0.D+0
      ENDDO
      RSD2MAX = 0.D+0
!
!--- Fixed species activity ---
!
      IF( NSTEP.EQ.0 ) THEN
        DO NSPKX = 1,NSPLK
          IF( ISPLK(14+NSPKX).LT.0 ) THEN
            NSPX = ABS(ISPLK(14+NSPKX))
            IF( NSPX.GT.1000 ) THEN
              NSPX = NSPX - 1000
              FACTV(NSPX) = SP_CX(NSPX)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
!
!---  For initial conditions, skip kinetic equations  ---
!
      IF( NSTEP.EQ.0 ) NEQR = NEQE+NEQC
!
!---  Top of Newton-Raphson loop  ---
!
      NC = 0
      L1: DO
        NC = NC + 1
        RSDMAX(NC) = 0.D+0
!
!---    Set index to compute both the Jacobian matrix
!       and residual vector  ---
!
        INDX = 0
        CALL ECKEJCB( ACTVX,AJM,BJM,CX,COX,SP_CX,N,INDX,IER )
!
!---    Solve linear system  ---
!
        IF( ISLC(78).EQ.1 ) THEN
          CALL LUDCMP( AJM,NEQR,LSPR,IJM,DJM )
          CALL LUBKSB( AJM,NEQR,LSPR,IJM,BJM )
!        ELSEIF( ISLC(78).EQ.2 ) THEN
!          CALL DGELG( BJM,AJM,LSPR,NEQR,EPSL,IER )
        ELSE
          CALL GAUSSJ( AJM,BJM,NEQR,LSPR,IER )
        ENDIF
        IF( IER.EQ.-1 ) EXIT L1
!
!---    Maximum residual  ---
!
        IC = 0
        DO NSP = 1,NSPR
          RSDX(NSP) = 0.D+0
          IF( ABS(SP_CX(NSP))/EPSL.GT.EPSL .OR.
     &      ABS(BJM(IEQ_S(NSP)))/EPSL.GT.EPSL ) THEN
            SP_CMX = MAX( ABS(SP_CX(NSP)),ABS(SP_C(N,NSP)),1.D-12 )
!
!---        Check for disappearing 1 or conservation 
!           species  --
!
            IF( NSP.LE.(NEQE+NEQC) ) THEN
              VARX = SP_CX(NSP)+BJM(IEQ_S(NSP))
              IF( VARX.LT.0.D+0 .AND. 
     &          ABS(VARX).GT.1.D+3*SP_CX(NSP) ) THEN
                RSDX(NSP) = 0.D+0
              ELSEIF( ABS(SP_CMX).GT.CMIN ) THEN
                RSDX(NSP) = ABS(BJM(IEQ_S(NSP)))/SP_CMX
              ELSE
                RSDX(NSP) = 0.D+0
              ENDIF
            ELSE
              IF( ABS(SP_CMX).GT.CMIN ) THEN
                RSDX(NSP) = ABS(BJM(IEQ_S(NSP)))/SP_CMX
              ELSE
                RSDX(NSP) = 0.D+0
              ENDIF
            ENDIF
            IF( INDEX(GETSPNM(NSP),'fix').NE.0 ) THEN
              RSDX(NSP) = 0.D+0
            ENDIF
            IF( RSDX(NSP).GT.RSDMAX(NC) ) THEN
              IF( ABS(SP_CX(NSP)).GT.1.D-16 ) THEN
                RSD2MAX = RSDMAX(NC)
                RSDMAX(NC) = RSDX(NSP)
                NSPMX = NSP
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        ISKIP = 0
        IF( NC.EQ.1 .AND. NSTEP-NRST.EQ.0 ) ISKIP = 1
!
!---    Update the species concentrations, mol/m^3  ---
!
        DO NSP = 1,NSPR
!
!---      Relaxation adjusted to residual  ---
!
          RLXF = 1.5D+0 - 7.5D-1*(RSDX(NSP)**1.D-1)
          RLXF = MAX( 1.D-2,MIN( 1.D+0,RLXF ) )
          IF( NC.EQ.1 ) RLXF = 1.D+0
          IF( NC.GT.16 ) RLXF = 6.D-1
!
!---      Relaxation/Check for NaNs  ---
!
          CCX = RLXF*BJM(IEQ_S(NSP))
          IF( CCX.NE.CCX )THEN
            IER = -1
            EXIT L1
          ENDIF
!
!---      Mineral species  ---
!
          IF( ISP_MN(NSP).EQ.1 ) THEN
!
!---        Negative adjustments to species concentrations  ---
!
            IF( CCX.LT.0.D+0 ) THEN
              NSP_MN = NSP - NSPL
              SP_CCX = MAX( SP_CMN(N,NSP_MN)+SP_CX(NSP),ABS(SP_CX(NSP)))
              SP_CCX = 0.6D+0*SP_CCX
              CCX = SIGN( MIN( SP_CCX,ABS(CCX) ),CCX )
            ENDIF
!
!---      Immobile species  ---
!
          ELSEIF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
!
!---        Negative adjustments to species concentrations  ---
!
            IF( CCX.LT.0.D+0 ) THEN
              SP_CCX = SP_CX(NSP)
              SP_CCX = 0.6D+0*SP_CX(NSP)
              CCX = SIGN( MIN( SP_CCX,ABS(CCX) ),CCX )
            ELSEIF( ISKIP.EQ.0 ) THEN
              SP_CCX = 1.D+1*SP_CX(NSP)
              CCX = SIGN( MIN( SP_CCX,ABS(CCX) ),CCX )
            ENDIF
!
!---      Mobile species  ---
!
          ELSE
!
!---        Negative adjustments to species concentrations  ---
!
            IF( CCX.LT.0.D+0 ) THEN
              SP_CCX = 0.6D+0*SP_CX(NSP)
              CCX = SIGN( MIN( SP_CCX,ABS(CCX) ),CCX )
            ELSEIF( ISKIP.EQ.0 ) THEN
              SP_CCX = 1.D+1*SP_CX(NSP)
              CCX = SIGN( MIN( SP_CCX,ABS(CCX) ),CCX )
            ENDIF
          ENDIF
!
!---      Update mineral species concentration  ---
!
          IF( ISP_MN(NSP).EQ.1 ) THEN
            SP_CX(NSP) = SP_CX(NSP)+CCX
            IF( ABS(SP_CX(NSP)).LT.CMIN ) SP_CX(NSP) = 0.D+0
!
!---      Update non-mineral species concentration  ---
!
          ELSE
!
!---      Fix species concentration
!
            IF( INDEX(GETSPNM(NSP),'fix').NE.0 ) THEN
              CCX = 0.D+0
              IF ( NSTEP-NRST.GT.0 ) THEN
!
!---            Convert reactive species from node volumetric, kmol/m^3
!               to node volumetric, mol/m^3  ---
!
                IF( MOD(IC_SP(N,NSP),10).EQ.1 ) THEN
                  SP_CX(NSP) = 1.D+3*SP_CI(N,NSP)
!
!---            Convert reactive species from aqueous volumetric,
!                kmol/m^3 to node volumetric, mol/m^3  ---
!
                ELSEIF( MOD(IC_SP(N,NSP),10).EQ.2 ) THEN
                  SP_CX(NSP) = 1.D+3*SP_CI(N,NSP)*SL(2,N)*PORD(2,N)
!
!---            Convert reactive species from aqueous molal, 
!               kmol/kg water to node volumetric, mol/m^3  ---
!
                ELSEIF( MOD(IC_SP(N,NSP),10).EQ.3 ) THEN
                  SP_CX(NSP) = 1.D+3*SP_CI(N,NSP)*SL(2,N)*PORD(2,N)*
     &              RHOL(2,N)*XLW(2,N)
!
!---            Convert reactive species from gas volumetric, kmol/m^3
!               to node volumetric, mol/m^3  ---
!
                ELSEIF( MOD(IC_SP(N,NSP),10).EQ.4 ) THEN
                  SP_CX(NSP) = 1.D+3*SP_CI(N,NSP)*SG(2,N)*PORD(2,N)
                ENDIF
              ENDIF
            ENDIF
            SP_CX(NSP) = SP_CX(NSP)+CCX
            IF( ABS(SP_CX(NSP)).LT.CMIN ) SP_CX(NSP) = 0.D+0
          ENDIF
        ENDDO
!
!---    Fixed species activity  ---
!
        DO NSPKX = 1,NSPLK
          IF( ISPLK(14+NSPKX).LT.0 ) THEN
            NSPX = ABS(ISPLK(14+NSPKX))
            IF( NSPX.GT.1000 ) THEN
              NSPX = NSPX - 1000
              IF( IACTV.EQ.3 ) THEN
               ACTVX(1,NSP) = ACTVC
              ELSEIF( IACTV.EQ.1 ) THEN
               CALL DAVIES( ACTVX,CLX,DCLX,SL(2,N),PORD(2,N),
     &          RHOL(2,N),T(2,N),XLW(2,N) )
              ELSE
               CALL BDOT( ACTVX,CLX,DCLX,SL(2,N),PORD(2,N),
     &           RHOL(2,N),T(2,N),XLW(2,N) )
              ENDIF
              SP_CX(NSPX) = FACTV(NSPX)/ACTVX(1,NSPX)
            ENDIF
          ENDIF
        ENDDO
!
!---    Unconverged species concentration exit or repeat
!       Newton-Raphson loop  ---
!
        IF( RSDMAX(NC).GT.1.D-6 ) THEN
!
!---      Maximum iteration count  ---
!
          IF( NC.EQ.NCX ) THEN
!
!---        Small second maximum and reasonably small maximum check  ---
!
            IF( RSD2MAX.GT.1.D-6 .OR. RSDMAX(NC).GT.1.D-1 ) THEN
              ECKE_ER = .TRUE.
            ENDIF 
          ELSE
            CYCLE L1
          ENDIF
!
!---    Converged solution  ---
!
        ELSE
          EXIT L1
        ENDIF
      ENDDO L1
      IF(IER.EQ.-1) ECKE_ER = .TRUE.
!
!---  Load activity coefficient for output  ---
!
      DO NSP = 1,NSPL
        ACTCX(NSP) = ACTVX(1,NSP)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ECKECN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ECKECN_R( ACTCX,CX,COX,SP_CX,N )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
! 
!     Conventional Newton scheme for solving the nonlinear
!     1, conservation, and kinetic chemistry equations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 3 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
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
      REAL*8 AJM(LSPR,LSPR),BJM(LSPR)
      REAL*8 SP_CX(LSPR),EQKX(LEQE)
      REAL*8 CX(LEQC+LEQK),COX(LEQC+LEQK)
      REAL*8 ACTCX(LSPL)
      REAL*8 ACTVX(LSPL+1,LSPL)
      REAL*8 RSDX(LSPR)
      REAL*8 RSDMAX(128)
      CHARACTER(32) :: CHMSG
      INTEGER IJM(LSPR)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ECKECN_R'
!
!---  Total number of equations  ---
!
      NEQR = NEQC+NEQK
      RLXF = 1.D+0
      NCX = 64
      IF( NSTEP-NRST.EQ.0 ) NCX = 128
      DO NC = 1,NCX
        RSDMAX(NC) = 0.D+0
      ENDDO
!
!---  Volumetric concentration to molality  ---
!
      VTOMX = 1.D+0/(SL(2,N)*PORD(2,N)*RHOL(2,N)*XLW(2,N))
!
!---  Equilibrium constant as a function of temperature  ---
!
      DO NEQ = 1,NEQE
        NS = IEQ_E(1,NEQ)
        IRCE = IEQ_E((NS+2),NEQ)
        CALL EQCN( EQKX(NEQ),T(2,N),IRCE,N )
      ENDDO
!
!---  For initial conditions, skip kinetic equations  ---
!     YFang -- Keep kinetic variables
!      IF( NSTEP.EQ.0 ) NEQR = NEQE+NEQC
!
!---  Guess specie concentrations, assuming activity
!     coefficients of 1.0  ---
!
      IF( NSTEP.EQ.NRST ) THEN
!
!---    Assign local values of old-time-step species concentration  ---
!
        DO NSP = 1,NSPR
          IF( ABS(SP_C(N,NSP)).LT.CMIN ) THEN
!            SP_CX(NSP) = 0.D+0
            SP_CX(NSP) = MAX(SP_CX(NSP),CMIN)
          ELSE
            SP_CX(NSP) = SP_C(N,NSP)
          ENDIF
        ENDDO
!
!---    Assign local values of component species concentration  ---
!
        DO NSP = 1,NEQC
          NSL = NSP + NSOLU
          CX(NSP) = C(N,NSL)
          COX(NSP) = C(N,NSL)
        ENDDO
!
!---    Assign local values of kinetic species concentration  ---
!
        DO NSP = 1,NEQK
          NSPX = NSP + NEQC
          NSL = NSPX + NSOLU
          CX(NSPX) = C(N,NSL)
          COX(NSPX) = C(N,NSL)
        ENDDO
        NSPH = MOD(ISPLK(1),1000)
        DO NEQ = 1,NEQC
          NSP = IEQ_C(2,NEQ)
          IF( NSP.EQ.NSPH ) CYCLE
          SP_CX(NSP) = MAX( CX(NEQ)*5.D-2,CMIN )
        ENDDO
!
!---  If uninitialized set pH to 1.D-7 molality, mol/kg water  ---
!
        IF( NSPH.GE.1 .AND. NSPH.LE.LSPR ) THEN
          IF( SP_CX(NSPH)/EPSL.LT.EPSL ) SP_CX(NSPH) = 1.D-7/VTOMX
        ENDIF
!
!---  Recalculate 1 species concentrations 
!
        DO NEQ = 1,NEQE
          NEQ_SP = IEQ_E(2,NEQ)
          NS = IEQ_E(1,NEQ)
!
!---      Equilibrium constant and exponent  ---
!
          IF( ISLC(60).EQ.0 ) THEN
            C_NEG = EQKX(NEQ)**EQ_E(NS,NEQ)
          ELSE
            C_NEG = EQ_E(NS,NEQ)*LOG(EQKX(NEQ))
          ENDIF
!
!---      Loop over species in 1 equation
!         skipping the 1 species  ---
!
          DO K = 1,NS
            NCM_SP = IEQ_E(K+1,NEQ)
!
!---        Skip 1 species  ---
!
            IF( NCM_SP.EQ.NEQ_SP ) CYCLE
!
!---        Convert conservation species concentration to
!           molality, mol/kg H2O  ---
!
            CMX = SP_CX(NCM_SP)*VTOMX
            IF( ISLC(60).EQ.0 ) THEN
              C_NEG = C_NEG*(CMX**EQ_E(K-1,NEQ))
            ELSE
              C_NEG = C_NEG+EQ_E(K-1,NEQ)*LOG(CMX)
            ENDIF
         ENDDO
         IF( ISLC(60).EQ.1 ) THEN
           C_NEG = EXP(C_NEG)
         ENDIF
!
!---    Convert 1 species concentration to
!       node volumetric, mol/m^3  ---
!
         SP_CX(NEQ_SP) = C_NEG/VTOMX
         SPCXX = MIN( 1.0D+30,SP_CX(NEQ_SP) )
         DO NCX = 1,NEQC
           L1: DO IX=1,IEQ_C(1,NCX)
             NSP=IEQ_C(IX+1,NCX)
             IF(NEQ_SP.EQ.NSP) THEN
               SPCXM=ABS(CX(NCX)/EQ_C(IX,NCX))
               SPCXX=MIN(SPCXX,SPCXM)
               EXIT L1
             ENDIF
           ENDDO L1
         ENDDO
         SP_CX(NEQ_SP) = MAX(SPCXX,CMIN)
        ENDDO
      ENDIF
!
!---  Top of Newton-Raphson loop  ---
!
      NC = 0
      L2: DO
        NC = NC + 1
!
!---    Set index to compute both the Jacobian matrix
!       and residual vector  ---
!
        INDX = 0
        CALL ECKEJCB( ACTVX,AJM,BJM,CX,COX,SP_CX,N,INDX,IER )
!
!---    Solve linear system  ---
!
        IF( ISLC(78).EQ.1 ) THEN
          CALL LUDCMP( AJM,NEQR,LSPR,IJM,DJM )
          CALL LUBKSB( AJM,NEQR,LSPR,IJM,BJM )
!        ELSEIF( ISLC(78).EQ.2 ) THEN
!          CALL DGELG( BJM,AJM,LSPR,NEQR,EPSL,IER )
        ELSE
          CALL GAUSSJ( AJM,BJM,NEQR,LSPR,IER )
        ENDIF
        IF( IER.EQ.-1 ) EXIT L2
!
!---    Maximum residual  ---
!
        RSDMAX(NC) = 1.D-20
        DO NSP = 1,NSPR
          RSDX(NSP) = 0.D+0
          NEQX = IEQ_S(NSP)
          IF( NEQX.GT.NEQE) THEN
            NROWX = NEQX-NEQE
            IF( SP_CX(NSP)/EPSL.GT.EPSL .OR.
     &        BJM(NROWX)/EPSL.GT.EPSL ) THEN
              SP_CMX = MAX( ABS(SP_CX(NSP)),ABS(SP_C(N,NSP)),CMIN )
              IF( ABS(SP_CMX).GT.1.D-16 ) THEN
                RSDX(NSP) = ABS(BJM(NROWX))/SP_CMX
              ELSE
                RSDX(NSP) = 0.D+0
              ENDIF
              IF( RSDX(NSP).GT.RSDMAX(NC) ) THEN
                IF( ABS(SP_CX(NSP)).GT.1.D-16 ) THEN
                  RSDMAX(NC) = RSDX(NSP)
                  NSPMX = NSP
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
!
!---    Update the species concentrations, mol/m^3  ---
!
        DO NSP = 1,NSPR
!
!---      Relaxation  ---
!
          NEQX = IEQ_S(NSP)
          IF( NEQX.GT.NEQE ) THEN
            NROWX = NEQX-NEQE
            CCX = RLXF*BJM(NROWX)
            IF( CCX.LT.0.D+0 ) THEN
!
!---          Mineral species  ---
!
              IF( ISP_MN(NSP).EQ.1 ) THEN
                NSP_MN = NSP - NSPL
                SP_CCX = SP_CMN(N,NSP_MN)+SP_CX(NSP)
                CCX = SIGN( MIN( SP_CCX,ABS(CCX) ),CCX )
!
!---          Immobile species  ---
!
              ELSEIF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
                SP_CCX = SP_CX(NSP)
                CCX = SIGN( MIN( SP_CCX,ABS(CCX) ),CCX )
!
!---          Mobile species  ---
!
              ELSE
                SP_CCX = 9.5D-1*SP_CX(NSP)
                CCX = SIGN( MIN( SP_CCX,ABS(CCX) ),CCX )
              ENDIF
            ENDIF
!
!---        Update mineral species concentration  ---
!
            IF( ISP_MN(NSP).EQ.1 ) THEN
              SP_CX(NSP) = SP_CX(NSP)+CCX
              IF( ABS(SP_CX(NSP)).LT.CMIN ) SP_CX(NSP) = 0.D+0
!
!---        Update non-mineral species concentration  ---
!
            ELSE
              SP_CX(NSP) = MAX( SP_CX(NSP)+CCX,0.D+0 )
            ENDIF
          ENDIF
        ENDDO
!
!---    Recalculate 1 species concentrations 
!
        DO NEQ = 1,NEQE
          NEQ_SP = IEQ_E(2,NEQ)
          NS = IEQ_E(1,NEQ)
!
!---      Equilibrium constant and exponent  ---
!
          IF( ISLC(60).EQ.0 ) THEN
            C_NEG = EQKX(NEQ)**EQ_E(NS,NEQ)
          ELSE
            C_NEG = EQ_E(NS,NEQ)*LOG(EQKX(NEQ))
          ENDIF
!
!---      Loop over species in 1 equation
!         skipping the 1 species  ---
!
          DO K = 1,NS
            NCM_SP = IEQ_E(K+1,NEQ)
!
!---        Skip 1 species  ---
!
            IF( NCM_SP.EQ.NEQ_SP ) CYCLE
!
!---        Convert conservation species concentration to
!           molality, mol/kg H2O  ---
!
            IF(NCM_SP.LE.NSPL) THEN
              ACX = ACTVX(1,NCM_SP)
            ELSE
              ACX = 1.D+0
            ENDIF
            CMX = SP_CX(NCM_SP)*VTOMX*ACX
            IF( ISLC(60).EQ.0 ) THEN
              C_NEG = C_NEG*(CMX**EQ_E(K-1,NEQ))
            ELSE
              C_NEG = C_NEG+EQ_E(K-1,NEQ)*LOG(CMX)
            ENDIF
          ENDDO
          IF( ISLC(60).EQ.1 ) THEN
            C_NEG = EXP(C_NEG)
          ENDIF
!
!---      Convert 1 species concentration to
!         node volumetric, mol/m^3  ---
!
          IF(NEQ_SP.LE.NSPL) THEN
            ACX = ACTVX(1,NEQ_SP)
          ELSE
            ACX = 1.D+0
          ENDIF
          SP_CX(NEQ_SP) = C_NEG/VTOMX/ACX
        ENDDO
!
!---    Unconverged species concentration exit or repeat
!       Newton-Raphson loop  ---
!
        IF( RSDMAX(NC).GT.1.D-6 ) THEN
          IF( NC.EQ.32 ) THEN
            RLXF = 6.D-1
            CYCLE
          ELSEIF( NC.EQ.200 ) THEN
            ECKE_ER = .TRUE.
          ELSE
            CYCLE
          ENDIF
        ENDIF
!
!---    Check mass balance for components
!
        DO NEQX=1,NEQC
          NSP =  ISP_S(NEQX+NEQE)
          IF( MOD(ISPLK(1),1000).NE.0 ) CYCLE
          DIFF = CX(NEQX)
          DO M=1,IEQ_C(1,NEQX)
            NSP=IEQ_C(M+1,NEQX)
            DIFF=DIFF-EQ_C(M,NEQX)*SP_CX(NSP)          
          ENDDO
          IF( CX(NEQX).EQ.0.D+0 .AND. DIFF.LE.1.D-8 ) CYCLE
          IF( DIFF .LE. 0.04*ABS(CX(NEQX)) ) CYCLE
          ECKE_ER = .TRUE.
        ENDDO
        EXIT L2
      ENDDO L2
      IF(IER.EQ.-1) ECKE_ER = .TRUE.
!
!---  Load activity coefficient for output  ---
!
      DO NSP = 1,NSPL
        ACTCX(NSP) = ACTVX(1,NSP)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ECKECN_R group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ECKEJCB( ACTVX,AJM,BJM,CX,COX,SP_CX,N,INDX,IER )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
! 
!     Numerical Recipes in Fortran 77, The Art of Scientific Computing
!     W.H. Press, B.P. Flannery, Saul A. Teukolsky, and W.T. Vetterling
!
!     INDX =  1 - This subroutine replaces the original funcv
!     INDX = -1 - This subroutine replaces the original fdjac
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 2 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AJM(LSPR,LSPR),BJM(LSPR)
      REAL*8 SP_CX(LSPR),DSP_CX(LSPR)
      REAL*8 CX(LEQC+LEQK),COX(LEQC+LEQK)
      REAL*8 DCLX(LSPL),CLX(LSPL),ACTVX(LSPL+1,LSPL)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ECKEJCB'
!
!---  Number of equations equals number of species  ---
!
      NEQR = NSPR
!
!---  Specie increments  ---
!
      DO NSP = 1,NSPR
        DSP_CX(NSP) = 1.D-5*SP_CX(NSP)
        IF( ABS(SP_CX(NSP)).LT.CMIN ) THEN
          IF( ISP_MN(NSP).EQ.1 ) THEN
            DSP_CX(NSP) = 1.D-6
          ELSE
            SP_CX(NSP) = CMIN
            DSP_CX(NSP) = 1.D-16
          ENDIF
        ENDIF
      ENDDO
!
!---  Load parameters for activity coefficient
!     calculations for aqueous species  ---
!
      DO NSP = 1,NSPL
        CLX(NSP) = SP_CX(NSP)
        DCLX(NSP) = DSP_CX(NSP)
      ENDDO
!
!---  Activity coefficients for aqueous species  ---
!
      IF( IACTV.EQ.3 ) THEN
        DO NSP = 1,NSPL
          DO M = 1,NSPL+1
            ACTVX(M,NSP) = ACTVC
          ENDDO
        ENDDO
      ELSEIF( IACTV.EQ.1 ) THEN
        CALL DAVIES( ACTVX,CLX,DCLX,SL(2,N),PORD(2,N),
     &    RHOL(2,N),T(2,N),XLW(2,N) )
      ELSEIF( IACTV.EQ.2 ) THEN
        CALL PITZER( ACTVX,CLX,DCLX,SL(2,N),PORD(2,N),
     &    RHOL(2,N),T(2,N),XLW(2,N) )
      ELSE
        CALL BDOT( ACTVX,CLX,DCLX,SL(2,N),PORD(2,N),
     &    RHOL(2,N),T(2,N),XLW(2,N) )
      ENDIF
!
!---  Initialize residual and residual partial
!     derivatives  ---
!
      DO NEQ = 1,NEQR
        BJM(NEQ) = 0.D+0
        DO M = 1,NEQR
          AJM(NEQ,M) = 0.D+0
        ENDDO
      ENDDO
!
!---  Loop over equations  ---
!
      IF( ISLC(57).EQ.0 ) THEN
        DO NEQ = 1,NEQR
!
!---      Equilibrium equation  ---
!
          IF( NEQ.LE.NEQE ) THEN
            CALL EECHEM( ACTVX,AJM,BJM,SP_CX,DSP_CX,N,NEQ,INDX )
!
!---      Conservation equation  ---
!
          ELSEIF( NEQ.LE.NEQE+NEQC ) THEN
            CALL CECHEM( ACTVX,AJM,BJM,CX,SP_CX,N,NEQ,INDX )
!
!---      Kinetic equation  ---
!
          ELSEIF( NEQ.LE.NEQE+NEQC+NEQK ) THEN
            CALL KECHEM( ACTVX,AJM,BJM,COX,SP_CX,DSP_CX,N,NEQ,INDX,IER )
          ENDIF
        ENDDO
      ELSE
        DO NEQ = 1,NEQR
!
!---      Equilibrium equation  ---
!         
          IF( NEQ.LE.NEQE ) THEN
            CYCLE
!         
!---      Conservation equation  ---
!         
          ELSEIF( NEQ.LE.NEQE+NEQC ) THEN
            CALL CECHEM_R( ACTVX,AJM,BJM,CX,SP_CX,N,NEQ,INDX )
!         
!---      Kinetic equation  ---
!         
          ELSEIF( NEQ.LE.NEQE+NEQC+NEQK ) THEN
            CALL KECHEM_R( ACTVX,AJM,BJM,COX,SP_CX,DSP_CX,N,NEQ,INDX )
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ECKEJCB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ELECTS(ZI,ZJ,IT,CPIX,APHI)
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
!     ECKEChemX
! 
!     This subroutine calculates higher order electrostatic functions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by A Felmy, from GMIN
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PTZRCOEF
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 A1(21),A2(21),A3(7),A4(9),A5(9),A6(9),P(10)
      REAL*8 J01,J02,J03,J11,J12,J13,J0,J1,J2,J21,J22,J23

      DATA A1/-.000000000010991D+0,-.000000000002563D+0
     & ,0.000000000001943D+0,0.000000000046333D+0,-.000000000050847D+0 
     & ,-.000000000821969D+0,0.000000001229405D+0,0.000000013522610D+0 
     & ,-.000000025267769D+0,-.000000202099617D+0,0.000000396566462D+0
     & ,0.000002937706971D+0,-.000004537895710D+0,-.000045036975204D+0
     & ,0.000036583601823D+0,0.000636874599598D+0,0.000388260636404D+0
     & ,-.007299499690937D+0,-.029779077456514D+0,-.060076477753119D+0
     & ,1.925154014814667D+0/

      DATA A2/0.000000000237816D+0,-.000000002849257D+0
     & ,-.000000006944757D+0,0.000000004558555D+0,0.000000080779570D+0
     & ,0.000000216991779D+0,-.000000250453880D+0,-.000003548684306D+0
     & ,-.000004583768938D+0,0.000034682122751D+0,0.000087294451594D+0
     & ,-.000242107641309D+0,-.000887171310131D+0,0.001130378079086D+0
     & ,0.006519840398744D+0,-.001668087945272D+0,-.036552745910311D+0
     & ,-.028796057604906D+0,0.150044637187895D+0,0.462762985338493D+0
     & ,0.628023320520852D+0/

      DATA A3/.000029308779366D+0,.000029869648486D+0
     & ,.000009838580455D+0,.000000827954226D+0,-.000000098522914D+0
     & ,.000000013943401D+0,-.000000005963131D+0/

      DATA A4/.018407992691D+0,.023626104695D+0,.005004621881D+0
     & ,-.000300844194D+0,-.000001533185D+0,.000009318246D+0
     & ,-.000004258797D+0,.000001509090D+0,-.000000766492D+0/

      DATA A5/3.9986000731D+0,3.7950588585D+0,-.3325673918D+0
     & ,-.0899335274D+0,.1338914658D+0,.1948882815D+0,.2368262404D+0
     & ,.1379406836D+0,.0478072558D+0/

      DATA A6/37.837805702D+0,24.470110234D+0,-3.719597578D+0
     & ,0.991895847D+0,-.327141729D+0,.121485594D+0,-.051539972D+0
     & ,.017779940D+0,-.011800766D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ELECTS'
      DSQI=DSQRT(CPIX)
      X1 = 6.D+0*APHI*DSQI
      P(1)=1.D+0
!
!--- Calculate integrals(see Harvie 1981)
!
      DO I = 1,3
        IF( I.EQ.1 ) XX = X1*ZI*ZI
        IF( I.EQ.2 ) XX = X1*ZJ*ZJ
        IF( I.EQ.3 ) XX = X1*ZI*ZJ
        BK = 0.D+0
        DK = 0.D+0
        BK1 = 0.D+0
        BK2 = 0.D+0
        DK1 = 0.D+0
        DK2 = 0.D+0
        IF( XX.LE.1.D+0 )THEN
          TT = 4.D+0*(XX**0.2D+0)
          DZ = TT/(5.D+0*XX)
          ZZ = TT-2.D+0
          IF( XX.LE.0.05D+0 )THEN
            Z2 = (XX/0.0245D+0)-1.040816D+0
            P(2) = Z2
            DO K = 3,7
              P(K)=2.D+0*Z2*P(K-1)-P(K-2)
            END DO
            J2 = 0.D+0
            DO K = 1,7
              J2 = J2+A3(K)*P(K)
            END DO
          ELSE
            Z2 = (XX/0.475D+0)-1.105263D+0
            P(2)=Z2
            DO K = 3,9
              P(K) = 2.D+0*Z2*P(K-1)-P(K-2)
            END DO
            J2 = 0.D+0
            DO K = 1,9
              J2 = J2+A4(K)*P(K)
            END DO
          END IF
          DO K = 1,21
            BK2 = BK1
            BK1 = BK
            BK = ZZ*BK1-BK2+A1(K)
            DK2 = DK1
            DK1 = DK
            DK = BK1+ZZ*DK1-DK2
          END DO
        ELSE
          TT = 4.444444444D+0*(XX**(-.1D+0))
          DZ = -TT/(10.D+0*XX)
          ZZ = TT-2.44444444D+0
          IF( XX.LE.50.D+0 )THEN
            Z2 = (XX/24.5D+0)-1.040816
            P(2) = Z2
            DO K = 3,9
              P(K) = 2.D+0*Z2*P(K-1)-P(K-2)
            END DO
            J2 = 0.D+0
            DO K = 1,9
              J2 = J2+A5(K)*P(K)
            END DO
        ELSE
            Z2 = (XX/425.D+0)-1.117647
            P(2) = Z2
            DO K = 3,9
              P(K) = 2.D+0*Z2*P(K-1)-P(k-2)
            END DO
            J2 = 0.D+0
            DO K = 1,9
              J2 = J2+A6(K)*P(K)
            END DO
          END IF
          DO K = 1,21
            BK2 = BK1
            BK1 = BK
            BK = ZZ*BK1-BK2+A2(K)
            DK2 = DK1
            DK1 = DK
            DK = BK1+ZZ*DK1-DK2
          END DO
        END IF
!
!---    Now calculate electrostatic functions
!
        J0 =0.25D+0*xx+0.5D+0*(BK-BK2)-1.D+0
        J1 = XX*(0.25D+0+0.5D+0*DZ*(DK-DK2))
        J2 = J2/XX
        J03 = J0
        J13 = J1
        J23 = J2
        IF( I.EQ.1 )THEN
          J01 = J0
          J11 = J1
          J21 = J2
        END IF
        IF( I.eq.2 )THEN
          J02 = J0
          J12 = J1
          J22 = J2
        END IF
      END DO
!
!---  Now calculate eth and ethp
!
      TMP = (ZI*ZJ)/(4.D+0*CPIX)
      ETH(IT) = TMP*(J03-0.5D+0*(J01+J02))
      ETHP(IT) = (TMP/(2.D+0*CPIX))*(J13-0.5D+0*(J11+J12))
     &  -ETH(IT)/CPIX
      ETHP2(IT) = -(1.D+0/(2.D+0*CPIX))*(ETH(IT)/CPIX+5.D+0*ETHP(IT))
     &  +(TMP/(4.D+0*CPIX*CPIX))*(J23-0.5D+0*(J21+J22))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ELECTS ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE EECHEM( ACTVX,AJM,BJM,SP_CX,DSP_CX,N,NEQ,INDX )
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
!     ECKEChemX
! 
!     Equilibrium Equation CHEMistry
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
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
      REAL*8 ACTVX(LSPL+1,LSPL)
      REAL*8 BJM(LSPR),AJM(LSPR,LSPR),SP_CX(LSPR),DSP_CX(LSPR)
      REAL*8 ACTEX(LESITE)
      LOGICAL FCHK
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/EECHEM'
!
!---  Skip for initial pH  ---
!
      IF( ISPLK(1).GT.1000 ) THEN
        IF( NEQ.EQ.IEQ_S(MOD(ISPLK(1),1000)) .AND.
     &    NSTEP.EQ.0 ) THEN
          AJM(NEQ,NEQ) = 1.D+0
          BJM(NEQ) = 0.D+0
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
      ENDIF
!
!--- Fixed species concentration or activity  ---
!
      DO NSLKX = 1,NSPLK
        IF( ISPLK(14+NSLKX).LT.0 ) THEN
          NSPX = ABS(ISPLK(14+NSLKX))
          IF( NSPX.GT.1000 ) NSPX = NSPX - 1000
          IF( NEQ.EQ.IEQ_S(NSPX) ) THEN
            AJM(NEQ,NEQ) = 1.D+0
            BJM(NEQ) = 0.D+0
            ISUB_LOG = ISUB_LOG-1
            RETURN
          ENDIF
        ENDIF
      ENDDO
!
!---  Total number of species  ---
!
      NSPR = NSPG + NSPL + NSPN + NSPS + NSPE
!
!---  Equilibrium constant as a function of temperature  ---
!
      NS = IEQ_E(1,NEQ)
      IRCE = IEQ_E((NS+2),NEQ)
      CALL EQCN( EQKX,T(2,N),IRCE,N )
!
!---  Volumetric concentration to molality  ---
!
      VTOMX = 1.D+0/(SL(2,N)*PORD(2,N)*RHOL(2,N)*XLW(2,N))
      VTOGX = 1.D+0/(MAX( 1.D-20,SG(2,N)*PORD(2,N)) )
      VTOLX = 1.D+0/(SL(2,N)*PORD(2,N))
      RHOWX = RHOL(2,N)*XLW(2,N)
!
!--- Prepare for exchanged species activity calculation
!
      IF( NSPE.NE.0 ) THEN
        DO ISITE = 1,NESITE
          ACTEX(ISITE) = 0.D+0
        ENDDO
!
!---    Gaines-Thomas convention
!
        IF( IACTEX.EQ.1 ) THEN
          DO NSP = 1,NSPE
            NCM_SP = NSPL + NSPS + NSP
            SP_CX(NCM_SP) = MAX(CMIN,SP_CX(NCM_SP))
            ISITE = ISP_E(NSP)
            ICAT = IEL_LK(NSP)
            ACTEX(ISITE) = ACTEX(ISITE)+SP_CX(NCM_SP)*SP_L(1,ICAT)
          ENDDO
        ENDIF
      ENDIF
!
!---  Base residual  ---
!
      IF( ISLC(60).EQ.0 ) THEN
      C_NEG = EQKX**EQ_E(NS,NEQ)
      ELSE
        C_NEG = EQ_E(NS,NEQ)*LOG(EQKX)
      ENDIF
!
!---  Loop over 1 species  ---
!
      DO M = 1,NS
        NCM_SP = IEQ_E(M+1,NEQ)
!
!---    Aqueous species,
!       concentration in molality, mol/kg H2O  ---
!
        IF( NCM_SP.LE.NSPL ) THEN
          CMX = SP_CX(NCM_SP)*VTOMX
          ACX = ACTVX(1,NCM_SP)
!
!---    Non-aqueous species,
!       concentration in kmolal, mol/kg H2O  ---
!
        ELSEIF(NCM_SP.LE.NSPL+NSPS) THEN
          CMX = SP_CX(NCM_SP)*VTOMX
          ACX = 1.D+0
!
!---    Exchanged species,
!       concentration in kmolal, mol/kg H2O  ---
!
        ELSEIF(NCM_SP.LE.NSPL+NSPS+NSPE) THEN
          ISITE = ISP_E(NCM_SP-NSPL-NSPS)
          CMX = SP_CX(NCM_SP)
          ICAT = IEL_LK(NCM_SP-NSPL-NSPS)
          ACX = CMX*SP_L(1,ICAT)/ACTEX(ISITE)
          CMX = 1.D+0
!
!---    Gas species,
!       convert from mol/node volume to mol/volume of gas,
!       use RHOWX for Henry's law  ---
!
        ELSEIF(NCM_SP.LE.NSPL+NSPS+NSPE+NSPG) THEN
          CMX = SP_CX(NCM_SP)*VTOGX/RHOWX
          ACX = 1.D+0
        ENDIF
!
!---    Equilibrium species  ---
!
        IF( M.EQ.1 ) THEN
          C_POS = ACX*CMX
!
!---    Conservation species  ---
!
        ELSE
          IF( ISLC(60).EQ.0 ) THEN
            C_NEG = C_NEG*((ACX*CMX)**EQ_E(M-1,NEQ))
          ELSE
            C_NEG = C_NEG+EQ_E(M-1,NEQ)*LOG(ACX*CMX)
          ENDIF
        ENDIF
      ENDDO
      IF( ISLC(60).EQ.1 ) THEN
        C_NEG = EXP(C_NEG)
      ENDIF
      BJM(NEQ) = C_POS - C_NEG
!
!---  Return residual vector  ---
!
      IF( INDX.EQ.1 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Incremented residuals, computed numerically  ---
!
      DO NSP = 1,NSPR
!
!---    Check whether specie is a 1 equation specie,
!         which affects the residual via the 1 equation
!       or an aqueous specie,
!         which affects the residual via the 1 equation,
!         via the 1 equation reactions
!         via the 1 constant,
!         via the activity coefficient,
!         via the ionic strength  ---
!
        FCHK = .FALSE.
        IF( NSP.LE.NSPL ) FCHK = .TRUE.
!
!---    Loop over 1 equation species  ---
!
        IF( .NOT.FCHK ) THEN
          DO M = 1,NS
            NSP_E = IEQ_E(M+1,NEQ)
            IF( NSP.EQ.NSP_E ) THEN
              FCHK = .TRUE.
              EXIT
            ENDIF
          ENDDO
        ENDIF
!
!---    Skip for initial pH  ---
!
        IF( ISPLK(1).GT.1000 ) THEN
          IF( NSP.EQ.MOD(ISPLK(1),1000) .AND.
     &      NSTEP.EQ.0 ) FCHK = .FALSE.
        ENDIF
!
!---    Specie is either an 1 equation specie
!       or an aqueous specie  ---
!
        IF( FCHK ) THEN
!
!---      Base residual  ---
!
          IF( ISLC(60).EQ.0 ) THEN
          C_NEG = EQKX**EQ_E(NS,NEQ)
          ELSE
            C_NEG = EQ_E(NS,NEQ)*LOG(EQKX)
          ENDIF
!
!---      Loop over 1 species  ---
!
          DO M = 1,NS
            NCM_SP = IEQ_E(M+1,NEQ)
!
!---        Aqueous species,
!           concentration in molality, mol/kg H2O  ---
!
            IF( NCM_SP.LE.NSPL ) THEN
!
!---          Incremented species  ---
!
              IF( NCM_SP.EQ.NSP )  THEN
                CMX = (SP_CX(NCM_SP)+DSP_CX(NCM_SP))*VTOMX
                ACX = ACTVX(NSP+1,NCM_SP)
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NCM_SP)*VTOMX
                ACX = ACTVX(1,NCM_SP)
              ENDIF
!
!---        Non-aqueous species,
!           concentration in kmolal, mol/kg H2O  ---
!
            ELSEIF(NCM_SP.LE.NSPL+NSPS) THEN
!
!---          Incremented species  ---
!
              IF( NCM_SP.EQ.NSP )  THEN
                CMX = (SP_CX(NCM_SP)+DSP_CX(NCM_SP))*VTOMX
                ACX = 1.D+0
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NCM_SP)*VTOMX
                ACX = 1.D+0
              ENDIF
!
!---        Exchanged species  ---
!
            ELSEIF(NCM_SP.LE.NSPL+NSPS+NSPE) THEN
!
!---          Incremented species  ---
!
              ISITE = ISP_E(NCM_SP-NSPL-NSPS)
              ICAT = IEL_LK(NCM_SP-NSPL-NSPS)
              CHARGX = SP_L(1,ICAT)
              IF( NCM_SP.EQ.NSP )  THEN
                CMX = SP_CX(NCM_SP)+DSP_CX(NCM_SP)
                ACX = CMX*CHARGX/(ACTEX(ISITE)+DSP_CX(NCM_SP)*CHARGX)
                CMX = 1.D+0
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NCM_SP)
                ACX = CMX*CHARGX/ACTEX(ISITE)
                CMX = 1.D+0
              ENDIF
!
!---        Gas species
!
            ELSEIF(NCM_SP.LE.NSPL+NSPS+NSPE+NSPG) THEN
!
!---          Incremented species  ---
!
              IF( NCM_SP.EQ.NSP )  THEN
                CMX = (SP_CX(NCM_SP)+DSP_CX(NCM_SP))*VTOGX/RHOWX
                ACX = 1.D+0
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NCM_SP)*VTOGX/RHOWX
                ACX = 1.D+0
              ENDIF
            ENDIF
!
!---        Equilibrium species  ---
!
            IF( M.EQ.1 ) THEN
              C_POS = ACX*CMX
!
!---        Conservation species  ---
!
            ELSE
              IF( ISLC(60).EQ.0 ) THEN
                C_NEG = C_NEG*((ACX*CMX)**EQ_E(M-1,NEQ))
              ELSE
                C_NEG = C_NEG+EQ_E(M-1,NEQ)*LOG(ACX*CMX)
              ENDIF
            ENDIF
          ENDDO
          IF( ISLC(60).EQ.1 ) THEN
            C_NEG = EXP(C_NEG)
          ENDIF
          AVX = C_NEG + BJM(NEQ)
          BVX = ABS(C_NEG) + ABS(BJM(NEQ))
          CVX = C_POS
          IF( BVX/EPSL.GT.EPSL ) THEN
            IF( ABS(AVX)/BVX.GT.EPSL ) THEN
              CVX = C_POS - AVX
            ENDIF
          ENDIF
          AJM(NEQ,IEQ_S(NSP)) = CVX/DSP_CX(NSP)
          IF( ABS(AJM(NEQ,IEQ_S(NSP))).LT.1.D-20 ) 
     &      AJM(NEQ,IEQ_S(NSP)) = 0.D+0
        ENDIF
      ENDDO
!
!---  Return Jacobian matrix and residual vector  ---
!
      IF( ABS(INDX).EQ.1 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      BJM(NEQ) = -BJM(NEQ)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of EECHEM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE EQCN( EQKX,TX,INDX,N )
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
!     ECKEChemX
! 
!     Equilibrium reaction constant as a function of temperature.
!
!     INDX > 0 : 1 reaction index  RC_E
!     INDX < 0 : kinetic reaction index RC_K
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
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
      SUB_LOG(ISUB_LOG) = '/EQCN'
!
!---  Absolute temperature  ---
!
      TKX = TX+TABS
!
!---  Equilibrium reaction  ---
!
      IF( INDX.GT.0 ) THEN
        EQKX = RC_E(1,INDX)*LOG(TKX) + RC_E(2,INDX) +
     &    RC_E(3,INDX)*TKX + RC_E(4,INDX)/TKX +
     &    RC_E(5,INDX)/(TKX**2)
        EQKX = EXP(TOLN*EQKX)
!
!---  Kinetic reaction  ---
!
      ELSEIF( INDX.LT.0 ) THEN
        JNDX = -INDX
        NSPRX = IRC_K(1,JNDX)
        NSPPX = IRC_K(2,JNDX)
        NSPKX = NSPRX+NSPPX
        N4 = MAX(1,N*IRCKN(NSPKX+4))
        N5 = MAX(1,N*IRCKN(NSPKX+5))
        N6 = MAX(1,N*IRCKN(NSPKX+6))
        N7 = MAX(1,N*IRCKN(NSPKX+7))
        N8 = MAX(1,N*IRCKN(NSPKX+8))
        EQKX = RC_K(NSPKX+4,N4,JNDX)*LOG(TKX) + RC_K(NSPKX+5,N5,JNDX) +
     &    RC_K(NSPKX+6,N6,JNDX)*TKX + RC_K(NSPKX+7,N7,JNDX)/TKX +
     &    RC_K(NSPKX+8,N8,JNDX)/(TKX**2)
        EQKX = EXP(TOLN*EQKX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of EQCN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE EQEQ( EQKX,SP_CX,VTOMX )
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
!     ECKEChemX
! 
!     Calculate 1 species concentrations according to
!     the conservation species concentrations, assuming
!     an activity of 1.0
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 SP_CX(LSPR)
      REAL*8 EQKX(LEQE)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/EQEQ'
!
!---  Set 1 species concentrations according to
!     the conservation species concentrations, assuming
!     an activity of 1.0  ---
!
      DO NEQ = 1,NEQE
        NEQ_SP = IEQ_E(2,NEQ)
        NS = IEQ_E(1,NEQ)
!
!---    Equilibrium constant and exponent  ---
!
        IF( ISLC(60).EQ.0 ) THEN
          C_NEG = EQKX(NEQ)**EQ_E(NS,NEQ)
        ELSE
         C_NEG = EQ_E(NS,NEQ)*LOG(EQKX(NEQ))
        ENDIF
!
!---    Loop over species in 1 equation
!       skipping the 1 species  ---
!
        DO K = 1,NS
          NCM_SP = IEQ_E(K+1,NEQ)
!
!---      Skip 1 species  ---
!
          IF( NCM_SP.EQ.NEQ_SP ) CYCLE
!
!---      Convert conservation species concentration to
!         molality, mol/kg H2O  ---
!
          SP_CX(NCM_SP) = MAX(CMIN,SP_CX(NCM_SP))
          CMX = SP_CX(NCM_SP)*VTOMX
          IF( ISLC(60).EQ.0 ) THEN
            C_NEG = C_NEG*(CMX**EQ_E(K-1,NEQ))
          ELSE
            C_NEG = C_NEG+EQ_E(K-1,NEQ)*LOG(CMX)
          ENDIF
        ENDDO
        IF( ISLC(60).EQ.1 ) THEN
          C_NEG = EXP(C_NEG)
        ENDIF
!
!---    Convert 1 species concentration to
!       node volumetric, mol/m^3  ---
!
        SP_CX(NEQ_SP) = C_NEG/VTOMX
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of EQEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE GAUSSJ( AX,BX,NX,NPX,IERX )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
! 
!
!     Numerical Recipes in Fortran 77, The Art of Scientific Computing
!     W.H. Press, B.P. Flannery, Saul A. Teukolsky, and W.T. Vetterling
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 3 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AX(NPX,NPX),BX(NPX)
      INTEGER IPIV(NPX),INDXR(NPX),INDXC(NPX)
!
!----------------------Executable Lines--------------------------------!
!
      IERX = 0
      EPSL = 1.D-14
      DO J = 1,NX
        IPIV(J) = 0
      ENDDO
      DO I = 1,NX
        BIGX = 0.D+0
        DO J = 1,NX
          IF( IPIV(J).NE.1 ) THEN
            DO K = 1,NX
              IF( IPIV(K).EQ.0 ) THEN
                IF( ABS(AX(J,K)).GE.BIGX ) THEN
                  BIGX = ABS(AX(J,K))
                  IROW = J
                  ICOL = K
                ENDIF
              ELSEIF( IPIV(K).GT.1 ) THEN
                M_ERR(1) = 'Execution Error: GAUSSJ: Pivot Error: '
                IF( N_DB.LT.0 ) THEN
                  M_ERR(2) = ' at Boundary Surface: '
                ELSE
                  M_ERR(2) = ' at Node: '
                ENDIF
                CALL PATH
                I_ERR(1) = ABS(N_DB)
                I_ERR(2) = 1
                I_ERR(3) = 2
                I_ERR(4) = ID
                IERX = -1
                RETURN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        IPIV(ICOL) = IPIV(ICOL) + 1
        IF( IROW.NE.ICOL ) THEN
          DO L = 1,NX
            DUMX = AX(IROW,L)
            AX(IROW,L) = AX(ICOL,L)
            AX(ICOL,L) = DUMX
          ENDDO
          DUMX = BX(IROW)
          BX(IROW) = BX(ICOL)
          BX(ICOL) = DUMX
        ENDIF
        INDXR(I) = IROW
        INDXC(I) = ICOL
        IF( ABS(AX(ICOL,ICOL))/EPSL.LT.EPSL ) THEN
          M_ERR(1) = 'Execution Error: GAUSSJ: Singular Matrix: '
          IF( N_DB.LT.0 ) THEN
            M_ERR(2) = ' at Boundary Surface: '
          ELSE
            M_ERR(2) = ' at Node: '
          ENDIF
          CALL PATH
          I_ERR(1) = ABS(N_DB)
          I_ERR(2) = 1
          I_ERR(3) = 2
          I_ERR(4) = ID
          IERX = -1
          RETURN
        ENDIF
        PIVINV = 1.D+0/AX(ICOL,ICOL)
        AX(ICOL,ICOL) = 1.D+0
        DO L = 1,NX
          AX(ICOL,L) = AX(ICOL,L)*PIVINV
        ENDDO
        BX(ICOL) = BX(ICOL)*PIVINV
        DO LL = 1,NX
          IF( LL.NE.ICOL ) THEN
            DUMX = AX(LL,ICOL)
            AX(LL,ICOL) = 0.D+0
            DO L = 1,NX
              AX(LL,L) = AX(LL,L) - AX(ICOL,L)*DUMX
            ENDDO
            BX(LL) = BX(LL) - BX(ICOL)*DUMX
          ENDIF
        ENDDO
      ENDDO
!
!---  End of GAUSSJ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      FUNCTION GETSPNM( INDX )
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
!     ECKEChemX
!
!     Get species name from global species index.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 1 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 GETSPNM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/GETSPNM'
!
!---  Aqueous species  ---
!
      IF( INDX.LE.NSPL ) THEN
        GETSPNM = SPNML(INDX)
!
!---  Solid species  ---
!
      ELSEIF( INDX.LE.(NSPL+NSPS) ) THEN
        GETSPNM = SPNMS(INDX-NSPL)
!
!---  Exchange species  ---
!
      ELSEIF( INDX.LE.(NSPL+NSPS+NSPE) ) THEN
        GETSPNM = SPNME(INDX-NSPL-NSPS)
!
!---  Gas species  ---
!
      ELSEIF( INDX.LE.(NSPL+NSPS+NSPE+NSPG) ) THEN
        GETSPNM = SPNMG(INDX-NSPL-NSPS-NSPE)







      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of GETSPNM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE GUSPCN( CX,SP_CX,N )
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
!     ECKEChemX
!
!     Guess specie concentrations using component rule of thumb,
!     1 equations, and pH.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
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
      REAL*8 CX(LEQC+LEQK),CIX(LEQC+LEQK)
      REAL*8 SP_CX(LSPR),SP_CIX(LSPR)
      REAL*8 DERRX(2,LEQC),EQKX(LEQE)
      INTEGER INEG(LEQC),IPOS(LEQC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/GUSPCN'
!
!---  Total number of equations and species  ---
!
      NSPR = NSPG + NSPL + NSPN + NSPS + NSPE
!
!---  Volumetric concentration to molality  ---
!
      VTOMX = 1.D+0/(SL(2,N)*PORD(2,N)*RHOL(2,N)*XLW(2,N))
!
!---  Set initial concentrations  ---
!
      DO NSP = 1,NSPR
        SP_CIX(NSP) = SP_CX(NSP)
      ENDDO
!
!---  Set conservation species concentrations to 1/50 of
!     the total component concentration  ---
!
      NSPH = MOD(ISPLK(1),1000)
      DO NEQ = 1,NEQC
        NSP = IEQ_C(2,NEQ)
        IF( NSP.EQ.NSPH ) CYCLE
        SP_CIX(NSP) = MAX( CX(NEQ)*5.D-2,0.D+0 )
      ENDDO
!
!---  If uninitialized set pH to 1.D-7 molality, mol/kg water  ---
!
      IF( NSPH.GE.1 .AND. NSPH.LE.LSPR ) THEN
        IF( SP_CX(NSPH)/EPSL.LT.EPSL ) SP_CIX(NSPH) = 1.D-7/VTOMX
      ENDIF
!
!---  Equilibrium constant as a function of temperature  ---
!
      DO NEQ = 1,NEQE
        NS = IEQ_E(1,NEQ)
        IRCE = IEQ_E((NS+2),NEQ)
        CALL EQCN( EQKX(NEQ),T(2,N),IRCE,N )
      ENDDO
!
!---  Recalculate 1 species concentrations according to
!     the conservation species concentrations, assuming
!     an activity of 1.0  ---
!
      CALL EQEQ( EQKX,SP_CIX,VTOMX )
!
!---  Total error in conservation species concentration  ---
!
      DO NEQ = 1,NEQC
        DERRX(1,NEQ) = 0.D+0
        CIX(NEQ) = 0.D+0
        DO M = 1,IEQ_C(1,NEQ)
          NSP = IEQ_C(M+1,NEQ)
          CIX(NEQ) = CIX(NEQ)+EQ_C(M,NEQ)*SP_CIX(NSP)
        ENDDO
        DERRX(1,NEQ) = ABS(CIX(NEQ)-CX(NEQ))/(CX(NEQ)+SMALL)
      ENDDO
!
!---  Adjust conservation species concentrations by 1.667 and search
!     for improvements in total error in component species
!     concentrations  ---
!
      NC = 0
      L1: DO
        NC = NC + 1
        IF( NC.GT.128 ) EXIT L1
!
!---    Loop over conservation species  ---
!
        L2: DO NEQX = 1,NEQC
          NSPX = IEQ_C(2,NEQX)
!
!---      Initialize conservation species direction flags  ---
!
          IPOS(NEQX) = 0
          INEG(NEQX) = 0
!
!---      Skip for initial pH  ---
!
          IF( NSPX.EQ.NSPH .AND. ISPLK(1).GT.1000 ) CYCLE L2
!
!---      Decrease conservation species concentration by 1.667  ---
!
          SP_CIX(NSPX) = SP_CIX(NSPX)/1.667D+0
!
!---      Recalculate 1 species concentrations according to
!         the conservation species concentrations, assuming
!         an activity of 1.0  ---
!
          CALL EQEQ( EQKX,SP_CIX,VTOMX )
!
!---      Total error in conservation species concentration  ---
!
          DO NEQ = 1,NEQC
            DERRX(2,NEQ) = 0.D+0
            CIX(NEQ) = 0.D+0
            DO M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              CIX(NEQ) = CIX(NEQ)+EQ_C(M,NEQ)*SP_CIX(NSP)
            ENDDO
            DERRX(2,NEQ) = ABS(CIX(NEQ)-CX(NEQ))/(CX(NEQ)+SMALL)
          ENDDO
!
!---      Check for improvement in total error in component species
!         concentration  ---
!
          CERRX = (DERRX(1,NEQX)-DERRX(2,NEQX))/(DERRX(1,NEQX)+SMALL)
!
!---      Improvement found  ---
!
          IF( CERRX.GT.1.D-2 ) THEN
            DO NEQ = 1,NEQC
              DERRX(1,NEQ) = DERRX(2,NEQ)
            ENDDO
            INEG(NEQX) = 1
!
!---      No improvement found, increase conservation species
!         concentration ---
!
          ELSE
            SP_CIX(NSPX) = SP_CIX(NSPX)*(1.667D+0**2)
!
!---        Recalculate 1 species concentrations according to
!           the conservation species concentrations, assuming
!           an activity of 1.0  ---
!
            CALL EQEQ( EQKX,SP_CIX,VTOMX )
!
!---        Total error in component species concentration  ---
!
            DO NEQ = 1,NEQC
              DERRX(2,NEQ) = 0.D+0
              CIX(NEQ) = 0.D+0
              DO M = 1,IEQ_C(1,NEQ)
                NSP = IEQ_C(M+1,NEQ)
                CIX(NEQ) = CIX(NEQ)+EQ_C(M,NEQ)*SP_CIX(NSP)
              ENDDO
              DERRX(2,NEQ) = ABS(CIX(NEQ)-CX(NEQ))/(CX(NEQ)+SMALL)
            ENDDO
!
!---        Check for improvement in total error in component species
!           concentration  ---
!
            CERRX = (DERRX(1,NEQX)-DERRX(2,NEQX))/(DERRX(1,NEQX)+SMALL)
!
!---        Improvement found  ---
!
            IF( CERRX.GT.1.D-2 ) THEN
              DO NEQ = 1,NEQC
                DERRX(1,NEQ) = DERRX(2,NEQ)
              ENDDO
              IPOS(NEQX) = 1
!
!---        No improvement found, reset conservation species
!           concentration ---
!
            ELSE
              SP_CIX(NSPX) = SP_CIX(NSPX)/1.667D+0
            ENDIF
          ENDIF
        ENDDO L2
!
!---    Loop over conservation species, checking for no
!       further improvements  ---
!
        DO NEQ = 1,NEQC
          IF( IPOS(NEQ).EQ.1 .OR. INEG(NEQ).EQ.1 ) CYCLE L1
        ENDDO
        EXIT L1
      ENDDO L1
!
!---  Set species concentration guesses  ---
!
      DO NSP = 1,NSPR
        SP_CX(NSP) = MAX( SP_CIX(NSP),0.D+0 )
      ENDDO
!
!---  End of GUSPCN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLHSP
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
!     ECKEChemX
! 
!     Convert reactive species initial condition concentrations into
!     mol/m^3 node volume
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNNL, 1 February 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
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
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 GETSPNM
      EXTERNAL GETSPNM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLHSP'
      NSPR = NSPG + NSPL + NSPN + NSPS + NSPE
!
!---  Define completely immobile conservation component species  ---
!
      DO NEQ = 1,NEQC
        IMMB(NEQ) = 1
!
!---    Loop over conservation-component species  ---
!
        DO M = 1,IEQ_C(1,NEQ)
          NSP = IEQ_C(M+1,NEQ)
!
!---      Mobile species ---
!
          IF( NSP.LE.NSPL .OR. NSP.GT.(NSPL+NSPS+NSPE) ) IMMB(NEQ) = 0
        ENDDO
      ENDDO
!
!---  Define completely immobile kinetic component species  ---
!
      DO NEQ = 1,NEQK
        IMMB(NEQ+NEQC) = 1
!
!---    Loop over kinetic-component species  ---
!
        DO M = 1,IEQ_K(1,NEQ)
          NSP = IEQ_K(M+1,NEQ)
!
!---      Mobile species ---
!
          IF( NSP.LE.NSPL .OR. NSP.GT.(NSPL+NSPS+NSPE) )
     &      IMMB(NEQ+NEQC) = 0
        ENDDO
      ENDDO
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = ND(N)
!
!---    Loop over reactive species  ---
!
        DO NSP = 1,NSPR
          IF( INDEX(GETSPNM(NSP),'fix').NE.0 ) THEN
            SP_CI(N,NSP) = SP_C(N,NSP)
          ENDIF
!
!---      Restart simulation  ---
!
          IF( IEO.EQ.2 .AND. ISLC(76).EQ.0 ) THEN
!
!---        Convert reactive species from node volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            IF( IC_SP(N,NSP).EQ.11 ) THEN
              SP_C(N,NSP) = 1.D+3*SP_C(N,NSP)
!
!---        Convert reactive species from aqueous volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( IC_SP(N,NSP).EQ.12 ) THEN
              SP_C(N,NSP) = 1.D+3*SP_C(N,NSP)*SL(2,N)*PORD(2,N)
!
!---        Convert reactive species from aqueous molal, kmol/kg water
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( IC_SP(N,NSP).EQ.13 ) THEN
              SP_C(N,NSP) = 1.D+3*SP_C(N,NSP)*SL(2,N)*PORD(2,N)*
     &          RHOL(2,N)*XLW(2,N)
!
!---        Convert reactive species from gas volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( IC_SP(N,NSP).EQ.14 ) THEN
              SP_C(N,NSP) = 1.D+3*SP_C(N,NSP)*SG(2,N)*PORD(2,N)
!
!---        Convert solid-species mineral volumetric fraction
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
              DO NRC = 1,NRCK
                IF( IRCKT(NRC).EQ.20 ) THEN
                  NSPX = IRC_K(1,NRC)
                ELSE
                  NSPRX = IRC_K(1,NRC)
                  NSPPX = IRC_K(2,NRC)
                  NSPX = IRC_K(3+NSPRX+NSPPX,NRC)
                ENDIF
                IF( NSP.EQ.NSPX ) THEN
                  NSPX = NSP-NSPL
!
!---              Restart simulation w/ lithology overwrite  ---
!
                  IF( SP_S(2,NSPX)/EPSL.GT.EPSL .AND.
     &              ISP_OW(NSPX,N).EQ.1 ) THEN
                    SP_CMN(N,NSPX) = 1.D+3*RS_S(2,NSPX,N)*
     &                SP_S(1,NSPX)/SP_S(2,NSPX)
                  ENDIF
                  EXIT
                ENDIF
              ENDDO
            ENDIF
!
!---      Normal simulation  ---
!
          ELSE
!
!---        Convert reactive species from node volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            IF( IC_SP(N,NSP).EQ.1 ) THEN
              SP_C(N,NSP) = 1.D+3*SP_C(N,NSP)
!
!---        Convert reactive species from aqueous volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( IC_SP(N,NSP).EQ.2 ) THEN
              SP_C(N,NSP) = 1.D+3*SP_C(N,NSP)*SL(2,N)*PORD(2,N)
!
!---        Convert reactive species from aqueous molal, kmol/kg water
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( IC_SP(N,NSP).EQ.3 ) THEN
              SP_C(N,NSP) = 1.D+3*SP_C(N,NSP)*SL(2,N)*PORD(2,N)*
     &          RHOL(2,N)*XLW(2,N)
!
!---        Convert reactive species from gas volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( IC_SP(N,NSP).EQ.4 ) THEN
              SP_C(N,NSP) = 1.D+3*SP_C(N,NSP)*SG(2,N)*PORD(2,N)
!
!---        Convert solid-species mineral volumetric fraction
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
              DO NRC = 1,NRCK
                IF( IRCKT(NRC).EQ.20 ) THEN
                  NSPX = IRC_K(1,NRC)
                ELSE
                  NSPRX = IRC_K(1,NRC)
                  NSPPX = IRC_K(2,NRC)
                  NSPX = IRC_K(3+NSPRX+NSPPX,NRC)
                ENDIF
                IF( NSP.EQ.NSPX ) THEN
                  NSPX = NSP-NSPL
                  IF( SP_S(2,NSPX)/EPSL.GT.EPSL ) THEN
                    SP_CMN(N,NSPX) = 1.D+3*RS_S(2,NSPX,N)*
     &                SP_S(1,NSPX)/SP_S(2,NSPX)
                  ENDIF
                  EXIT
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
!
!---    Component species concentration  ---
!
        DO NEQ = 1,NEQC
          NSL = NEQ + NSOLU
!
!---      Restart simulation  ---
!
          IF( IEO.EQ.2 .AND. ISLC(76).EQ.0 ) THEN
!
!---        Convert reactive species from node volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            IF( ICT(N,NSL).EQ.11 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)
!
!---        Convert reactive species from aqueous volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.12 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SL(2,N)*PORD(2,N)
!
!---        Convert reactive species from aqueous molal, kmol/kg water
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.13 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SL(2,N)*PORD(2,N)*
     &          RHOL(2,N)*XLW(2,N)
!
!---        Convert reactive species from gas volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.14 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SG(2,N)*PORD(2,N)
            ELSE
              IEDL(NSL) = ISP_IEDL
              IF( IEQW.GT.0 ) THEN
                SMDL(NSL) = SP_MDL
                IF( ISP_IEDL.EQ.2 .OR. ISP_IEDL.EQ.4 ) THEN
                  SDCL(1,N,NSL) = SP_SDCL(1)
                  SDCL(2,N,NSL) = SP_SDCL(2)
                  SDCL(3,N,NSL) = SP_SDCL(3)
                ENDIF
              ENDIF
              IF( IEQA.GT.0 .OR. IOM.EQ.30 .OR. IOM.EQ.40 .OR. 
     &          IOM.EQ.43 ) SMDG(NSL) = SP_MDG



              C(N,NSL) = 0.D+0
              DO NSP = 1,IEQ_C(1,NEQ)
                C(N,NSL) = C(N,NSL)+EQ_C(NSP,NEQ)*
     &            SP_C(N,IEQ_C(NSP+1,NEQ))
              ENDDO
            ENDIF
            CO(N,NSL) = C(N,NSL)
!
!---      Normal simulation  ---
!
          ELSE
!
!---        Convert reactive species from node volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            IF( ICT(N,NSL).EQ.1 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)
!
!---        Convert reactive species from aqueous volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.2 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SL(2,N)*PORD(2,N)
!
!---        Convert reactive species from aqueous molal, kmol/kg water
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.3 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SL(2,N)*PORD(2,N)*
     &          RHOL(2,N)*XLW(2,N)
!
!---        Convert reactive species from gas volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.4 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SG(2,N)*PORD(2,N)
            ELSE
              IEDL(NSL) = ISP_IEDL
              IF( IEQW.GT.0 ) THEN
                SMDL(NSL) = SP_MDL
                IF( ISP_IEDL.EQ.2 .OR. ISP_IEDL.EQ.4 ) THEN
                  SDCL(1,N,NSL) = SP_SDCL(1)
                  SDCL(2,N,NSL) = SP_SDCL(2)
                  SDCL(3,N,NSL) = SP_SDCL(3)
                ENDIF
              ENDIF
              IF( IEQA.GT.0 .OR. IOM.EQ.30 .OR. IOM.EQ.40 .OR. 
     &          IOM.EQ.43 ) SMDG(NSL) = SP_MDG



              C(N,NSL) = 0.D+0
              DO NSP = 1,IEQ_C(1,NEQ)
                C(N,NSL) = C(N,NSL)+EQ_C(NSP,NEQ)*
     &            SP_C(N,IEQ_C(NSP+1,NEQ))
              ENDDO
            ENDIF
            CO(N,NSL) = C(N,NSL)
          ENDIF
        ENDDO
!
!---    Kinetic species concentration  ---
!
        DO NEQ = 1,NEQK
          NSL = NEQ + NEQC + NSOLU
!
!---      Restart simulation  ---
!
          IF( IEO.EQ.2 .AND. ISLC(76).EQ.0 ) THEN
!
!---        Convert reactive species from node volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            IF( ICT(N,NSL).EQ.11 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)
!
!---        Convert reactive species from aqueous volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.12 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SL(2,N)*PORD(2,N)
!
!---        Convert reactive species from aqueous molal, kmol/kg water
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.13 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SL(2,N)*PORD(2,N)*
     &          RHOL(2,N)*XLW(2,N)
!
!---        Convert reactive species from gas volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.14 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SG(2,N)*PORD(2,N)
            ELSE
              IEDL(NSL) = ISP_IEDL
              IF( IEQW.GT.0 ) THEN
                SMDL(NSL) = SP_MDL
                IF( ISP_IEDL.EQ.2 .OR. ISP_IEDL.EQ.4 ) THEN
                  SDCL(1,N,NSL) = SP_SDCL(1)
                  SDCL(2,N,NSL) = SP_SDCL(2)
                  SDCL(3,N,NSL) = SP_SDCL(3)
                ENDIF
              ENDIF
              IF( IEQA.GT.0 .OR. IOM.EQ.30 .OR. IOM.EQ.40 .OR. 
     &          IOM.EQ.43 ) SMDG(NSL) = SP_MDG



              C(N,NSL) = 0.D+0
              DO NSP = 1,IEQ_K(1,NEQ)
                C(N,NSL) = C(N,NSL)+EQ_K(NSP,NEQ)*
     &            SP_C(N,IEQ_K(NSP+1,NEQ))
              ENDDO
              CO(N,NSL) = C(N,NSL)
            ENDIF
!
!---      Normal simulation  ---
!
          ELSE
!
!---        Convert reactive species from node volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            IF( ICT(N,NSL).EQ.1 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)
!
!---        Convert reactive species from aqueous volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.2 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SL(2,N)*PORD(2,N)
!
!---        Convert reactive species from aqueous molal, kmol/kg water
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.3 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SL(2,N)*PORD(2,N)*
     &          RHOL(2,N)*XLW(2,N)
!
!---        Convert reactive species from gas volumetric, kmol/m^3
!           to node volumetric, mol/m^3  ---
!
            ELSEIF( ICT(N,NSL).EQ.4 ) THEN
              C(N,NSL) = 1.D+3*C(N,NSL)*SG(2,N)*PORD(2,N)
            ELSE
              IEDL(NSL) = ISP_IEDL
              IF( IEQW.GT.0 ) THEN
                SMDL(NSL) = SP_MDL
                IF( ISP_IEDL.EQ.2 .OR. ISP_IEDL.EQ.4 ) THEN
                  SDCL(1,N,NSL) = SP_SDCL(1)
                  SDCL(2,N,NSL) = SP_SDCL(2)
                  SDCL(3,N,NSL) = SP_SDCL(3)
                ENDIF
              ENDIF
              IF( IEQA.GT.0 .OR. IOM.EQ.30 .OR. IOM.EQ.40 .OR. 
     &          IOM.EQ.43 ) SMDG(NSL) = SP_MDG



              C(N,NSL) = 0.D+0
              DO NSP = 1,IEQ_K(1,NEQ)
                C(N,NSL) = C(N,NSL)+EQ_K(NSP,NEQ)*
     &            SP_C(N,IEQ_K(NSP+1,NEQ))
              ENDDO
              CO(N,NSL) = C(N,NSL)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = ND(N)
!
!---    Loop over reactive species  ---
!
        DO NSP = 1,NSPR
          SP_CO(N,NSP) = SP_C(N,NSP)
        ENDDO
      ENDDO
!
!---  Assign boundary solute concentrations for initial condition
!     type boundary conditions  ---
!
      DO NSP = 1,NSPR
        DO NB = 1,NBC(ID+1)
          N = IBCN(NB)
          SP_CBO(NB,NSP) = SP_C(N,NSP)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLHSP group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HOMIX(CPIX,APHI)
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
!     ECKEChemX
! 
!     This subroutine calculates higher order mixing terms for
!     unsymmetrical electrolyte mixings.
!
!----------------------Authors-----------------------------------------!
!
!     Written by A. Felmy, from GMIN
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
      USE PTZRCOEF
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
      SUB_LOG(ISUB_LOG) = '/HOMIX'
      MAXETH = LMCG**2
      DO I = 1,MAXETH
        ETH(I)=0.D+0
        ETHP(I)=0.D+0
        ETHP2(I)=0.D+0
      ENDDO
      IF( NCC_PZ.GE.2 )THEN
        NT = 1
        DO I = 1,NCC_PZ
          DO J = I+1,NCC_PZ
            IT = 1
            ZI = SP_L(1,JPC(I))
            ZJ = SP_L(1,JPC(J))
            IF( ZI.NE.ZJ ) THEN
              IT = INT(ZI*ZJ)
              IF( ETH(IT).EQ.0.D+0 ) CALL ELECTS(ZI,ZJ,IT,CPIX,APHI)
            ENDIF
            CTCPH(NT) = TCC(NT)+ETH(IT)+CPIX*ETHP(IT)
            CTC(NT) = TCC(NT)+ETH(IT)
            CTCPR(NT) = ETHP(IT)
            CTCPPR(NT) = ETHP2(IT)
            NT = NT+1
          ENDDO
        ENDDO
      ENDIF
      IF( NA_PZ.GE.2 )THEN
        NT = 1
        DO I = 1,NA_PZ
          DO J = I+1,NA_PZ
            IT = 1
            ZI = SP_L(1,JPA(I))
            ZJ = SP_L(1,JPA(J))
            IF( ZI.NE.ZJ )THEN
              IT = INT(ZI*ZJ)
              IF( ETH(IT).EQ.0.D+0) CALL ELECTS(ZI,ZJ,IT,CPIX,APHI)
            ENDIF
            CTAPH(NT) = TAA(NT)+ETH(IT)+CPIX*ETHP(IT)
            CTA(NT) = TAA(NT)+ETH(IT)
            CTAPR(NT) = ETHP(IT)
            CTAPPR(NT) = ETHP2(IT)
            NT = NT+1
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HOMIX ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IMOBCF( NEQ )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
!
!     Add immobile species concentrations into component species
!     concentrations.
!
!     C(N,NSL) component species concentration (kmol/m^3 node)
!     SP_C(N,NSP) species concentration (kmol/m^3 node)
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 2 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
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
      SUB_LOG(ISUB_LOG) = '/IMOBCF'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      L1: DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE L1
        N_DB = ND(N)
        NSL = NSOLU + NEQ
!
!---    Linked aqueous air   ---
!
        IF( ISPLK(4).EQ.NSL ) THEN
          COX = 1.D+3*XLA(1,N)*RHOL(1,N)*SL(1,N)*PORD(1,N)/WTMA
          CX = 1.D+3*XLA(2,N)*RHOL(2,N)*SL(2,N)*PORD(2,N)/WTMA
          CO(N,NSL) = COX
          C(N,NSL) = COX + (CX-COX)*DT/DT_RST
          CYCLE L1
        ENDIF
!
!---    Linked aqueous CO2   ---
!
        IF( ISPLK(6).EQ.NSL ) THEN
          COX = 1.D+3*XLA(1,N)*RHOL(1,N)*SL(1,N)*PORD(1,N)/WTMA
          CX = 1.D+3*XLA(2,N)*RHOL(2,N)*SL(2,N)*PORD(2,N)/WTMA
          CO(N,NSL) = COX
          C(N,NSL) = COX + (CX-COX)*DT/DT_RST

        ENDIF
!
!---    Loop over conservation-component species  ---
!
        DO M = 1,IEQ_C(1,NEQ)
          NSP = IEQ_C(M+1,NEQ)
!
!---      Solid and exchange species ---
!
          IF( NSP.GT.NSPL .AND. NSP.LE.(NSPL+NSPS+NSPE) ) THEN
            IF( ABS(SP_C(N,NSP)).LT.CMIN ) THEN
              SP_CX = 0.D+0
            ELSE
              SP_CX = SP_C(N,NSP)
            ENDIF
            C(N,NSL) = C(N,NSL) + EQ_C(M,NEQ)*SP_CX
          ENDIF
        ENDDO
!
!---    Update old time step conservation-component species ---
!
        CO(N,NSL) = C(N,NSL)
      ENDDO L1
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IMOBCF group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IMOBKF( NEQ )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
!
!     Mobile kinetic component fractions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 1 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
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
      SUB_LOG(ISUB_LOG) = '/IMOBKF'
!
!---  Convert to global component indices  ---
!
      NEQX = NEQ + NEQC
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      L1: DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE L1
        N_DB = ND(N)
        NSL = NSOLU + NEQX
!
!---    Linked aqueous air   ---
!
        IF( ISPLK(4).EQ.NSL ) THEN
          COX = 1.D+3*XLA(1,N)*RHOL(1,N)*SL(1,N)*PORD(1,N)/WTMA
          CX = 1.D+3*XLA(2,N)*RHOL(2,N)*SL(2,N)*PORD(2,N)/WTMA
          CO(N,NSL) = COX
          C(N,NSL) = COX + (CX-COX)*DT/DT_RST
          CYCLE L1
        ENDIF
!
!---    Linked aqueous CO2   ---
!
        IF( ISPLK(6).EQ.NSL ) THEN
          COX = 1.D+3*XLA(1,N)*RHOL(1,N)*SL(1,N)*PORD(1,N)/WTMA
          CX = 1.D+3*XLA(2,N)*RHOL(2,N)*SL(2,N)*PORD(2,N)/WTMA
          CO(N,NSL) = COX
          C(N,NSL) = COX + (CX-COX)*DT/DT_RST
          CYCLE L1
        ENDIF

!
!---    Loop over kinetic-component species  ---
!
        DO M = 1,IEQ_K(1,NEQ)
          NSP = IEQ_K(M+1,NEQ)
!
!---      Solid and exchange species ---
!
          IF( NSP.GT.NSPL .AND. NSP.LE.(NSPL+NSPS+NSPE) ) THEN
            IF( ABS(SP_C(N,NSP)).LT.CMIN ) THEN
              SP_CX = 0.D+0
            ELSE
              SP_CX = SP_C(N,NSP)
            ENDIF
            C(N,NSL) = C(N,NSL) + EQ_K(M,NEQ)*SP_CX
          ENDIF
        ENDDO
!
!---    Update old time step kinetic-component species ---
!
        CO(N,NSL) = C(N,NSL)
      ENDDO L1
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IMOBKF group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_PTZRCOEF
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
!     ECKEChemX
!
!     Initialize Pitzer coefficients
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 February 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE PTZRCOEF
      USE GRID
      USE GLB_PAR
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
      SUB_LOG(ISUB_LOG) = '/INTLZ_PTZRCOEF'
      NCC_PZ = 0
      NA_PZ = 0
      NNN_PZ = 0
      DO L = 1,LANI
        JPA(L) = 0
      ENDDO
      DO L = 1,LCAT
        JPC(L) = 0
      ENDDO
      DO L = 1,LNEU
        JPN(L) = 0
      ENDDO
      DO L = 1,LANI
        DO K = 1,LCAT
          B0(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LANI
        DO K = 1,LCAT
          B1(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LANI
        DO K = 1,LCAT
          B2(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LANI
        DO K = 1,LCAT
          CMXX(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LNCF
        TCC(L) = 0.D+0
      ENDDO
      DO L = 1,LNAF
        TAA(L) = 0.D+0
      ENDDO
      DO L = 1,LANI
        DO K = 1,LNCF
          PSIC(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LCAT
        DO K = 1,LNAF
          PSIA(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LANI
        DO K = 1,LNEU
          ALAMB(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LCAT
        DO K = 1,LNEU
          CLAMB(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LNEU
        DO K = 1,LNEU
          ELAMB(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LANI
        DO K = 1,LNNF
          HOLAMB(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LANI
        DO K = 1,LCAT
          BPPR(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LANI
        DO K = 1,LCAT
          BPHI(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LANI
        DO K = 1,LCAT
          BPR(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LANI
        DO K = 1,LCAT
          BMMX(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO M = 1,8
        DO L = 1,LANI
          DO K = 1,LCAT
            ATB0(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO M = 1,8
        DO L = 1,LANI
          DO K = 1,LCAT
            ATB1(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO M = 1,8
        DO L = 1,LANI
          DO K = 1,LCAT
            ATB2(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO M = 1,8
        DO L = 1,LANI
          DO K = 1,LCAT
            ATCMX(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO M = 1,6
        DO L = 1,LNEU
          DO K = 1,LNEU
            ATNLAM(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO M = 1,6
        DO L = 1,LCAT
          DO K = 1,LNEU
            ATCLAM(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO M = 1,6
        DO L = 1,LANI
          DO K = 1,LNEU
            ATALAM(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO M = 1,6
        DO L = 1,LANI
          DO K = 1,LNNF
            ATHLAM(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,6
        DO K = 1,LNCF
          ATTC(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO M = 1,6
        DO L = 1,LANI
          DO K = 1,LNCF
            ATPC(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,6
        DO K = 1,LNAF
          ATTA(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO M = 1,6
        DO L = 1,LCAT
          DO K = 1,LNAF
            ATPA(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,LNCF
        CTCPH(L) = 0.D+0
      ENDDO
      DO L = 1,LNCF
        CTC(L) = 0.D+0
      ENDDO
      DO L = 1,LNCF
        CTCPR(L) = 0.D+0
      ENDDO
      DO L = 1,LNCF
        CTCPPR(L) = 0.D+0
      ENDDO
      DO L = 1,LNAF
        CTAPH(L) = 0.D+0
      ENDDO
      DO L = 1,LNAF
        CTA(L) = 0.D+0
      ENDDO
      DO L = 1,LNAF
        CTAPR(L) = 0.D+0
      ENDDO
      DO L = 1,LNAF
        CTAPPR(L) = 0.D+0
      ENDDO
      DO L = 1,LMCG*LMCG
        ETH(L) = 0.D+0
      ENDDO
      DO L = 1,LMCG*LMCG
        ETHP(L) = 0.D+0
      ENDDO
      DO L = 1,LMCG*LMCG
        ETHP2(L) = 0.D+0
      ENDDO
      DO L = 1,LSPL
        IDD_PZ(L) = 0
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_PTZRCOEF group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KECHEM( ACTVX,AJM,BJM,COX,SP_CX,DSP_CX,N,NEQ,INDX,IER )
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
!     ECKEChemX
!
!     Kinetic Equation CHEMistry
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
      USE PROP
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
      REAL*8 ACTVX(LSPL+1,LSPL)
      REAL*8 BJM(LSPR),AJM(LSPR,LSPR)
      REAL*8 SP_CX(LSPR),DSP_CX(LSPR)
      REAL*8 COX(LEQC+LEQK)
      REAL*8 EQKX(LREK),RRBX(LREK),RRCX(LREK)
      LOGICAL FCHK
      CHARACTER*64 GETSPNM
      EXTERNAL GETSPNM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/KECHEM'
!
!---  Skip for initial pH  ---
!
      IF( ISPLK(1).GT.1000 ) THEN
        IF( NEQ.EQ.IEQ_S(MOD(ISPLK(1),1000)) .AND.
     &    NSTEP.EQ.0 ) THEN
          AJM(NEQ,NEQ) = 1.D+0
          BJM(NEQ) = 0.D+0
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
      ENDIF
!
!---  Fixed species concentration or activity  ---
!
      DO NSLKX = 1,NSPLK
        IF( ISPLK(14+NSLKX).LT.0 ) THEN
          NSPX = ABS(ISPLK(14+NSLKX))
          IF( NSPX.GT.1000 ) NSPX = NSPX - 1000
          IF( NEQ.EQ.IEQ_S(NSPX) ) THEN
            AJM(NEQ,NEQ) = 1.D+0
            BJM(NEQ) = 0.D+0
            ISUB_LOG = ISUB_LOG-1
            RETURN
          ENDIF
        ENDIF
      ENDDO
!
!---  Volumetric concentration to molality, mol/m^3 -> mol/kg aqu  ---
!
      VTOMX = 1.D+0/(SL(2,N)*PORD(2,N)*RHOL(2,N)*XLW(2,N))
!
!---  Volumetric concentration to aqueous concentration,
!     mol/m^3 -> mol/m^3 aqu  ---
!
      VTOLX = 1.D+0/(SL(2,N)*PORD(2,N))
!
!---  Volumetric concentration to sorbed concentration,
!     mol/m^3 -> mol/kg sol  ---
!
      VTOSX = 1.D+0/((1.D+0-PORT(2,N))*RHOS(N))
!
!---  Total number of species  ---
!
      NEQX = NEQ - NEQE - NEQC
      NSPR = NSPG + NSPL + NSPN + NSPS + NSPE
!
!---  Number of species and reactions in kinetic equation  ---
!
      NS = IEQ_K(1,NEQX)
      NR = IEQ_K(NS+2,NEQX)
      IF( NSTEP.EQ.0 ) THEN
        NR = 0
        NSPRX = 0
        NSPPX = 0
        NSPKX = 0
      ENDIF
!
!---  Find aqueous silica
!
      DO M = NSPG+1,NSPG+NSPL
        IF (GETSPNM(M) == 'sio2(aq)') NSPSI = M
      END DO
!
!---  Loop over kinetic reactions in kinetic equation to
!     determine rate constants  ---
!
      DO M = 1,NR
!
!---    Reaction index, number of reactants, number of products  ---
!
        IRCX = IEQ_K(NS+2+M,NEQX)
        NSPRX = IRC_K(1,IRCX)
        NSPPX = IRC_K(2,IRCX)
        NSPKX = NSPRX+NSPPX
!
!---    Dissolution-precipitation kinetic reaction  ---
!
        IF( (IRCKT(IRCX).GE.10 .AND. IRCKT(IRCX).LE.12) .OR.
     &    (IRCKT(IRCX).EQ.14) .OR.
     &    (IRCKT(IRCX).GE.5 .AND. IRCKT(IRCX).LE.9) ) THEN
!
!---      Equilibrium constants as a function of temperature  ---
!
          IRCX = -IRCX
          CALL EQCN( EQKX(M),T(2,N),IRCX,N )
          IRCX = -IRCX
!
!---      Reaction rate constants as a function of temperature
!         mol/m^2 s  ---
!
          TKX = T(2,N)+TABS
          N1 = MAX(1,N*IRCKN(NSPKX+1))
          N2 = MAX(1,N*IRCKN(NSPKX+2))
          N3 = MAX(1,N*IRCKN(NSPKX+3))
          TKRX = RC_K(NSPKX+3,N3,IRCX)+TABS
          RRCX(M) = RC_K(NSPKX+1,N1,IRCX)*EXP( -RC_K(NSPKX+2,N2,IRCX)*
     &      ((1.D+0/TKX)-(1.D+0/TKRX))/(1.D-3*RCU) )
        ENDIF
      ENDDO
!
!---  Base residual  ---
!
      BJM(NEQ) = 0.D+0
      RSBX = 0.D+0
!
!---  Loop over kinetic reactions  ---
!
      DO M = 1,NR
!
!---    Reaction index, number of reactants, number of products  ---
!
        IRCX = IEQ_K(NS+2+M,NEQX)
        NSPRX = IRC_K(1,IRCX)
        NSPPX = IRC_K(2,IRCX)
        NSPKX = NSPRX+NSPPX
!
!---    TST kinetic reaction  ---
!
        IF( (IRCKT(IRCX).GE.10 .AND. IRCKT(IRCX).LE.12) .OR.
     &    (IRCKT(IRCX).EQ.14) .OR.
     &    (IRCKT(IRCX).GE.5 .AND. IRCKT(IRCX).LE.9) ) THEN
!
!---      Ion activity product mol/kg water, loop over species in
!         kinetic reaction  ---
!
          QX = 1.D+0
!
!---      Glass 1 dependent on aqueous silica only
!
          IF (IRCKT(IRCX).EQ.14) THEN
            CMX = SP_CX(NSPSI)*VTOMX
            ACX = ACTVX(1,NSPSI)
            QX = QX*(CMX*ACX)
          ELSE
!
!---        Loop over species in kinetic reaction  ---
!
            DO L = 1,NSPKX
              NSPX = IRC_K(L+2,IRCX)
!
!---          Aqueous species,
!             concentration in molality, mol/kg H2O  ---
!
              IF( NSPX.LE.NSPL ) THEN
                CMX = SP_CX(NSPX)*VTOMX
                ACX = ACTVX(1,NSPX)
!
!---            Reactants  ---
!
                N1 = MAX(1,N*IRCKN(L))
                IF( L.LE.NSPRX ) THEN
                  QX = QX*((CMX*ACX)**RC_K(L,N1,IRCX))
!
!---            Products  ---
!
                ELSE
                  QX = QX/((CMX*ACX)**RC_K(L,N1,IRCX))
                ENDIF
!
!---            CFMX is scaling factor to translate between pore-scale
!               and macro-scale simulations.  Default = 1
!
                IF (ISLC(58).EQ.1) THEN
                  QX = CFMX(N)*QX
                ENDIF
!
!---          Solid species, skip  ---
!
              ELSEIF( NSPX.LE.NSPL+NSPS ) THEN
                CYCLE
              ENDIF
            ENDDO
          ENDIF
!
!---      Initial reactive surface area, initial mineral volume 
!         fraction, current mineral volume fraction, minimum current 
!         mineral volume fraction allows re-precipitation of dissolved 
!         primary minerals NSP_M - mineral species number  ---
!
          NSPX = IRC_K(3+NSPKX,IRCX)
          NSP_M =  NSPX - NSPL
!
!---      Primary mineral  ---
!
          IF( RS_S(2,NSP_M,N).GT.EPSL ) THEN
            AOX = RS_S(1,NSP_M,N)*VOL(N)*RS_S(2,NSP_M,N)*SP_S(1,NSP_M)
            VFMOX = RS_S(2,NSP_M,N)
            IF( ISP_MN(NSPX).EQ.1 ) THEN
              VFMX = 1.D-3*(SP_CMN(N,NSP_M)+SP_CX(NSPX))
     &          *SP_S(2,NSP_M)/SP_S(1,NSP_M)
            ELSE
              VFMX = 1.D-3*SP_CX(NSPX)*SP_S(2,NSP_M)/SP_S(1,NSP_M)
            ENDIF
            IF( (IRCKT(IRCX).NE.7) .AND. (IRCKT(IRCX).NE.9) .AND.
     &          (IRCKT(IRCX).NE.12) .AND. (IRCKT(IRCX).NE.14) )
     &      VFMX = MAX( VFMX,1.D-5 )
     
!
!---      Secondary mineral, initial reactive surface area
!         for seconary minerals is set to 0.25 m^2/dm^3  ---
!
          ELSE
            IF( RS_S(1,NSP_M,N).GT.EPSL ) THEN
              VFMOX = 1.D-5
              AOX = RS_S(1,NSP_M,N)*VOL(N)*VFMOX*SP_S(1,NSP_M)
            ELSE
              AOX = 0.25D+3*VOL(N)
              VFMOX = 1.D-2
            ENDIF
            VFMX = 1.D-3*SP_CX(NSPX)*SP_S(2,NSP_M)/SP_S(1,NSP_M)
            IF( (IRCKT(IRCX).NE.7) .AND. (IRCKT(IRCX).NE.9) .AND.
     &          (IRCKT(IRCX).NE.12) .AND. (IRCKT(IRCX).NE.14) )
     &      VFMX = MAX( VFMX,1.D-5 )
          ENDIF
          VFMX = MAX( VFMX,0.D+0 )
!
!---      Reactive surface area  ---
!
          AX = AOX*(((POR0(2,N)*VFMX)/
     &      (POR(2,N)*VFMOX))**(2.D+0/3.D+0))
          IF( ISLC(56).EQ.1 ) AX = AX * SL(2,N)
          IF( ISLC(56).EQ.2 ) AX = AOX
          IF( ISLC(58).EQ.1 ) AX = 1.0D+0
!
!---      Reaction rate, mol/s  ---
!
          RRBX(M) = -AX*RRCX(M)*(1.D+0-(QX/EQKX(M)))
          IF( RRBX(M).NE.RRBX(M) ) THEN
            IER = -1
            RETURN
          ENDIF
!
!---      pH dependence  ---
!
          IF( ((IRCKT(IRCX).GE.5 .AND. IRCKT(IRCX).LE.9)
     &    .OR. (IRCKT(IRCX).EQ.14))
     &      .AND. ISPLK(1).NE.0 ) THEN
            NSP_PHX = MOD(ISPLK(1),1000)
            PHX = -LOG10(1.D-3*SP_CX(NSP_PHX)*VTOLX)
            IF( IRCKT(IRCX).GE.8 .AND. IRCKT(IRCX).LE.9 ) THEN
              RRBX(M) = RRBX(M)*MAX( 0.D+0,
     &          (7.9201D-1 - 1.3479D-1*PHX + 5.2D-3*(PHX**2)))
            ELSE
              N9 = MAX(1,N*IRCKN(NSPKX+9))
              RRBX(M) = RRBX(M)*(1.D+1**(-RC_K(NSPKX+9,N9,IRCX)*PHX))
            ENDIF
          ENDIF
!
!---      Reaction rate, mol/m^3 aqu s  ---
!
          RRBX(M) = RRBX(M)*VTOLX/VOL(N)
          SP_AREA(N,NSP_M) = AX
          NSP_MIN = IEQ_K(2,NEQX)-NSPL
          IF (M.EQ.1) THEN
            SP_RATE(N,NSP_MIN) = RRBX(M)/VTOLX*VOL(N) * EQ_K(M+1,NEQX)
          ELSE
            SP_RATE(N,NSP_MIN) = SP_RATE(N,NSP_MIN) 
     &        + RRBX(M)/VTOLX*VOL(N) * EQ_K(M+1,NEQX)
          ENDIF
!
!---      Direction limited  ---
!
          IF( IRCKT(IRCX).EQ.6 .OR. IRCKT(IRCX).EQ.8
     &      .OR. IRCKT(IRCX).EQ.11 ) THEN
            RRBX(M) = MAX( RRBX(M),0.D+0 )
          ELSEIF( IRCKT(IRCX).EQ.7 .OR. IRCKT(IRCX).EQ.9
     &      .OR. IRCKT(IRCX).EQ.12 .OR. IRCKT(IRCX).EQ.14 ) THEN
            RRBX(M) = MIN( RRBX(M),0.D+0 )
          ENDIF
!
!---    Forward-backward kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.1 ) THEN
          N1 = MAX(1,N*IRCKN(NSPKX+1))
          N2 = MAX(1,N*IRCKN(NSPKX+2))
          FRRX = RC_K(NSPKX+1,N1,IRCX)
          BRRX = RC_K(NSPKX+2,N2,IRCX)
!
!---      Loop over reactants  ---
!
          DO L = 1,NSPRX
!
!---        Global species index  ---
!
            NSPX = IRC_K(L+2,IRCX)
            CMX = SP_CX(NSPX)*VTOLX
            N1 = MAX(1,N*IRCKN(L))
            FRRX = FRRX*(CMX**RC_K(L,N1,IRCX))
          ENDDO
!
!---      Loop over products  ---
!
          DO L = 1,NSPPX
!
!---        Global species index  ---
!
            NSPX = IRC_K(L+2+NSPRX,IRCX)
            CMX = SP_CX(NSPX)*VTOLX
            N1 = MAX(1,N*IRCKN(L+NSPRX))
            BRRX = BRRX*(CMX**RC_K(L+NSPRX,N1,IRCX))
          ENDDO
          RRBX(M) = FRRX - BRRX
!
!---    Valocchi-Monod kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.2 ) THEN
!
!---      Concentration of biomass in mol/m^3 aqu  ---
!
          NSPX = IRC_K(5,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
!
!---      Partial rate of donor degradation  ---
!
          N6 = MAX(1,N*IRCKN(6))
          RRBX(M) = RC_K(6,N6,IRCX)*CMX
!
!---      Concentration of donor in mol/kg water  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Partial rate of donor degradation  ---
!
          N4 = MAX(1,N*IRCKN(4))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---      Concentration of acceptor in mol/kg water  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Rate of donor degradation  ---
!
          N5 = MAX(1,N*IRCKN(5))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---    Schroth-Monod kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.17 ) THEN
!
!---      Concentration of biomass in mol/m^3 aqu  ---
!
          NSPX = IRC_K(5,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
!
!---      Partial rate of donor degradation  ---
!
          N6 = MAX(1,N*IRCKN(6))
          RRBX(M) = RC_K(6,N6,IRCX)*CMX
!
!---      Concentration of donor in mol/kg water  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Partial rate of donor degradation  ---
!
          N4 = MAX(1,N*IRCKN(4))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---      Concentration of acceptor in mol/kg water  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Rate of donor degradation  ---
!
          N5 = MAX(1,N*IRCKN(5))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---      Summation of limiter concentrations  ---
!
          SCMX = 0.D+0
          DO L = 1,IRC_K(5,IRCX)
            NSPX = IRC_K(5+L,IRCX)
            SCMX = SCMX + SP_CX(NSPX)*VTOMX
          ENDDO
!
!---      Rate of total limiter control  ---
!
          N6 = MAX(1,N*IRCKN(6))
          RRBX(M) = RRBX(M)*(1.D+0-(SCMX/(RC_K(6,N6,IRCX)+SCMX)))
!
!---    Valocchi-Sorption kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.3 ) THEN
!
!---      Concentration of sorbed species in mol/kg soil  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOSX
          N4 = MAX(1,N*IRCKN(4))
          RRBX(M) = CMX/RC_K(4,N4,IRCX)
!
!---      Concentration of aqueous species in mol/m^3 water  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
!
!---      Rate of sorption in mol/m^3 aqu s  ---
!
          N3 = MAX(1,N*IRCKN(3))
          RRBX(M) = RC_K(3,N3,IRCX)*(CMX-RRBX(M))
!
!---    Langmuir-Sorption kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.13 ) THEN
!
!---      Concentration of sorbed species in mol/kg soil  ---
!
          NSPX = IRC_K(4,IRCX)
          CSX = SP_CX(NSPX)*VTOSX
!
!---      Concentration of aqueous species in mol/m^3 water  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
!
!---      Rate of sorption in mol/m^3 aqu s  ---
!
          N3 = MAX(1,N*IRCKN(3))
          N4 = MAX(1,N*IRCKN(4))
          N5 = MAX(1,N*IRCKN(5))
          RRBX(M) = RC_K(3,N3,IRCX)*CMX*(RC_K(5,N5,IRCX)-CSX)-
     &      RC_K(4,N4,IRCX)*CSX
!
!---    Valocchi-Biomass kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.4 ) THEN
!
!---      Concentration of biomass in mol/m^3 aqu  ---
!
          NSPX = IRC_K(5,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
!
!---      Partial rate of donor degradation  ---
!
          N6 = MAX(1,N*IRCKN(6))
          RRBX(M) = RC_K(6,N6,IRCX)*CMX
!
!---      Concentration of donor in mol/m^3 aqu  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Partial rate of donor degradation  ---
!
          N4 = MAX(1,N*IRCKN(4))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---      Concentration of acceptor in mol/m^3 aqu  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Rate of donor degradation  ---
!
          N5 = MAX(1,N*IRCKN(5))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---      Concentration of biomass in mol/m^3 aqu  ---
!
          NSPX = IRC_K(5,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
!
!---      Rate of biomass production in mol/m^3 aqu s  ---
!
          N7 = MAX(1,N*IRCKN(7))
          N8 = MAX(1,N*IRCKN(8))
          RRBX(M) = RC_K(7,N7,IRCX)*RRBX(M) - RC_K(8,N8,IRCX)*CMX
!
!---    Emulsion- or oil-sorption kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.15 ) THEN
!
!---      Neighborhood-particle dimensionless number  ---
!
          GAMMAX = (MAX( 1.D+0-PORD(2,N),0.D+0 ))**(THIRD)
          ASX = 2.D+0*(1.D+0-(GAMMAX**5))/(2.D+0 - 3.D+0*GAMMAX
     &      + 3.D+0*(GAMMAX**5) - 2.D+0*(GAMMAX**6))
!
!---      Hamaker constant (J)  ---
!
          CHX = 1.D-20
!
!---      Bolztmann constant (kg m^2/s^2 K)  ---
!
          CBX = 1.38D-23
!
!---      Average aqueous velocity (m/s)  ---
!
          ULAVX = 0.D+0
          NC = 0
          IF( ABS(UL(1,1,N)).GT.EPSL ) THEN
            NC = NC+1
            ULAVX = ULAVX + UL(1,1,N)
          ENDIF
          IF( ABS(UL(1,2,N)).GT.EPSL ) THEN
            NC = NC+1
            ULAVX = ULAVX + UL(1,2,N)
          ENDIF
          IF( ABS(VL(1,1,N)).GT.EPSL ) THEN
            NC = NC+1
            ULAVX = ULAVX + VL(1,1,N)
          ENDIF
          IF( ABS(VL(1,2,N)).GT.EPSL ) THEN
            NC = NC+1
            ULAVX = ULAVX + VL(1,2,N)
          ENDIF
          IF( ABS(WL(1,1,N)).GT.EPSL ) THEN
            NC = NC+1
            ULAVX = ULAVX + WL(1,1,N)
          ENDIF
          IF( ABS(WL(1,2,N)).GT.EPSL ) THEN
            NC = NC+1
            ULAVX = ULAVX + WL(1,2,N)
          ENDIF
          ULAVX = ULAVX/REAL(NC)
!
!---      London - van der Waals dimensionless number  ---
!
          N6 = MAX(1,N*IRCKN(6))
          DNLOX = 4.D+0*CHX/
     &      (9.D+0*GPI*VISL(2,N)*(RC_K(6,N6,IRCX)**2)*ULAVX)
!
!---      Inception dimensionless number  ---
!
          N3 = MAX(1,N*IRCKN(3))
          DNRX = RC_K(6,N6,IRCX)/RC_K(3,N3,IRCX)
!
!---      Sedimentation dimensionless number  ---
!
          N7 = MAX(1,N*IRCKN(7))
          DNGX = GRAV*(RHOL(2,N)-RC_K(7,N7,IRCX))*(RC_K(6,N6,IRCX)**2)
     &      /(1.8D+1*VISL(2,N)*ULAVX)
!
!---      Diffusion dimensionless number  ---
!
          DNPEX = 3.D+0*GPI*VISL(2,N)*ULAVX*RC_K(3,N3,IRCX)
     &      *RC_K(6,N6,IRCX)/(CBX*(T(2,N)+TABS))
!
!---      Single collector efficiency  ---
!
          ETAX = ASX*((DNLOX**(1.D+0/8.D+0))*(DNRX**(1.5D+1/8.D+0)) +
     &      3.38D-3*(DNGX**1.2D+0)*(DNRX**(-4.D-1)) +
     &      4.D+0*(DNPEX**(-2.D+0/3.D+0)))
!
!---      Concentration of immobile oil (kg oil/kg soil)  ---
!
          NSPX = IRC_K(3,IRCX)
          WTMX = SP_S(2,(NSPX-NSPL))
          CIMX = 1.D-3*SP_CX(NSPX)*VTOSX*WTMX
!
!---      Concentration of mobile oil (kg oil/m^3 aqu)  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = 1.D-3*SP_CX(NSPX)*VTOLX*WTMX
!
!---      Rate of immobile-oil production in kg oil/kg soil s  ---
!
          N4 = MAX(1,N*IRCKN(4))
          N5 = MAX(1,N*IRCKN(5))
          RRBX(M) = 3.D+0*ULAVX*ETAX*MAX(1.D+0-PORD(2,N),0.D+0)*
     &      RC_K(4,N4,IRCX)*MAX(RC_K(5,N5,IRCX)-CIMX,0.D+0)*CMX/
     &      (2.D+0*RC_K(3,N3,IRCX)*RC_K(5,N5,IRCX))
!
!---      Rate of immobile-oil production in mol/m^3 aqu s  ---
!
          RRBX(M) = 1.D+3*RRBX(M)*VTOLX/(VTOSX*WTMX)
!
!---    Multirate  ---
!
        ELSEIF( IRCKT(IRCX).EQ.20 ) THEN
!
!---      Neutral reaction rate, mol/m^2 s  ---
!
          TKX = T(2,N)+TABS
          N1 = MAX(1,N*IRCKN(1))
          N2 = MAX(1,N*IRCKN(2))
          N3 = MAX(1,N*IRCKN(3))
          TKRX = RC_K(3,N3,IRCX)+TABS
          RRCX(M) = RC_K(1,N1,IRCX)*EXP( -RC_K(2,N2,IRCX)*
     &      ((1.D+0/TKX)-(1.D+0/TKRX))/(1.D-3*RCU) )
!
!---      Loop over mechanisms  ---
!
          DO NKRMX = 1,IRC_K(2,IRCX)
            IX = 3+((NKRMX-1)*6)
!
!---        Ion activity product mol/kg water, loop over species  ---
!
            QX = 1.D+0
            DO L = 1,IRC_K(IX,IRCX)
              IX = 3+((NKRMX-1)*6)+L
              NX = MAX(1,N*IRCKN(IX))
              NSPX = IRC_K(IX,IRCX)
!
!---          Aqueous species,
!             concentration in molality, mol/kg H2O  ---
!
              IF( NSPX.LE.NSPL ) THEN
                CMX = SP_CX(NSPX)*VTOMX
                ACX = ACTVX(1,NSPX)
                IX = 6+((NKRMX-1)*8)+L
                QX = QX*((CMX*ACX)**RC_K(IX,NX,IRCX))
              ENDIF
            ENDDO
            IX = 4+((NKRMX-1)*8)
            NX = MAX(1,N*IRCKN(IX))
            N1 = MAX(1,N*IRCKN(IX+1))
            N2 = MAX(1,N*IRCKN(IX+2))
            TKRX = RC_K(IX+2,N2,IRCX)+TABS
            RRCX(M) = RRCX(M) + RC_K(IX,NX,IRCX)
     &        *EXP( -RC_K(IX+1,N1,IRCX)*
     &      ((1.D+0/TKX)-(1.D+0/TKRX))/(1.D-3*RCU) )*QX
          ENDDO
!
!---      Initial reactive surface area, initial mineral volume 
!         fraction, current mineral volume fraction, minimum current
!         mineral volume fraction allows re-precipitation of dissolved
!         primary minerals NSP_M - mineral species number  ---
!
          NSPX = IRC_K(1,IRCX)
          NSP_M =  NSPX - NSPL
!
!---      Primary mineral  ---
!
          IF( RS_S(2,NSP_M,N).GT.EPSL ) THEN
            AOX = RS_S(1,NSP_M,N)*VOL(N)*RS_S(2,NSP_M,N)*SP_S(1,NSP_M)
            VFMOX = RS_S(2,NSP_M,N)
            IF( ISP_MN(NSPX).EQ.1 ) THEN
              VFMX = 1.D-3*(SP_CMN(N,NSP_M)+SP_CX(NSPX))
     &          *SP_S(2,NSP_M)/SP_S(1,NSP_M)
            ELSE
              VFMX = 1.D-3*SP_CX(NSPX)*SP_S(2,NSP_M)/SP_S(1,NSP_M)
            ENDIF
            IF( (IRCKT(IRCX).NE.7) .AND. (IRCKT(IRCX).NE.9) .AND.
     &          (IRCKT(IRCX).NE.12) .AND. (IRCKT(IRCX).NE.14) )
     &      VFMX = MAX( VFMX,1.D-5 )
!
!---      Secondary mineral, initial reactive surface area
!         for seconary minerals is set to 0.25 m^2/dm^3  ---
!
          ELSE
            IF( RS_S(1,NSP_M,N).GT.EPSL ) THEN
              VFMOX = 1.D-5
              AOX = RS_S(1,NSP_M,N)*VOL(N)*VFMOX*SP_S(1,NSP_M)
            ELSE
              AOX = 0.25D+3*VOL(N)
              VFMOX = 1.D-2
            ENDIF
            VFMX = 1.D-3*SP_CX(NSPX)*SP_S(2,NSP_M)/SP_S(1,NSP_M)
            IF( (IRCKT(IRCX).NE.7) .AND. (IRCKT(IRCX).NE.9) .AND.
     &          (IRCKT(IRCX).NE.12) .AND. (IRCKT(IRCX).NE.14) )
     &      VFMX = MAX( VFMX,1.D-5 )
          ENDIF
          VFMX = MAX( VFMX,0.D+0 )
!
!---      Reactive surface area  ---
!
          AX = AOX*(((POR0(2,N)*VFMX)/
     &      (POR(2,N)*VFMOX))**(2.D+0/3.D+0))
          IF( ISLC(56).EQ.1 ) AX = AX * SL(2,N)
          IF( ISLC(56).EQ.2 ) AX = AOX
!
!---      Reaction rate, mol/s  ---
!
          RRBX(M) = -AX*RRCX(M)
!
!---      Reaction rate, mol/m^3 aqu s  ---
!
          RRBX(M) = RRBX(M)*VTOLX/VOL(N)
          SP_AREA(N,NSP_M) = AX
          NSP_MIN = IEQ_K(2,NEQX)-NSPL
          IF (M.EQ.1) THEN
            SP_RATE(N,NSP_MIN) = RRBX(M)/VTOLX*VOL(N) * EQ_K(M+1,NEQX)
          ELSE
            SP_RATE(N,NSP_MIN) = SP_RATE(N,NSP_MIN) 
     &                         + RRBX(M)/VTOLX*VOL(N) * EQ_K(M+1,NEQX)
          ENDIF
!
!---    Monod kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.22 ) THEN
          JCX = 0
          RRBX(M) = 1.D+0
!
!---      Loop over the number of reactants, less one  ---
!
          DO NSP = 1,NSPRX-1
!
!---        Concentration of reactant in mol/m^3 aqu  ---
!
            NSPX = IRC_K(NSP+3,IRCX)
            CMX = SP_CX(NSPX)*VTOLX
!
!---        Partial rate of donor degradation  ---
!
            JCX = JCX+1
            N1 = MAX(1,N*IRCKN(JCX))
            RRBX(M) = RRBX(M)*(CMX/(RC_K(JCX,N1,IRCX)+CMX))
          ENDDO
!
!---      Concentration of biomass in mol/m^3 aqu  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
!
!---      Rate of reactant degradation in mol/m^3 aqu s  ---
!
          JCX = JCX+1
          N1 = MAX(1,N*IRCKN(JCX))
          RRBX(M) = -RRBX(M)*RC_K(JCX,N1,IRCX)*CMX
!
!---    Biomass kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.24 ) THEN
          JCX = 0
          RRBX(M) = 0.D+0
!
!---      Loop over the number of reactants, less one  ---
!
          DO NSP = 1,NSPRX-1
!
!---        Concentration of reactant in mol/m^3 water  ---
!
            NSPX = IRC_K(NSP+3,IRCX)
            CMX = SP_CX(NSPX)*VTOLX
!
!---        Partial rate of reactant degradation  ---
!
            JCX = JCX+1
            N1 = MAX(1,N*IRCKN(JCX))
            RRBXX = (CMX/(RC_K(JCX,N1,IRCX)+CMX))
!
!---        Concentration of biomass in mol/m^3 aqu  ---
!
            NSPX = IRC_K(3,IRCX)
            CMX = SP_CX(NSPX)*VTOLX
!
!---        Partial rate of biomass degradation  ---
!
            JCX = JCX+1
            N1 = MAX(1,N*IRCKN(JCX))
            RRBXX = RRBXX*RC_K(JCX,N1,IRCX)*CMX
            RRBX(M) = RRBX(M) + RRBXX
          ENDDO
!
!---      Microbial specific yield coefficient  ---
!
          JCX = JCX+1
          N1 = MAX(1,N*IRCKN(JCX))
          RRBX(M) = RRBX(M)*RC_K(JCX,N1,IRCX)
!
!---      Concentration of reactant in mol/kg water  ---
!
          NSPX = IRC_K(NSP+3,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Partial rate of microbial degradation  ---
!
          JCX = JCX+1
          N1 = MAX(1,N*IRCKN(JCX))
          RRBX(M) = RRBX(M) - RC_K(JCX,N1,IRCX)*CMX
!
!---    Dual Monod kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.35 ) THEN
!
!---      Partial rate of donor degradation  ---
!
          N5 = MAX(1,N*IRCKN(5))
          RRBX(M) = RC_K(5,N5,IRCX)
!
!---      Concentration of donor in mol/kg water  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Partial rate of donor degradation  ---
!
          N3 = MAX(1,N*IRCKN(3))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(3,N3,IRCX)+CMX))
!
!---      Concentration of acceptor in mol/kg water  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Rate of donor degradation  ---
!
          N4 = MAX(1,N*IRCKN(4))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---    Single Monod kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.36 ) THEN
!
!---      Partial rate of donor degradation  ---
!
          N4 = MAX(1,N*IRCKN(4))
          RRBX(M) = RC_K(4,N4,IRCX)
!
!---      Concentration of donor in mol/kg water  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Partial rate of donor degradation  ---
!
          N3 = MAX(1,N*IRCKN(3))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(3,N3,IRCX)+CMX))
!
!---      Concentration of acceptor in mol/kg water  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Rate of donor degradation  ---
!
          RRBX(M) = RRBX(M)*CMX
!
!---    Dual Monod with Inhibitor kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.37 ) THEN
!
!---      Mass transfer coefficient  ---
!
          N6 = MAX(1,N*IRCKN(6))
          RRBX(M) = RC_K(6,N6,IRCX)
!
!---      Concentration of donor in mol/kg water  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Half-saturation constant for electron donor  ---
!
          N3 = MAX(1,N*IRCKN(3))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(3,N3,IRCX)+CMX))
!
!---      Concentration of acceptor in mol/kg water  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Half-saturation constant for electron acceptor  ---
!
          N4 = MAX(1,N*IRCKN(4))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---      Concentration of inhibitor in mol/kg water  ---
!
          NSPX = IRC_K(5,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Half-saturation constant for inhibitor  ---
!
          N5 = MAX(1,N*IRCKN(5))
          RRBX(M) = RRBX(M)*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---    Liu's multi-rate kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.41 ) THEN
          N1 = MAX(1,N*IRCKN(NSPKX+1))
          N2 = MAX(1,N*IRCKN(NSPKX+2))
          N3 = MAX(1,N*IRCKN(NSPKX+3))
          N4 = MAX(1,N*IRCKN(NSPKX+4))
          N5 = MAX(1,N*IRCKN(NSPKX+5))
          RMX = RC_K(NSPKX+1,N1,IRCX)
          SDENX = RC_K(NSPKX+2,N2,IRCX)*VTOLX/VTOSX
          PFRCX = RC_K(NSPKX+3,N3,IRCX)
          XLGK1 = RC_K(NSPKX+4,N4,IRCX)
          XLGK2 = RC_K(NSPKX+5,N5,IRCX)
          FRRX = 1.D+0
          BRRX = 1.D+0
!
!---      Loop over reactants  ---
!
          DO L = 1,NSPRX
!
!---        Global species index  ---
!
            NSPX = IRC_K(L+2,IRCX)
            CMX = SP_CX(NSPX)*VTOMX
            IF( NSPX.LE.NSPL ) THEN
              ACX = ACTVX(1,NSPX)
            ELSE
              ACX = 1.D+0
            ENDIF
            N1 = MAX(1,N*IRCKN(L))
            FRRX = FRRX*((CMX*ACX)**RC_K(L,N1,IRCX))
          ENDDO
!
!---      Loop over products  ---
!
          DO L = 1,NSPPX-1
!
!---        Global species index  ---
!
            NSPX = IRC_K(L+2+NSPRX,IRCX)
            CMX = SP_CX(NSPX)*VTOMX
            IF( NSPX.LE.NSPL ) THEN
              ACX = ACTVX(1,NSPX)
            ELSE
              ACX = 1.D+0
            ENDIF
            N1 = MAX(1,N*IRCKN(L+NSPRX))
            BRRX = BRRX*((CMX*ACX)**RC_K(L+NSPRX,N1,IRCX))
          ENDDO
          L = NSPPX
          NSPX = IRC_K(L+2+NSPRX,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
          FRRX = FRRX*(1.D+1**XLGK1)
          BRRX = BRRX*(1.D+1**XLGK2)
          RRBX(M) = RMX*(SDENX*PFRCX*FRRX/(1.D+0+FRRX+BRRX)-CMX)
!
!---    Liu's dual domain kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.42 ) THEN
          N1 = MAX(1,N*IRCKN(NSPKX+1))
          RMX = RC_K(NSPKX+1,N1,IRCX)
!
!---      Loop over reactants  ---
!
          N3 = MAX(1,N*IRCKN(NSPKX+3))
          PFRCX = RC_K(NSPKX+3,N3,IRCX)
          FRRX = 0.D+0
!
!---      Global species index  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
          FRRX = FRRX + (CMX**VTOLX)
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
          FRRX = FRRX - (CMX**VTOLX)
          RRBX(M) = RMX*FRRX
        ENDIF
!
!---    Residual contribution  ---
!
        RSBX = RSBX + RRBX(M)*EQ_K(NS+M,NEQX)
      ENDDO
!
!---  Loop over kinetic species in kinetic equation  ---
!
      CX = 0.D+0
      CTOX = 0.D+0
      DO M = 1,NS
        NSPX = IEQ_K(M+1,NEQX)
        CX = CX + SP_CX(NSPX)*EQ_K(M,NEQX)
        IF( ISP_MN(NSPX).EQ.1 ) THEN
          NSP_M =  NSPX - NSPL
          CTOX = CTOX + (SP_CO(N,NSPX)+SP_CMN(N,NSP_M))*EQ_K(M,NEQX)
        ELSE
          CTOX = CTOX + SP_CO(N,NSPX)*EQ_K(M,NEQX)
        ENDIF
        DO NSPKX = 1,NSPLK
          IF(ISPLK(14+NSPKX).LT.-1000) THEN
            NSPXX=ABS(ISPLK(14+NSPKX))-1000
            IF(NSPX.EQ.NSPXX) THEN
             IF( NSPX.LE.NSPL ) THEN
               ACX = ACTVX(1,NSPX)
             ELSE
               ACX = 1.D0
             ENDIF
             CTOX = CTOX-SP_CO(N,NSPX)*EQ_K(M,NEQX)
     &            + SP_CO(N,NSPX)/ACX*EQ_K(M,NEQX)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      CMX = CX*VTOLX
      CMOX = COX(NEQ-NEQE)*VTOLX
      CMTOX = CTOX*VTOLX
!
!---  Check for complete consumption  ---
!
      IF( RSBX.LT.0.D+0 ) THEN
        RSBX = MAX( RSBX,(-CMTOX*DTI) )
      ENDIF
      BJM(NEQ) = (CMX-CMOX)*DTI - RSBX
!
!---  Return residual vector  ---
!
      IF( INDX.EQ.1 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Incremented residuals  ---
!
      DO NSP = 1,NSPR
!
!---  Check whether specie is a kinetic equation specie,
!     which affects the residual via the kinetic equation
!     or a kinetic reaction specie,
!     which affects the residual via the kinetic equation,
!     via the kinetic equation reactions  ---
!
        FCHK = .FALSE.
        IFIND = 0
!
!---    Fixed species concentration  ---
!
        DO NSPKX = 1,NSPLK
          IF( ISPLK(14+NSPKX).LT.0 ) THEN
            NSPXX = ABS(ISPLK(14+NSPKX))
            IF( NSP.EQ.NSPXX ) THEN
              FCHK = .FALSE.
              IFIND = 1
              EXIT
            ENDIF
          ENDIF
        ENDDO
!
!---    Loop over species in kinetic equation  ---
!
        IF( IFIND.EQ.0 ) THEN
          DO M = 1,NS
            NSPX = IEQ_K(M+1,NEQX)
!
!---        Specie is a kinetic reaction specie  ---
!
            IF( NSP.EQ.NSPX ) THEN
              FCHK = .TRUE.
              IFIND = 1
              EXIT
            ENDIF
          ENDDO
        ENDIF
!
!---    Loop over kinetic reaction reactants  ---
!
        IF( IFIND.EQ.0 ) THEN
          DO L = 1,NSPRX
!
!---        Global species index  ---
!
            NSPX = IRC_K(L+2,IRCX)
!
!---        Specie is a kinetic reaction specie  ---
!
            IF( NSP.EQ.NSPX ) THEN
              FCHK = .TRUE.
              IFIND = 1
              EXIT
            ENDIF
          ENDDO
        ENDIF
!
!---    Skip for initial pH  ---
!
        IF( ISPLK(1).GT.1000 .AND. IFIND.EQ.0 ) THEN
          IF( NSP.EQ.MOD(ISPLK(1),1000) .AND.
     &      NSTEP.EQ.0 ) FCHK = .FALSE.
        ENDIF
!
!---    Loop over kinetic reaction products  ---
!
        IF( IFIND.EQ.0 ) THEN
          DO L = 1,NSPPX
!
!---        Global species index  ---
!
            NSPX = IRC_K(L+2+NSPRX,IRCX)
!
!---        Specie is a kinetic reaction specie  ---
!
            IF( NSP.EQ.NSPX ) THEN
              FCHK = .TRUE.
              IFIND = 1
              EXIT
            ENDIF
          ENDDO
        ENDIF
!
!---    Specie is a kinetic equation specie,
!       or a kinetic equation reaction specie  ---
!
        IF( FCHK ) THEN
          RSX = 0.D+0
!
!---      Loop over kinetic reactions  ---
!
          DO M = 1,NR
!
!---        Reaction index, number of reactants, number of products  ---
!
            IRCX = IEQ_K(NS+2+M,NEQX)
            NSPRX = IRC_K(1,IRCX)
            NSPPX = IRC_K(2,IRCX)
            NSPKX = NSPRX+NSPPX
!
!---        Dissolution-precipitation kinetic reaction  ---
!
            IF( (IRCKT(IRCX).GE.10 .AND. IRCKT(IRCX).LE.12) .OR.
     &        (IRCKT(IRCX).EQ.14) .OR.
     &        (IRCKT(IRCX).GE.5 .AND. IRCKT(IRCX).LE.9) ) THEN
!
!---          Ion activity product mol/kg water, loop over species in
!             kinetic reaction  ---
!
              QX = 1.D+0
!
!---          Glass 1 dependent on aqueous silica only
!
              IF (IRCKT(IRCX).EQ.14) THEN
!
!---            Incremented species  ---
!
                IF( NSPSI.EQ.NSP ) THEN
                  SP_CXX = SP_CX(NSPSI)+DSP_CX(NSPSI)
!
!---            Unincremented species  ---
!
                ELSE
                  SP_CXX = SP_CX(NSPSI)
                ENDIF
                  CMX = SP_CXX*VTOMX
                  ACX = ACTVX(1,NSPSI)
                  QX = QX*(CMX*ACX)
              ELSE
!
!---            Loop over species in kinetic reaction  ---
!
                DO L = 1,NSPKX
                  NSPX = IRC_K(L+2,IRCX)
!
!---              Incremented species  ---
!
                  IF( NSPX.EQ.NSP ) THEN
                    SP_CXX = SP_CX(NSPX)+DSP_CX(NSPX)
!
!---              Unincremented species  ---
!
                  ELSE
                    SP_CXX = SP_CX(NSPX)
                  ENDIF
!
!---              Aqueous species,
!                 concentration in molality, mol/kg H2O  ---
!
                  IF( NSPX.LE.NSPL ) THEN
                    CMX = SP_CXX*VTOMX
                    ACX = ACTVX(1,NSPX)
!
!---                Reactants  ---
!
                    N1 = MAX(1,N*IRCKN(L))
                    IF( L.LE.NSPRX ) THEN
                      QX = QX*((CMX*ACX)**RC_K(L,N1,IRCX))
!
!---                Products  ---
!
                    ELSE
                      QX = QX/((CMX*ACX)**RC_K(L,N1,IRCX))
                    ENDIF
!
!---                CFMX is scaling factor to translate between
!                   pore-scale and macro-scale simulations.  Default = 1
!
                    IF (ISLC(58).EQ.1) THEN
                      QX = CFMX(N)*QX
                    ENDIF
!
!---              Solid species, skip  ---
!
                  ELSEIF( NSPX.LE.NSPL+NSPS ) THEN
                    CYCLE
                  ENDIF
                ENDDO
              ENDIF
!
!---          Reactive surface area
!             NSP_M - mineral species number  ---
!
              NSPX = IRC_K(3+NSPKX,IRCX)
              NSP_M = NSPX - NSPL
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                SP_CXX = SP_CX(NSPX)+DSP_CX(NSPX)
!
!---          Unincremented species  ---
!
              ELSE
                SP_CXX = SP_CX(NSPX)
              ENDIF
!
!---          Primary mineral  ---
!
              IF( RS_S(2,NSP_M,N).GT.EPSL ) THEN
                AOX = RS_S(1,NSP_M,N)*VOL(N)*
     &            RS_S(2,NSP_M,N)*SP_S(1,NSP_M)
                VFMOX = RS_S(2,NSP_M,N)
                IF( ISP_MN(NSPX).EQ.1 ) THEN
                  VFMX = 1.D-3*(SP_CMN(N,NSP_M)+SP_CXX)
     &              *SP_S(2,NSP_M)/SP_S(1,NSP_M)
                ELSE
                  VFMX = 1.D-3*SP_CXX*SP_S(2,NSP_M)/SP_S(1,NSP_M)
                ENDIF
                IF( (IRCKT(IRCX).NE.7) .AND. (IRCKT(IRCX).NE.9) .AND.
     &              (IRCKT(IRCX).NE.12) .AND. (IRCKT(IRCX).NE.14) )
     &          VFMX = MAX( VFMX,1.D-5 )
!
!---          Secondary mineral, initial reactive surface area
!             for seconary minerals is set to 0.25 m^2/dm^3  ---
!
              ELSE
                IF( RS_S(1,NSP_M,N).GT.EPSL ) THEN
                  VFMOX = 1.D-5
                  AOX = RS_S(1,NSP_M,N)*VOL(N)*VFMOX*SP_S(1,NSP_M)
                ELSE
                  AOX = 0.25D+3*VOL(N)
                  VFMOX = 1.D-2
                ENDIF
                VFMX = 1.D-3*SP_CXX*SP_S(2,NSP_M)/SP_S(1,NSP_M)
                IF( (IRCKT(IRCX).NE.7) .AND. (IRCKT(IRCX).NE.9) .AND.
     &              (IRCKT(IRCX).NE.12) .AND. (IRCKT(IRCX).NE.14) )
     &          VFMX = MAX( VFMX,1.D-5 )
              ENDIF
              VFMX = MAX( VFMX,0.D+0 )
!
!---          Reactive surface area  ---
!
              AX = AOX*(((POR0(2,N)*VFMX)/
     &          (POR(2,N)*VFMOX))**(2.D+0/3.D+0))
              IF( ISLC(56).EQ.1 ) AX = AX * SL(2,N)
              IF( ISLC(56).EQ.2 ) AX = AOX
              IF( ISLC(58).EQ.1 ) AX = 1.0D+0
!
!---          Reaction rate, mol/s  ---
!
              RRX = -AX*RRCX(M)*(1.D+0-(QX/EQKX(M)))
!
!---          pH dependence  ---
!
              IF( ((IRCKT(IRCX).GE.5 .AND. IRCKT(IRCX).LE.9)
     &          .OR. (IRCKT(IRCX).EQ.14))
     &          .AND. ISPLK(1).NE.0 ) THEN
                NSP_PHX = MOD(ISPLK(1),1000)
!
!---            Incremented species  ---
!
                IF( NSP_PHX.EQ.NSP ) THEN
                  SP_CXX = SP_CX(NSP_PHX)+DSP_CX(NSP_PHX)
!
!---            Unincremented species  ---
!
                ELSE
                  SP_CXX = SP_CX(NSP_PHX)
                ENDIF
                PHX = -LOG10(1.D-3*SP_CXX*VTOLX)
                IF( IRCKT(IRCX).GE.8 .AND. IRCKT(IRCX).LE.9 ) THEN
                  RRX = RRX*MAX( 0.D+0,
     &              (7.9201D-1 - 1.3479D-1*PHX + 5.2D-3*(PHX**2)))
                ELSE
                  N9 = MAX(1,N*IRCKN(NSPKX+9))
                  RRX = RRX*(1.D+1**(-RC_K(NSPKX+9,N9,IRCX)*PHX))
                ENDIF
              ENDIF
!
!---          Reaction rate, mol/m^3 aqu s  ---
!
              RRX = RRX*VTOLX/VOL(N)
!
!---          Direction limited  ---
!
              IF( IRCKT(IRCX).EQ.6 .OR. IRCKT(IRCX).EQ.8
     &          .OR. IRCKT(IRCX).EQ.11 ) THEN
                RRX = MAX( RRX,0.D+0 )
              ELSEIF( IRCKT(IRCX).EQ.7 .OR. IRCKT(IRCX).EQ.9
     &          .OR. IRCKT(IRCX).EQ.12 .OR. IRCKT(IRCX).EQ.14 ) THEN
                RRX = MIN( RRX,0.D+0 )
              ENDIF
!
!---        Forward-backward kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.1 ) THEN
              N1 = MAX(1,N*IRCKN(NSPKX+1))
              N2 = MAX(1,N*IRCKN(NSPKX+2))
              FRRX = RC_K(NSPKX+1,N1,IRCX)
              BRRX = RC_K(NSPKX+2,N2,IRCX)
!
!---          Loop over reactants  ---
!
              DO L = 1,NSPRX
!
!---            Global species index  ---
!
                NSPX = IRC_K(L+2,IRCX)
!
!---            Incremented species  ---
!
                IF( NSPX.EQ.NSP ) THEN
                  CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---            Unincremented species  ---
!
                ELSE
                  CMX = SP_CX(NSPX)*VTOLX
                ENDIF
                N1 = MAX(1,N*IRCKN(L))
                FRRX = FRRX*(CMX**RC_K(L,N1,IRCX))
              ENDDO
!
!---          Loop over products  ---
!
              DO L = 1,NSPPX
!
!---            Global species index  ---
!
                NSPX = IRC_K(L+2+NSPRX,IRCX)
!
!---            Incremented species  ---
!
                IF( NSPX.EQ.NSP ) THEN
                  CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---            Unincremented species  ---
!
                ELSE
                  CMX = SP_CX(NSPX)*VTOLX
                ENDIF
                N1 = MAX(1,N*IRCKN(L+NSPRX))
                BRRX = BRRX*(CMX**RC_K(L+NSPRX,N1,IRCX))
              ENDDO
              RRX = FRRX - BRRX
!
!---        Valocchi-Monod kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.2 ) THEN
!
!---          Concentration of biomass in mol/m^3 aqu  ---
!
              NSPX = IRC_K(5,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOLX
              ENDIF
!
!---          Partial rate of donor degradation  ---
!
              N6 = MAX(1,N*IRCKN(6))
              RRX = RC_K(6,N6,IRCX)*CMX
!
!---          Concentration of donor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(3,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Partial rate of donor degradation  ---
!
              N4 = MAX(1,N*IRCKN(4))
              RRX = RRX*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---          Concentration of acceptor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(4,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Rate of donor degradation, mol/m^3 aqu s  ---
!
              N5 = MAX(1,N*IRCKN(5))
              RRX = RRX*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---        Schroth-Monod kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.17 ) THEN
!
!---          Concentration of biomass in mol/m^3 aqu  ---
!
              NSPX = IRC_K(5,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOLX
              ENDIF
!
!---          Partial rate of donor degradation  ---
!
              N6 = MAX(1,N*IRCKN(6))
              RRX = RC_K(6,N6,IRCX)*CMX
!
!---          Concentration of donor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(3,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Partial rate of donor degradation  ---
!
              N4 = MAX(1,N*IRCKN(4))
              RRX = RRX*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---          Concentration of acceptor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(4,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Rate of donor degradation, mol/m^3 aqu s  ---
!
              N5 = MAX(1,N*IRCKN(5))
              RRX = RRX*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---          Summation of limiter concentrations  ---
!
              SCMX = 0.D+0
              DO L = 1,IRC_K(5,IRCX)
                NSPX = IRC_K(5+L,IRCX)
!
!---            Incremented species  ---
!
                IF( NSPX.EQ.NSP ) THEN
                  CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---            Unincremented species  ---
!
                ELSE
                  CMX = SP_CX(NSPX)*VTOMX
                ENDIF
                SCMX = SCMX + CMX
              ENDDO
!
!---          Rate of total limiter control  ---
!
              N6 = MAX(1,N*IRCKN(6))
              RRX = RRX*(1.D+0-(SCMX/(RC_K(6,N6,IRCX)+SCMX)))
!
!---        Valocchi-Sorption kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.3 ) THEN
!
!---          Concentration of sorbed species in mol/gm soil  ---
!
              NSPX = IRC_K(4,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
!               CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOSX*1.D-3
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOSX
!
!---          Unincremented species  ---
!
              ELSE
!               CMX = SP_CX(NSPX)*VTOSX*1.D-3
                CMX = SP_CX(NSPX)*VTOSX
              ENDIF
              N4 = MAX(1,N*IRCKN(4))
              RRX = CMX/RC_K(4,N4,IRCX)
!
!---          Concentration of aqueous species in mol/m^3 aqu  ---
!
              NSPX = IRC_K(3,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
!               CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---          Unincremented species  ---
!
              ELSE
!               CMX = SP_CX(NSPX)*VTOMX
                CMX = SP_CX(NSPX)*VTOLX
              ENDIF
!
!---          Rate of sorption in mol/m^3 aqu s  ---
!
!             RRX = RC_K(3,IRCX)*(CMX-RRX)*1.D+3
              N3 = MAX(1,N*IRCKN(3))
              RRX = RC_K(3,N3,IRCX)*(CMX-RRX)
!
!---        Langmuir-Sorption kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.13 ) THEN
!
!---          Concentration of sorbed species in mol/gm soil  ---
!
              NSPX = IRC_K(4,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CSX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOSX
!
!---          Unincremented species  ---
!
              ELSE
                CSX = SP_CX(NSPX)*VTOSX
              ENDIF
!
!---          Concentration of aqueous species in mol/m^3 aqu  ---
!
              NSPX = IRC_K(3,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOLX
              ENDIF
!
!---          Rate of sorption in mol/m^3 aqu s  ---
!
              N3 = MAX(1,N*IRCKN(3))
              N4 = MAX(1,N*IRCKN(4))
              N5 = MAX(1,N*IRCKN(5))
              RRX = RC_K(3,N3,IRCX)*CMX*(RC_K(5,N5,IRCX)-CSX)
     &          -RC_K(4,N4,IRCX)*CSX
!
!---        Valocchi-Biomass kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.4 ) THEN
!
!---          Concentration of biomass in mol/m^3 aqu  ---
!
              NSPX = IRC_K(5,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOLX
              ENDIF
!
!---          Partial rate of donor degradation  ---
!
              N6 = MAX(1,N*IRCKN(6))
              RRX = RC_K(6,N6,IRCX)*CMX
!
!---          Concentration of donor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(3,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Partial rate of donor degradation  ---
!
              N4 = MAX(1,N*IRCKN(4))
              RRX = RRX*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---          Concentration of acceptor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(4,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Rate of donor degradation, mol/m^3 aqu s  ---
!
              N5 = MAX(1,N*IRCKN(5))
              RRX = RRX*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---          Concentration of biomass in mol/m^3 aqu  ---
!
              NSPX = IRC_K(5,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOLX
              ENDIF
!
!---          Rate of biomass production in mol/m^3 aqu s  ---
!
              N7 = MAX(1,N*IRCKN(7))
              N8 = MAX(1,N*IRCKN(8))
              RRX = RC_K(7,N7,IRCX)*RRX - RC_K(8,N8,IRCX)*CMX
!
!---        Emulsion- or oil-sorption kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.15 ) THEN
!
!---          Concentration of immobile oil (kg oil/kg soil)  ---
!
              NSPX = IRC_K(3,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CIMX = 1.D-3*(SP_CX(NSPX)+DSP_CX(NSPX))*VTOSX*WTMX
!
!---          Unincremented species  ---
!
              ELSE
                CIMX = 1.D-3*SP_CX(NSPX)*VTOSX*WTMX
              ENDIF
!
!---          Concentration of mobile oil (kg oil/m^3 aqu)  ---
!
              NSPX = IRC_K(4,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = 1.D-3*(SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX*WTMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = 1.D-3*SP_CX(NSPX)*VTOLX*WTMX
              ENDIF
!
!---          Rate of immobile-oil production in kg oil/kg soil s  ---
!
              N3 = MAX(1,N*IRCKN(3))
              N4 = MAX(1,N*IRCKN(4))
              N5 = MAX(1,N*IRCKN(5))
              RRX = 3.D+0*ULAVX*ETAX*MAX(1.D+0-PORD(2,N),0.D+0)*
     &          RC_K(4,N4,IRCX)*MAX(RC_K(5,N5,IRCX)-CIMX,0.D+0)*CMX/
     &          (2.D+0*RC_K(3,N3,IRCX)*RC_K(5,N5,IRCX))
!
!---          Rate of immobile-oil production in mol/m^3 aqu s  ---
!
              RRX = 1.D+3*RRX*VTOLX/(VTOSX*WTMX)
!
!---        Monod kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.22 ) THEN
              JCX = 0
              RRX = 1.D+0
!
!---          Loop over the number of reactants, less one  ---
!
              DO L = 1,NSPRX-1
!
!---            Concentration of reactant in mol/m^3 water  ---
!
                NSPX = IRC_K(L+3,IRCX)
!
!---            Incremented species  ---
!
                IF( NSPX.EQ.NSP ) THEN
                  CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---            Unincremented species  ---
!
                ELSE
                  CMX = SP_CX(NSPX)*VTOLX
                ENDIF
!
!---            Partial rate of reactant degradation  ---
!
                JCX = JCX+1
!               RRX = RRX*(CMX/(RC_K(JCX,IRCX)+CMX))
                N1 = MAX(1,N*IRCKN(JCX))
                RRX = RRX*(CMX**RC_K(JCX,N1,IRCX))*VTOLX/VTOSX*1.D-3
              ENDDO
!
!---          Concentration of biomass in mol/m^3 aqu  ---
!
              NSPX = IRC_K(3,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOLX
              ENDIF
!
!---          Partial rate of biomass degradation  ---
!
              JCX = JCX+1
              N1 = MAX(1,N*IRCKN(JCX))
              RRX = -RRX*RC_K(JCX,N1,IRCX)*CMX
!
!---        Biomass kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.24 ) THEN
              JCX = 0
              RRX = 0.D+0
!
!---          Loop over the number of reactants, less one  ---
!
              DO L = 1,NSPRX-1
!
!---            Concentration of reactant in mol/m^3 water  ---
!
                NSPX = IRC_K(L+3,IRCX)
!
!---            Incremented species  ---
!
                IF( NSPX.EQ.NSP ) THEN
                  CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---            Unincremented species  ---
!
                ELSE
                  CMX = SP_CX(NSPX)*VTOLX
                ENDIF
!
!---            Partial rate of reactant degradation  ---
!
                JCX = JCX+1
                N1 = MAX(1,N*IRCKN(JCX))
                RRBXX = (CMX/(RC_K(JCX,N1,IRCX)+CMX))
!
!---            Concentration of biomass in mol/m^3 aqu  ---
!
                NSPX = IRC_K(3,IRCX)
                CMX = SP_CX(NSPX)*VTOLX
!
!---            Partial rate of biomass degradation  ---
!
                JCX = JCX+1
                N1 = MAX(1,N*IRCKN(JCX))
                RRBXX = RRBXX*RC_K(JCX,N1,IRCX)*CMX
                RRX = RRX + RRBXX
              ENDDO
!
!---          Microbial specific yield coefficient  ---
!
              JCX = JCX+1
              N1 = MAX(1,N*IRCKN(JCX))
              RRX = RRX*RC_K(JCX,N1,IRCX)
!
!---          Concentration of biomass in mol/m^3 aqu  ---
!
              NSPX = IRC_K(3,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOLX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOLX
              ENDIF
!
!---          Partial rate of microbial degradation  ---
!
              JCX = JCX+1
              N1 = MAX(1,N*IRCKN(JCX))
              RRX = RRX - RC_K(JCX,N1,IRCX)*CMX
!
!---        Dual-Monod kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.35 ) THEN
!
!---          Partial rate of donor degradation  ---
!
              N5 = MAX(1,N*IRCKN(5))
              RRX = RC_K(5,N5,IRCX)
!
!---          Concentration of donor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(3,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Partial rate of donor degradation  ---
!
              N3 = MAX(1,N*IRCKN(3))
              RRX = RRX*(CMX/(RC_K(3,N3,IRCX)+CMX))
!
!---          Concentration of acceptor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(4,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Rate of donor degradation, mol/m^3 aqu s  ---
!
              N4 = MAX(1,N*IRCKN(4))
              RRX = RRX*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---        Single-Monod kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.36 ) THEN
!
!---          Partial rate of donor degradation  ---
!
              N4 = MAX(1,N*IRCKN(4))
              RRX = RC_K(4,N4,IRCX)
!
!---          Concentration of donor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(3,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Partial rate of donor degradation  ---
!
              N3 = MAX(1,N*IRCKN(3))
              RRX = RRX*(CMX/(RC_K(3,N3,IRCX)+CMX))
!
!---          Concentration of acceptor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(4,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Rate of donor degradation, mol/m^3 aqu s  ---
!
              RRX = RRX*CMX
!
!---        Dual-Monod_Inhibitor kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.37 ) THEN
!
!---          Partial rate of donor degradation  ---
!
              N6 = MAX(1,N*IRCKN(6))
              RRX = RC_K(6,N6,IRCX)
!
!---          Concentration of donor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(3,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Partial rate of donor degradation  ---
!
              N3 = MAX(1,N*IRCKN(3))
              RRX = RRX*(CMX/(RC_K(3,N3,IRCX)+CMX))
!
!---          Concentration of acceptor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(4,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Rate of inhibitor degradation, mol/m^3 aqu s  ---
!
              N5 = MAX(1,N*IRCKN(5))
              RRX = RRX*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---          Concentration of inhibitor in mol/m^3 aqu  ---
!
              NSPX = IRC_K(5,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOMX
              ENDIF
!
!---          Rate of inhibitor degradation, mol/m^3 aqu s  ---
!
              N5 = MAX(1,N*IRCKN(5))
              RRX = RRX*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---        Liu's multi-rate kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.41 ) THEN
              N1 = MAX(1,N*IRCKN(NSPKX+1))
              N2 = MAX(1,N*IRCKN(NSPKX+2))
              N3 = MAX(1,N*IRCKN(NSPKX+3))
              N4 = MAX(1,N*IRCKN(NSPKX+4))
              N5 = MAX(1,N*IRCKN(NSPKX+5))
              RMX = RC_K(NSPKX+1,N1,IRCX)
              SDENX = RC_K(NSPKX+2,N2,IRCX)*VTOLX/VTOSX
              PFRCX = RC_K(NSPKX+3,N3,IRCX)
              XLGK1 = RC_K(NSPKX+4,N4,IRCX)
              XLGK2 = RC_K(NSPKX+5,N5,IRCX)
              FRRX = 1.D+0
              BRRX = 1.D+0
!
!---          Loop over reactants  ---
!
              DO L = 1,NSPRX
!
!---            Global species index  ---
!
                NSPX = IRC_K(L+2,IRCX)
                IF( NSPX.EQ.NSP ) THEN
                  CMX = (SP_CX(NSPX)+DSP_CX(NSPX))*VTOMX
                  IF(NSPX.LE.NSPL) THEN
                    ACX = ACTVX(NSP+1,NSPX)
                  ELSE
                    ACX = 1.D+0
                  ENDIF
!
!---            Unincremented species  ---
!
                ELSE
                  CMX = SP_CX(NSPX)*VTOMX
                  IF( NSPX.LE.NSPL ) THEN
                    ACX = ACTVX(1,NSPX)
                  ELSE
                    ACX = 1.D+0
                  ENDIF
                ENDIF
                N1 = MAX(1,N*IRCKN(L))
                FRRX = FRRX*((CMX*ACX)**RC_K(L,N1,IRCX))
              ENDDO
!
!---          Loop over products  ---
!
              DO L = 1,NSPPX-1
!
!---            Global species index  ---
!
                NSPX = IRC_K(L+2+NSPRX,IRCX)
!
!---            Incremented species  ---
!
                IF( NSPX.EQ.NSP ) THEN
                  CMX = (SP_CX(NSPX) + DSP_CX(NSPX))*VTOMX
                  IF( NSPX.LE.NSPL ) THEN
                    ACX = ACTVX(NSP+1,NSPX)
                  ELSE
                    ACX = 1.D+0
                  ENDIF
!
!---            Unincremented species  ---
!
                ELSE
                  CMX = SP_CX(NSPX)*VTOMX
                  IF( NSPX.LE.NSPL ) THEN
                    ACX = ACTVX(1,NSPX)
                  ELSE
                    ACX = 1.D+0
                  ENDIF
               ENDIF
                N1 = MAX(1,N*IRCKN(L+NSPRX))
                BRRX = BRRX*((CMX*ACX)**RC_K(L+NSPRX,N1,IRCX))
              ENDDO
              L = NSPPX
              NSPX = IRC_K(L+2+NSPRX,IRCX)
!
!---          Incremented species  ---
!
              IF( NSPX.EQ.NSP ) THEN
                CMX = (SP_CX(NSPX) + DSP_CX(NSPX))*VTOLX
!
!---          Unincremented species  ---
!
              ELSE
                CMX = SP_CX(NSPX)*VTOLX
              ENDIF
              FRRX = FRRX*(1.D+1**XLGK1)
              BRRX = BRRX*(1.D+1**XLGK2)
              RRX = RMX*(SDENX*PFRCX*FRRX/(1.D+0+FRRX+BRRX)-CMX)
!
!---        Liu's dual domain kinetic reaction  ---
!
            ELSEIF( IRCKT(IRCX).EQ.42 ) THEN
              N1 = MAX(1,N*IRCKN(NSPKX+1))
              N3 = MAX(1,N*IRCKN(NSPKX+3))
              RMX = RC_K(NSPKX+1,N1,IRCX)
              PFRCX = RC_K(NSPKX+3,N3,IRCX)
              FRRX = 0.D+0
!
!---        Global species index  ---
!
              NSPX = IRC_K(3,IRCX)
              IF( NSPX.EQ.NSP ) THEN
                SP_CXX = SP_CX(NSPX) + DSP_CX(NSPX)
              ELSE
                SP_CXX = SP_CX(NSPX)
              ENDIF               
              CMX = SP_CXX*VTOLX
              FRRX = FRRX + (CMX**VTOLX)
              NSPX = IRC_K(4,IRCX)
              IF( NSPX.EQ.NSP ) THEN
                SP_CXX = SP_CX(NSPX) + DSP_CX(NSPX)
              ELSE
                SP_CXX = SP_CX(NSPX)
              ENDIF               
              CMX = SP_CXX*VTOLX
              FRRX = FRRX - (CMX**VTOLX)
              RRX = RMX*FRRX
            ENDIF
!
!---        Liu's dual domain kinetic reaction  ---
!
            IF( IRCKT(IRCX).EQ.42 ) THEN
              EQ_KX = EQ_K(NS+M,NEQX)
              IF( IMMB(NEQC+NEQX).EQ.1 ) THEN
                N2 = MAX(1,N*IRCKN(NSPKX+2))
                EQ_KX = EQ_K(NS+M,NEQX)*RC_K(NSPKX+2,N2,IRCX)
              ENDIF
!
!---        Residual contribution  ---
!
            ELSE
              RSX = RSX + RRX*EQ_K(NS+M,NEQX)
            ENDIF
          ENDDO
!
!---      Loop over kinetic species in kinetic equation  ---
!
          L1: DO M = 1,NS
            NSPX = IEQ_K(M+1,NEQX)
!
!---        Incremented species  ---
!
            IF( NSPX.EQ.NSP ) THEN
              AJM(NEQ,IEQ_S(NSP)) = EQ_K(M,NEQX)*VTOLX*DTI
!
!---          Fixed species activity  ---
!
              DO NSPKX = 1,NSPLK
                IF( ISPLK(14+NSPKX).LT.-1000 ) THEN
                  NSPXX = ABS(ISPLK(14+NSPKX))-1000
                  IF( NSP.EQ.NSPXX ) THEN
                    IF( NSP.LE.NSPL ) THEN
                      ACX = ACTVX(1,NSP)
                    ELSE
                      ACX = 1.D+0
                    ENDIF
                    AJM(NEQ,IEQ_S(NSP))=EQ_K(M,NEQX)/ACX*VTOLX*DTI
                  ENDIF
                ENDIF
              ENDDO
              EXIT L1
            ENDIF
          ENDDO L1
!
!---      Check for complete consumption  ---
!
          IF( RSX.LT.0.D+0 ) THEN
            RSX = MAX( RSX,(-CMTOX*DTI) )
          ENDIF
          AJM(NEQ,IEQ_S(NSP)) = AJM(NEQ,IEQ_S(NSP)) -
     &      (RSX-RSBX)/DSP_CX(NSP)
        ENDIF
      ENDDO
!
!---  Return residual vector and Jacobian matrix  ---
!
      IF( ABS(INDX).EQ.1 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      BJM(NEQ) = -BJM(NEQ)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KECHEM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KECHEM_R( ACTVX,AJM,BJM,COX,SP_CX,DSP_CX,N,NEQ,INDX )
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
!     Kinetic Equation CHEMistry
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 January 2005.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
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
      REAL*8 ACTVX(LSPL+1,LSPL)
      REAL*8 BJM(LSPR),AJM(LSPR,LSPR)
      REAL*8 SP_CX(LSPR),DSP_CX(LSPR)
      REAL*8 COX(LEQC+LEQK)
      REAL*8 EQKX(LREK),RRCX(LREK),RSBX(LREK)
      SAVE NR
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/KECHEM_R'
!
!---  Kinetic equation location in Jacobian matrix  ---
!
      NROW = NEQ-NEQE
!
!---  Skip for initial pH  ---
!
      IF( ISPLK(1).GT.1000 ) THEN
        IF( NEQ.EQ.IEQ_S(MOD(ISPLK(1),1000)) .AND.
     &    NSTEP.EQ.0 ) THEN
          AJM(NROW,NROW) = 1.D+0
          BJM(NROW) = 0.D+0
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
      ENDIF
!
!--- fixed species concentration or activity
!
      DO NSLKX = 1,NSPLK
        IF( ISPLK(14+NSLKX).LT.0 ) THEN
          NSPX = ABS(ISPLK(14+NSLKX))
          IF( NSPX.GT.1000 ) NSPX = NSPX - 1000
          IF( NEQ.EQ.IEQ_S(NSPX) ) THEN
            AJM(NROW,NROW) = 1.D+0
            BJM(NROW) = 0.D+0
            ISUB_LOG = ISUB_LOG-1
            RETURN
          ENDIF
        ENDIF
      ENDDO
!
!---  Volumetric concentration to molality, mol/m^3 -> mol/kg aqu  ---
!
      VTOMX = 1.D+0/(SL(2,N)*PORD(2,N)*RHOL(2,N)*XLW(2,N))
!
!---  Volumetric concentration to aqueous concentration,
!     mol/m^3 -> mol/m^3 aqu  ---
!
      VTOLX = 1.D+0/(SL(2,N)*PORD(2,N))
!
!---  Volumetric concentration to sorbed concentration,
!     mol/m^3 -> mol/kg sol  ---
!
      VTOSX = 1.D+0/((1.D+0-PORT(2,N))*RHOS(N))
!
!---  Total number of species  ---
!
      NEQX = NEQ - NEQE - NEQC
      NSPR = NSPG + NSPL + NSPN + NSPS
!
!---  Number of species and reactions in kinetic equation  ---
!
      NS = IEQ_K(1,NEQX)
      NR = IEQ_K(NS+2,NEQX)
      IF( NSTEP.EQ.0 ) THEN
        NR = 0
        NSPRX = 0
        NSPPX = 0
        NSPKX = 0
      ENDIF
!
!---  Loop over kinetic reactions in kinetic equation to
!     determine rate constants  ---
!
      DO M = 1,NR
!
!---    Reaction index, number of reactants, number of products  ---
!
        IRCX = IEQ_K(NS+2+M,NEQX)
        NSPRX = IRC_K(1,IRCX)
        NSPPX = IRC_K(2,IRCX)
        NSPKX = NSPRX+NSPPX
!
!---    Dissolution-precipitation kinetic reaction  ---
!
        IF( (IRCKT(IRCX).GE.10 .AND. IRCKT(IRCX).LE.12) .OR.
     &    (IRCKT(IRCX).GE.5 .AND. IRCKT(IRCX).LE.9) ) THEN
!
!---      Equilibrium constants as a function of temperature  ---
!
          IRCX = -IRCX
          CALL EQCN( EQKX(M),T(2,N),IRCX,N )
          IRCX = -IRCX
!
!---      Reaction rate constants as a function of temperature
!         mol/m^2 s  ---
!
          TKX = T(2,N)+TABS
          N1 = MAX(1,N*IRCKN(NSPKX+1))
          N2 = MAX(1,N*IRCKN(NSPKX+2))
          N3 = MAX(1,N*IRCKN(NSPKX+3))
          TKRX = RC_K(NSPKX+3,N3,IRCX)+TABS
          RRCX(M) = RC_K(NSPKX+1,N1,IRCX)*EXP( -RC_K(NSPKX+2,N2,IRCX)*
     &      ((1.D+0/TKX)-(1.D+0/TKRX))/(1.D-3*RCU) )
        ENDIF
      ENDDO
!
!---  Base residual  ---
!
      BJM(NROW) = 0.D+0
      RSBXX = 0.D0
!
!---  Loop over kinetic reactions  ---
!
      DO M = 1,NR
        RSBX(M) = 0.D+0
!
!---    Reaction index, number of reactants, number of products  ---
!
        IRCX = IEQ_K(NS+2+M,NEQX)
        NSPRX = IRC_K(1,IRCX)
        NSPPX = IRC_K(2,IRCX)
        NSPKX = NSPRX+NSPPX
        CALL RATER( NEQX,M,NS,N,IRCX,NSPKX,NSPRX,NSPPX,SP_CX,ACTVX,
     &      VTOMX,VTOLX,VTOSX,EQKX(M),RRCX(M),RSBX(M) )
        RSBXX = RSBXX + RSBX(M)
      ENDDO
!
!---  Loop over kinetic species in kinetic equation  ---
!
      CX = 0.D+0
      CTOX = 0.D+0
      DO M = 1,NS
        NSPX = IEQ_K(M+1,NEQX)
        CX = CX + SP_CX(NSPX)*EQ_K(M,NEQX)
        IF( ISP_MN(NSPX).EQ.1 ) THEN
          NSP_M =  NSPX - NSPL
          CTOX = CTOX + (SP_CO(N,NSPX)+SP_CMN(N,NSP_M))*EQ_K(M,NEQX)
        ELSE
          CTOX = CTOX + SP_CO(N,NSPX)*EQ_K(M,NEQX)
        ENDIF
        DO NSPKX = 1,NSPLK
          IF(ISPLK(14+NSPKX).LT.-1000) THEN
            NSPXX=ABS(ISPLK(14+NSPKX))-1000
            IF(NSPX.EQ.NSPXX) THEN
             IF( NSPX.LE.NSPL ) THEN
               ACX = ACTVX(1,NSPX)
             ELSE
               ACX = 1.D0
             ENDIF             
             CTOX = CTOX-SP_CO(N,NSPX)*EQ_K(M,NEQX)
     &            + SP_CO(N,NSPX)/ACX*EQ_K(M,NEQX)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      CMX = CX*VTOLX
      CMOX = COX(NEQ-NEQE)*VTOLX
      CMTOX = CTOX*VTOLX
!
!---  Check for complete consumption  ---
!
      BJM(NROW) = (CMX-CMOX)*DTI - RSBXX
!
!---  Return residual vector  ---
!
      IF( INDX.EQ.1 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Loop over kinetic reactions  ---
!
      DO M = 1,NR
!
!---    Reaction index, number of reactants, number of products  ---
!
        IRCX = IEQ_K(NS+2+M,NEQX)
        NSPRX = IRC_K(1,IRCX)
        NSPPX = IRC_K(2,IRCX)
        NSPKX = NSPRX+NSPPX
        DO L = 1,NSPKX
          RSX = 0.D0
          NSPX = IRC_K(L+2,IRCX)
          CMMX = SP_CX(NSPX)
          IF( NSPX.LE.NSPL ) ACX = ACTVX(1,NSPX)
!          DSP_CX(NSPX) = MAX(1.D-10,DSP_CX(NSPX))
          SP_CX(NSPX) = CMMX+DSP_CX(NSPX)
          IF( NSPX.LE.NSPL ) ACTVX(1,NSPX) = ACTVX(NSPX+1,NSPX)
          CALL RATER( NEQX,M,NS,N,IRCX,NSPKX,NSPRX,NSPPX,SP_CX,ACTVX,
     &       VTOMX,VTOLX,VTOSX,EQKX(M),RRCX(M),RSX )
          SP_CX(NSPX) = CMMX
          IF( NSPX.LE.NSPL ) ACTVX(1,NSPX) = ACX
          NEQXX=IEQ_S(NSPX)
          IF( NEQXX.GT.NEQE) THEN
            NCOL = NEQXX-NEQE
            AJM(NROW,NCOL) = AJM(NROW,NCOL) - (RSX-RSBX(M))/DSP_CX(NSPX)
          ELSE
!
!---        Equilibrium mass action
!
            NSE = IEQ_E(1,NEQXX)
!
!---        Loop over 1 species  ---
!
            DO MSX = 2,NSE
              NCM_SP = IEQ_E(MSX+1,NEQXX)
              NCOL = IEQ_S(NCM_SP)-NEQE
              AJM(NROW,NCOL) = AJM(NROW,NCOL)-(RSX-RSBX(M))/
     &          DSP_CX(NSPX)*EQ_E(MSX-1,NEQXX)*SP_CX(NSPX)/SP_CX(NCM_SP)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
!
!---  Loop over kinetic species in kinetic equation  ---
!
      DO M = 1,NS
        NSPX = IEQ_K(M+1,NEQX)
!
!---    Incremented species  ---
!
        NEQXX = IEQ_S(NSPX)
        IF( NEQXX.GT.NEQE ) THEN
          NCOL = IEQ_S(NSPX)-NEQE
          AJM(NROW,NCOL) = AJM(NROW,NCOL)+EQ_K(M,NEQX)*VTOLX*DTI
!
!---      Fixed species activity
!
          DO NSPKX = 1,NSPLK
            IF( ISPLK(14+NSPKX).LT.-1000 ) THEN
              NSPXX = ABS(ISPLK(14+NSPKX))-1000
              IF( NSPX.EQ.NSPXX ) THEN
                IF( NSPX.LE.NSPL ) THEN
                  ACX = ACTVX(1,NSPX)
                ELSE
                  ACX = 1.D0
                ENDIF
                AJM(NROW,NCOL)=AJM(NROW,NCOL)+EQ_K(M,NEQX)/ACX*VTOLX*DTI
              ENDIF
            ENDIF
           ENDDO
        ELSE
!
!---      Equilibrium mass action
!
          NSE = IEQ_E(1,NEQXX)
!
!---      Loop over 1 species  ---
!
          DO MSX = 2,NSE
            NCM_SP = IEQ_E(MSX+1,NEQXX)
            NCOL = IEQ_S(NCM_SP)-NEQE
            AJM(NROW,NCOL) = AJM(NROW,NCOL)+EQ_K(M,NEQX)*VTOLX*DTI*
     &        EQ_E(MSX-1,NEQXX)*SP_CX(NSPX)/SP_CX(NCM_SP)
          ENDDO
        ENDIF
      ENDDO
!
!---  Return residual vector and Jacobian matrix  ---
!
      IF( ABS(INDX).EQ.1 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      BJM(NROW) = -BJM(NROW)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KECHEM_R group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE MOBCF( NEQ )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
!
!     Mobile conservation component fractions.
!
!     YSPLX aqueous fraction of component species NEQ at node N
!     YSPGX gas fraction of component species NEQ at node N
!     YSPNX NAPL fraction of component species NEQ at node N
!     C(N,NSL) component species concentration (kmol/m^3 node)
!     SP_C(N,NSP) species concentration (kmol/m^3 node)
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 1 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
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
      SUB_LOG(ISUB_LOG) = '/MOBCF'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = ND(N)
        YSPLX = 0.D+0
        YSPGX = 0.D+0
        YSPNX = 0.D+0
        YSPLZ = 0.D+0
        YSPGZ = 0.D+0
        YSPNZ = 0.D+0
        NSL = NSOLU + NEQ
        IF( ICT(N,NSL).NE.0 .AND. NSTEP.EQ.0 ) CYCLE
        C(N,NSL) = 0.D+0
!
!---    Loop over conservation-component species  ---
!
        DO M = 1,IEQ_C(1,NEQ)
          NSP = IEQ_C(M+1,NEQ)
          IF( ABS(SP_C(N,NSP)).LT.CMIN ) THEN
            SP_CX = 0.D+0
          ELSE
            SP_CX = SP_C(N,NSP)
          ENDIF
!
!---      Aqueous species ---
!
          IF( NSP.LE.NSPL .AND. IEQW.GT.0 ) THEN
            YSPLX = YSPLX + EQ_C(M,NEQ)*SP_CX
            YSPLZ = YSPLZ + EQ_C(M,NEQ)
            C(N,NSL) = C(N,NSL) + EQ_C(M,NEQ)*SP_CX
!
!---      Gas species ---
!
          ELSEIF( NSP.GT.(NSPL+NSPS+NSPE) .AND.
     &      NSP.LE.(NSPL+NSPS+NSPE+NSPG) .AND. (IEQA.GT.0 .OR. 
     &      IOM.EQ.30 .OR. IOM.EQ.40 .OR. IOM.EQ.43) ) THEN
            YSPGX = YSPGX + EQ_C(M,NEQ)*SP_CX
            YSPGZ = YSPGZ + EQ_C(M,NEQ)
            C(N,NSL) = C(N,NSL) + EQ_C(M,NEQ)*SP_CX

          ENDIF
        ENDDO
!
!---    Update old time step conservation-component species ---
!
        CO(N,NSL) = C(N,NSL)
        YSPZ = YSPLZ+YSPGZ+YSPNZ
!
!---    Aqueous species ---
!
        IF( IEQW.GT.0 ) THEN
!
!---      Zero mobile species  ---
!
          IF( ABS(YSPZ)/EPSL.LT.EPSL ) THEN
            YSPLX = 0.D+0
!
!---      Zero species concentration  ---
!
          ELSEIF( ABS(C(N,NSL))/EPSL.LT.EPSL ) THEN
            YSPLX = YSPLZ/YSPZ
!            YSPLX = 0.D+0
!
!---      Non-zero species concentration  ---
!
          ELSE
            YSPLX = YSPLX/C(N,NSL)
          ENDIF
          YL(N,NSL) = YSPLX
!
!---      pH link ---
!
          IF( ISPLK(1).EQ.NSL ) THEN
            YL(N,NSL) = 1.D+0
!
!---      Air link ---
!
          ELSEIF( ISPLK(4).EQ.NSL ) THEN
            YL(N,NSL) = 1.D+0
!
!---      CO2 link ---
!
          ELSEIF( ISPLK(6).EQ.NSL ) THEN
            YL(N,NSL) = 1.D+0
          ENDIF
!
!---      Component link ---
!
          ISHIFTX=0
          IF( IOM.EQ.43 ) ISHIFTX=2
          DO IGC = 1,NGC+ISHIFTX
          IF( ISPLK(14+NSPLK+IGC).EQ.NSL ) THEN
            YL(N,NSL) = 1.D+0
          ENDIF
          ENDDO
        ENDIF
!
!---    Gas species ---
!
        IF( IEQA.GT.0 .OR. IOM.EQ.30 .OR. IOM.EQ.40 .OR. 
     &    IOM.EQ.43 ) THEN
!
!---      Zero mobile species  ---
!
          IF( ABS(YSPZ)/EPSL.LT.EPSL ) THEN
            YSPGX = 0.D+0
!
!---      Zero species concentration  ---
!
          ELSEIF( ABS(C(N,NSL))/EPSL.LT.EPSL ) THEN
!            YSPGX = 0.D+0
            YSPGX = YSPGZ/(YSPZ)
!
!---      Non-zero species concentration  ---
!
          ELSE
            YSPGX = YSPGX/C(N,NSL)
          ENDIF
          YG(N,NSL) = YSPGX
!
!---      Component link ---
!
          ISHIFTX=0
          IF( IOM.EQ.43 ) ISHIFTX=2
          DO IGC = 1,NGC+ISHIFTX
            IF( ISPLK(14+NSPLK+NGC+ISHIFTX+IGC).EQ.NSL ) THEN
              YG(N,NSL) = 1.D+0
            ENDIF
          ENDDO
        ENDIF

      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of MOBCF group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE MOBKF( NEQ )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
!
!     Mobile kinetic component fractions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 1 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
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
      SUB_LOG(ISUB_LOG) = '/MOBKF'
!
!---  Convert to global component indices  ---
!
      NEQX = NEQ + NEQC
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = ND(N)
        YSPLX = 0.D+0
        YSPGX = 0.D+0
        YSPNX = 0.D+0
        YSPLZ = 0.D+0
        YSPGZ = 0.D+0
        YSPNZ = 0.D+0
        NSL = NSOLU + NEQX
        IF( ICT(N,NSL).NE.0 .AND. NSTEP.EQ.0 ) CYCLE
        C(N,NSL) = 0.D+0
!
!---    Loop over kinetic-component species  ---
!
        DO M = 1,IEQ_K(1,NEQ)
          NSP = IEQ_K(M+1,NEQ)
          IF( ABS(SP_C(N,NSP)).LT.CMIN ) THEN
            SP_CX = 0.D+0
          ELSE
            SP_CX = SP_C(N,NSP)
          ENDIF
!
!---      Aqueous species ---
!
          IF( NSP.LE.NSPL ) THEN
            YSPLX = YSPLX + EQ_K(M,NEQ)*SP_CX
            YSPLZ = YSPLZ + EQ_K(M,NEQ)
            C(N,NSL) = C(N,NSL) + EQ_K(M,NEQ)*SP_CX
!
!---      Gas species ---
!
          ELSEIF( NSP.GT.(NSPL+NSPS+NSPE) .AND.
     &      NSP.LE.(NSPL+NSPS+NSPE+NSPG) ) THEN
            YSPGX = YSPGX + EQ_K(M,NEQ)*SP_CX
            YSPGZ = YSPGZ + EQ_K(M,NEQ)
            C(N,NSL) = C(N,NSL) + EQ_K(M,NEQ)*SP_CX
!
!---      NAPL species ---
!
          ELSEIF( NSP.GT.(NSPL+NSPS+NSPE+NSPG) .AND.
     &      NSP.LE.(NSPL+NSPS+NSPE+NSPG+NSPN) ) THEN
            YSPNX = YSPNX + EQ_K(M,NEQ)*SP_CX
            YSPNZ = YSPNZ + EQ_K(M,NEQ)
            C(N,NSL) = C(N,NSL) + EQ_K(M,NEQ)*SP_CX
          ENDIF
        ENDDO
!
!---    Update old time step kinetic-component species ---
!
        CO(N,NSL) = C(N,NSL)
        YSPZ = YSPLZ+YSPGZ+YSPNZ
!
!---    Aqueous species ---
!
        IF( IEQW.GT.0 ) THEN
!
!---      Zero mobile species ---
!
          IF( ABS(YSPZ)/EPSL.LT.EPSL ) THEN
            YSPLX = 0.D+0
!
!---      Zero species concentration  ---
!
          ELSEIF( ABS(C(N,NSL))/EPSL.LT.EPSL ) THEN
!            YSPLX = 0.D+0
            YSPLX = YSPLZ/YSPZ
!
!---      Non-zero species concentration  ---
!
          ELSE
            YSPLX = YSPLX/C(N,NSL)
          ENDIF
          YL(N,NSL) = YSPLX
        ENDIF
!
!---    Gas species ---
!
        IF( IEQA.GT.0 .OR. IOM.EQ.30 .OR. IOM.EQ.40 .OR. 
     &    IOM.EQ.43 ) THEN
!
!---      Zero mobile species ---
!
          IF( ABS(YSPZ)/EPSL.LT.EPSL ) THEN
            YSPGX = 0.D+0
!
!---      Zero species concentration  ---
!
          ELSEIF( ABS(C(N,NSL))/EPSL.LT.EPSL ) THEN
!            YSPGX = 0.D+0
            YSPGX = YSPGZ/YSPZ
!
!---      Non-zero species concentration  ---
!
          ELSE
            YSPGX = YSPGX/C(N,NSL)
          ENDIF
          YG(N,NSL) = YSPGX
        ENDIF

      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of MOBKF group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE NMNSP
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
!
!     Normalize mineral species.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 11 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
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
      SUB_LOG(ISUB_LOG) = '/NMNSP'
!
!---  Loop over active nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    Loop over solid species  ---
!
        DO NSPX = 1,NSPS
          NSP = NSPL + NSPX
!
!---      Mineral species  ---
!
          IF( ISP_MN(NSP).EQ.1 ) THEN
            IF( IEO.EQ.2 ) THEN
              SP_C(N,NSP) = SP_C(N,NSP) - SP_CMN(N,NSPX)
            ELSE
              SP_CMN(N,NSPX) = SP_C(N,NSP)
              SP_C(N,NSP) = 0.D+0
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of NMNSP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PITZER( ACTVX,SP_CX,DSP_CX,SLX,PORDX,RHOLX,TX,XLWX )
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
!     This subroutine computes activity and osmotic coefficients for
!     electrolyte solutions: Pitzer's equations for mixed 
!     electrolyte solutions.
!
!----------------------Authors-----------------------------------------!
!
!     ECKEChemX
!
!     Written by A. Felmy, from GMIN.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
      USE PTZRCOEF
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!

      REAL*8 ACTVX(LSPL+1,LSPL),SP_CX(LSPL),DSP_CX(LSPL)
      REAL*8 CPIX(LSPL+1),CAPZ(LSPL+1),LNG(LSPL+1,LSPL),TFUNC
      REAL*8 FF(LSPL+1),G4M(LSPL+1)
      EXTERNAL TFUNC
      REAL*8 LNA, G(4),GP(4),GPP(4),ALPHA(4),ATTMP(8)
!
!----------------------Data Statements---------------------------------!
!
!      SAVE TSMI
!      DATA TSMI / -1.D+3 /
      DATA ALPHA / 2.D+00,1.4D+00,1.2D+01,5.D+01 /
      DATA BB  /1.2D+00 /

!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PITZER'
!
!--- Initialize Variables
!
      DO I = 1,4
        G(I) = 0.D+0
        GP(I) = 0.D+0
        GPP(I) = 0.D+0
      ENDDO
      
      DO I = 1,8
        ATTMP(I) = 0.D+0
      ENDDO

      DO M = 1,NSPL+1
        CAPZ(M) = 0.D+0
        CPIX(M) = 0.D+0
          FF(M) = 0.D+0
          G4M(M) = 0.D+0
        DO NSP = 1,NSPL
          LNG(M,NSP) = 0.D+0
        ENDDO
      ENDDO
      
      IKH2O = 0

!
!---  Recalculate Pitzer parameters for non-isothermal solution.
!
      TK = TX + 273.15D+0

      APHI = 0.336901532D+0-6.3210043D-04*TK+9.14252359D+0/
     &       TK-1.35143986D-02*DLOG(TK)+2.26089488D-03/(TK-263.D+0)+
     &       1.92118597D-06*TK*TK+4.52586464D+01/(680.D+0-TK)

!  
!---  Re-calculate binary parameters
!
        DO I = 1,NCC_PZ
          DO J = 1,NA_PZ
            DO k = 1,8
              ATTMP(K) = ATB0(I,J,K)
            ENDDO
            B0(I,J) = TFUNC(ATTMP,TK)
          ENDDO
        ENDDO
        
        DO I = 1,NCC_PZ
          DO J = 1,NA_PZ
            DO K = 1,8
              ATTMP(K) = ATB1(I,J,K)
            ENDDO
            B1(I,J) = TFUNC(ATTMP,TK)
          ENDDO
        ENDDO
 
        DO I = 1,NCC_PZ
          DO J = 1,NA_PZ
            DO K = 1,8
             ATTMP(K) = ATB2(I,J,K)
            ENDDO
            B2(I,J) = TFUNC(ATTMP,TK)
          ENDDO
        ENDDO
  
        DO I = 1,NCC_PZ
          DO J = 1,NA_PZ
            DO K = 1,8
              ATTMP(K) = ATCMX(I,J,K)
            ENDDO
            CMXX(I,J) = TFUNC(ATTMP,TK)
            CMXX(I,J)=CMXX(I,J)/(2.*SQRT(ABS(SP_L(1,JPC(I))*
     &                SP_L(1,JPA(J)))))
          ENDDO
        ENDDO
!  
!---  Recalculate theta and psi.
!  
        II = 0
        IF( NCC_PZ.GE.2 )THEN
          DO I = 1,NCC_PZ
            DO J = I+1,NCC_PZ
              II = II+1
              DO K = 1,6
                ATTMP(K) = ATTC(II,K)
              ENDDO
              TCC(II)=TFUNC(ATTMP,TK)
              DO KK = 1,NA_PZ
                DO K = 1,6
                  ATTMP(K) = ATPC(II,KK,K)
                ENDDO
                PSIC(II,KK)=TFUNC(ATTMP,TK)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
 
        II = 0
        IF( NA_PZ.GE.2 )THEN
          DO I = 1,NA_PZ
            DO J = I+1,NA_PZ
              II = II+1
              DO K = 1,6
                ATTMP(K) = ATTA(II,K)
              ENDDO
              TAA(II)=TFUNC(ATTMP,TK)
              DO KK = 1,NCC_PZ
                DO K = 1,6
                  ATTMP(K) = ATPA(II,KK,K)
                ENDDO
                PSIA(II,KK)=TFUNC(ATTMP,TK)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!  
!--- Recalculate ternary parameters.
!  
        DO I = 1,NNN_PZ
          DO J = 1,NNN_PZ
            DO K = 1,6
              ATTMP(K) = ATNLAM(I,J,K)
            ENDDO
            ELAMB(I,J) = TFUNC(ATTMP,TK)
          ENDDO
          
          DO J = 1,NCC_PZ
            DO K = 1,6
              ATTMP(K) = ATCLAM(I,J,K)
            ENDDO
            CLAMB(I,J) = TFUNC(ATTMP,TK)
          ENDDO
          
          DO J = 1,NA_PZ
            DO K = 1,6
              ATTMP(K) = ATALAM(I,J,K)
            ENDDO
            ALAMB(I,J) = TFUNC(ATTMP,TK)
          ENDDO
        ENDDO
        
        II = 0
        DO I = 1,NNN_PZ
          DO J = 1,NCC_PZ
            II = II+1
            DO KK = 1,NA_PZ
              DO K = 1,6
               ATTMP(K) = ATHLAM(II,KK,K)
              ENDDO
              HOLAMB(II,KK) = TFUNC(ATTMP,TK)
            ENDDO
          ENDDO
        ENDDO
!
!--- End temperature-dependent calculation of Pitzer parameters
!
!
!---  Ionic strength of the aqueous solution  ---
!
      DO M = 1,NSPL+1
        CPIX(M) = 0.D+0
        SUM_M = 0.D+0
        DO NSP = 1,NSPL
          IF(IDD_PZ(NSP).GE.300000) CYCLE
          IF(SPNML(NSP).EQ.'h2o') THEN
            IKH2O = NSP
            CYCLE
          ENDIF
          IF( NSP.EQ.(M-1) ) THEN
            CLX = SP_CX(NSP) + DSP_CX(NSP)
          ELSE
            CLX = SP_CX(NSP)
          ENDIF
!       
!---      Molarity in mol solute/m^3 aqueous
!         or mol solute/l aqueous  ---
!       
          CMX = CLX/(SLX*PORDX)
!       
!---      Molality in mol solute/kg water  ---
!       
          CMX = CMX/(RHOLX*XLWX)
          SUM_M = SUM_M + CMX
          CPIX(M) = CPIX(M) + CMX*(SP_L(1,NSP)**2)
          CAPZ(M) = CAPZ(M) + CMX*DABS(SP_L(1,NSP))
        ENDDO
        CPIX(M) = CPIX(M)*5.D-1
      ENDDO

!
!--- Calculate Pitzer Activities
!
      DO M = 1,NSPL+1

        PHI1 = 0.D+0
        PHI2 = 0.D+0
        PHI3 = 0.D+0
        PHI4 = 0.D+0
        PHI5 = 0.D+0
        PHI6 = 0.D+0
        PHI7 = 0.D+0
        PHI8 = 0.D+0
        
        F1   = 0.D+0
        F2   = 0.D+0
        F3   = 0.D+0
        F4   = 0.D+0
        F5   = 0.D+0
        TMA  = 0.D+0
        TMCX = 0.D+0
        TMN  = 0.D+0
        FPR  = 0.D+0
!
!--- Calculate g functions
!
        DO I = 1,4
          X1 = ALPHA(I)*DSQRT(CPIX(M))
          X2 = X1*X1
          DEX = DEXP(-X1)
          G(I) = 2.D+0*(1.D+0-(1.D+0+X1)*DEX)/X2
          GP(I) = -2.D+0*(1.D+0-(1.D+0+X1+X2/2.D+0)*DEX)/X2
          GPP(I) = -(2.D+0*GP(I)+(X1/2.D+0)*DEX)/(CPIX(M)*CPIX(M))
        ENDDO
!
!--- Calculate b functions
!
        DO I = 1,NCC_PZ
          DO J = 1,NA_PZ
            K1 = 1
            K2 = 3
            IF( SP_L(1,JPC(I)).GE.2.D+0.AND.
     &        DABS(SP_L(1,(JPA(J)))).GE.2.D+0 )THEN
              K1 = 2
              K2 = 4
            ENDIF
            IF(SP_L(1,JPC(I)).EQ.2.D+0.AND.
     &        DABS(SP_L(1,JPA(J))).EQ.2.D+0 )THEN
              K1 = 2
              K2 = 3
            ENDIF
            X1 = -ALPHA(K1)*DSQRT(CPIX(M))
            X2 = -ALPHA(K2)*DSQRT(CPIX(M))
            BMMX(I,J) = B0(I,J)+B1(I,J)*G(K1)+B2(I,J)*G(K2)
            BPHI(I,J) = B0(I,J)+B1(I,J)*DEXP(X1)+B2(I,J)*DEXP(X2)
            BPR(I,J) = B1(I,J)*GP(K1)/CPIX(M)+B2(I,J)*GP(K2)/CPIX(M)
            BPPR(I,J) = B1(I,J)*GPP(K1)+B2(I,J)*GPP(K2)
          ENDDO
        ENDDO
!
!---    Calculate higher order mixing functions
!
        CALL HOMIX( CPIX(M),APHI )
!
!---    Start calculations for activity and osmotic coefficients
!
        TMP = 1.D+0+BB*DSQRT(CPIX(M))
        F1 = -APHI*DSQRT(CPIX(M))/TMP
        F2 = -APHI*(2.D+0/BB)*DLOG(TMP)
        PHI1 = F1*CPIX(M)
        FPR = -APHI*(TMP+0.5D+0)/(DSQRT(CPIX(M))*TMP*TMP)
!
!---    Start first major loop
!---    The terms are labeled in roughly the order they appear
!---    in the sums given by Felmy and Weare

        F3 = 0.D+0
        PHI2 = 0.D+0
        HBPP = 0.D+0

        DO I = 1,NCC_PZ
          DO J = 1,NA_PZ
            IF( JPC(I).EQ.(M-1) ) THEN
              CLXC = SP_CX(JPC(I)) + DSP_CX(JPC(I))
            ELSE
              CLXC = SP_CX(JPC(I))
            ENDIF
            IF( JPA(J).EQ.(M-1) ) THEN
              CLXA = SP_CX(JPA(J)) + DSP_CX(JPA(J))
            ELSE
              CLXA = SP_CX(JPA(J))
            ENDIF
            CMXC = CLXC/(SLX*PORDX)
            CMXC = CMXC/(RHOLX*XLWX)
            CMXA = CLXA/(SLX*PORDX)
            CMXA = CMXA/(RHOLX*XLWX)
            TMP = CMXC*CMXA
            F3 = F3+TMP*BPR(I,J)
            G4M(M) = G4M(M)+TMP*CMXX(I,J)
            PHI2 = PHI2+TMP*BPHI(I,J)
            HBPP = HBPP+TMP*BPPR(I,J)
            TMP1 = 2.D+0*BMMX(I,J)+CAPZ(M)*CMXX(I,J)
            LNG(M,JPC(I)) = LNG(M,JPC(I))+CMXA*TMP1
            LNG(M,JPA(J)) = LNG(M,JPA(J))+CMXC*TMP1
          ENDDO
        ENDDO

        
        PHI2 = PHI2+CAPZ(M)*G4M(M)

!
!---    Ternary electrolyte terms
!
        F4 = 0.D+0
        PHI3 = 0.D+0
        HTC = 0.D+0
        
        IF( NCC_PZ.GE.2 )THEN
          NT = 1
          DO I = 1,NCC_PZ-1
            DO J = I+1,NCC_PZ
              IF( JPC(I).EQ.(M-1) ) THEN
                CLXC = SP_CX(JPC(I)) + DSP_CX(JPC(I))
              ELSE
                CLXC = SP_CX(JPC(I))
              ENDIF
              IF( JPC(J).EQ.(M-1) ) THEN
                CLXC2 = SP_CX(JPC(J)) + DSP_CX(JPC(J))
              ELSE
                CLXC2 = SP_CX(JPC(J))
              ENDIF
              CMXC = CLXC/(SLX*PORDX)
              CMXC = CMXC/(RHOLX*XLWX)
              CMXC2 = CLXC2/(SLX*PORDX)
              CMXC2 = CMXC2/(RHOLX*XLWX)
              TMP = CMXC*CMXC2
              F4 = F4+TMP*CTCPR(NT)
              PHI3 = PHI3+TMP*CTCPH(NT)
              HTC = HTC+TMP*CTCPPR(NT)
!      
!---          Now g2m
!      
              DO K = 1,NCC_PZ
                N1 = 0
                IF( JPC(I).EQ.JPC(K) ) N1=J
                IF( JPC(J).EQ.JPC(K) ) N1=I
                IF( N1.NE.0 )THEN
                  IF( JPC(N1).EQ.(M-1) ) THEN
                    CLXC = SP_CX(JPC(N1)) + DSP_CX(JPC(N1))
                  ELSE
                    CLXC = SP_CX(JPC(N1))
                  ENDIF
                  CMXC = CLXC/(SLX*PORDX)
                  CMXC = CMXC/(RHOLX*XLWX)
                  TMP1 = CMXC
                  LNG(M,JPC(K)) = LNG(M,JPC(K))+2.D+0*TMP1*CTC(NT)
                  DO N = 1,NA_PZ
                    IF( JPA(N).EQ.(M-1) ) THEN
                      CLXA = SP_CX(JPA(N)) + DSP_CX(JPA(N))
                    ELSE
                      CLXA = SP_CX(JPA(N))
                    ENDIF
                    CMXA = CLXA/(SLX*PORDX)
                    CMXA = CMXA/(RHOLX*XLWX)
                    LNG(M,JPC(K)) = LNG(M,JPC(K))+TMP1*CMXA*PSIC(NT,N)
                  ENDDO
                ENDIF
              ENDDO
             
              DO K = 1,NA_PZ
                TMP1 = TMP*PSIC(NT,K)
                LNG(M,JPA(K)) = LNG(M,JPA(K))+TMP1
                IF( JPA(K).EQ.(M-1) ) THEN
                  CLXA = SP_CX(JPA(K)) + DSP_CX(JPA(K))
                ELSE
                  CLXA = SP_CX(JPA(K))
                ENDIF
                CMXA = CLXA/(SLX*PORDX)
                CMXA = CMXA/(RHOLX*XLWX)
                PHI3 = PHI3+TMP1*CMXA
              ENDDO
              NT=NT+1
            ENDDO
          ENDDO
        ENDIF

        F5 = 0.D+0
        PHI4 = 0.D+0
        HTA = 0.D+0

        IF( NA_PZ.GE.2 )THEN
          NT = 1
          DO I = 1,NA_PZ-1
            DO J = I+1,NA_PZ
              IF( JPA(I).EQ.(M-1) ) THEN
                CLXA = SP_CX(JPA(I)) + DSP_CX(JPA(I))
              ELSE
                CLXA = SP_CX(JPA(I))
              ENDIF
              IF( JPA(J).EQ.(M-1) ) THEN
                CLXA2 = SP_CX(JPA(J)) + DSP_CX(JPA(J))
              ELSE
                CLXA2 = SP_CX(JPA(J))
              ENDIF
              CMXA = CLXA/(SLX*PORDX)
              CMXA = CMXA/(RHOLX*XLWX)
              CMXA2 = CLXA2/(SLX*PORDX)
              CMXA2 = CMXA2/(RHOLX*XLWX)
              TMP = CMXA*CMXA2
              F5 = F5+TMP*CTAPR(NT)
              PHI4 = PHI4+TMP*CTAPH(NT)
              HTA = HTA+TMP*CTAPPR(NT)
!      
!---          Now g2x
!      
              DO K = 1,NA_PZ
                N1 = 0
                IF( JPA(I).EQ.JPA(K) )N1=J
                IF( JPA(J).EQ.JPA(K) )N1=I
                IF( N1.NE.0 )THEN
                  IF( JPA(N1).EQ.(M-1) ) THEN
                    CLXA = SP_CX(JPA(N1)) + DSP_CX(JPA(N1))
                  ELSE
                    CLXA = SP_CX(JPA(N1))
                  ENDIF
                  CMXA = CLXA/(SLX*PORDX)
                  CMXA = CMXA/(RHOLX*XLWX)
                  TMP1 = CMXA
                  LNG(M,JPA(K)) = LNG(M,JPA(K))+2.D+0*TMP1*CTA(NT)
                  DO N = 1,NCC_PZ
                    IF( JPC(N).EQ.(M-1) ) THEN
                      CLXC = SP_CX(JPC(N)) + DSP_CX(JPC(N))
                    ELSE
                      CLXC = SP_CX(JPC(N))
                    ENDIF
                    CMXC = CLXC/(SLX*PORDX)
                    CMXC = CMXC/(RHOLX*XLWX)
                    LNG(M,JPA(K)) = LNG(M,JPA(K))+TMP1*CMXC*PSIA(NT,N)
                  ENDDO
                ENDIF
              ENDDO
          
              DO K = 1,NCC_PZ
                TMP1 = TMP*PSIA(NT,K)
                LNG(M,JPC(K)) = LNG(M,JPC(K))+TMP1
                IF( JPC(K).EQ.(M-1) ) THEN
                  CLXC = SP_CX(JPC(K)) + DSP_CX(JPC(K))
                ELSE
                  CLXC = SP_CX(JPC(K))
                ENDIF
                CMXC = CLXC/(SLX*PORDX)
                CMXC = CMXC/(RHOLX*XLWX)
                PHI4 = PHI4+TMP1*CMXC
              ENDDO
              NT = NT+1
            ENDDO
          ENDDO
        ENDIF


        PHI5 = 0.D+0
        PHI6 = 0.D+0
        PHI7 = 0.D+0
        PHI8 = 0.D+0
        
        IF( NNN_PZ.GT.0 )THEN
          NT = 1
          DO I = 1,NNN_PZ
            IF( JPN(I).EQ.(M-1) ) THEN
              CLXN = SP_CX(JPN(I)) + DSP_CX(JPN(I))
            ELSE
              CLXN = SP_CX(JPN(I))
            ENDIF
            CMXN = CLXN/(SLX*PORDX)
            CMXN = CMXN/(RHOLX*XLWX)
            TMN = CMXN
            DO J = 1,NCC_PZ
              IF( JPC(J).EQ.(M-1) ) THEN
                CLXC = SP_CX(JPC(J)) + DSP_CX(JPC(J))
              ELSE
                CLXC = SP_CX(JPC(J))
              ENDIF
              CMXC = CLXC/(SLX*PORDX)
              CMXC = CMXC/(RHOLX*XLWX)
              TMCX = CMXC
              PHI5 = PHI5+TMN*TMCX*CLAMB(I,J)
              LNG(M,JPC(J)) = LNG(M,JPC(J))+2.D+0*TMN*CLAMB(I,J)
              LNG(M,JPN(I)) = LNG(M,JPN(I))+2.D+0*TMCX*CLAMB(I,J)
              DO K = 1,NA_PZ
                IF( JPA(K).EQ.(M-1) ) THEN
                  CLXA = SP_CX(JPA(K)) + DSP_CX(JPA(K))
                ELSE
                  CLXA = SP_CX(JPA(K))
                ENDIF
                CMXA = CLXA/(SLX*PORDX)
                CMXA = CMXA/(RHOLX*XLWX)
                TMA = CMXA
                PHI7 = PHI7+TMN*TMCX*TMA*HOLAMB(NT,K)
                LNG(M,JPN(I)) = LNG(M,JPN(I))+TMCX*TMA*HOLAMB(NT,K)
                LNG(M,JPC(J)) = LNG(M,JPC(J))+TMN*TMA*HOLAMB(NT,K)
                LNG(M,JPA(K)) = LNG(M,JPA(K))+TMN*TMCX*HOLAMB(NT,K)
              ENDDO
              NT = NT+1
            ENDDO
            DO K = 1,NA_PZ
              LNG(M,JPA(K)) = LNG(M,JPA(K))+2.D+0*TMN*ALAMB(I,K)
              IF( JPA(K).EQ.(M-1) ) THEN
                CLXA = SP_CX(JPA(K)) + DSP_CX(JPA(K))
              ELSE
                CLXA = SP_CX(JPA(K))
              ENDIF
              CMXA = CLXA/(SLX*PORDX)
              CMXA = CMXA/(RHOLX*XLWX)
              LNG(M,JPN(I)) = LNG(M,JPN(I))+2.D+0*CMXA*ALAMB(I,K)
              PHI6 = PHI6+TMN*CMXA*ALAMB(I,K)
            ENDDO
            DO K = 1,NNN_PZ
              IF( JPN(K).EQ.(M-1) ) THEN
                CLXN = SP_CX(JPN(K)) + DSP_CX(JPN(K))
              ELSE
                CLXN = SP_CX(JPN(K))
              ENDIF
              CMXN = CLXN/(SLX*PORDX)
              CMXN = CMXN/(RHOLX*XLWX)
              PHI8 = PHI8+CMXN*TMN*ELAMB(I,K)
              LNG(M,JPN(I)) = LNG(M,JPN(I))+2.D+0*CMXN*ELAMB(I,K)
            ENDDO
          ENDDO
        ENDIF
!
!---    Sum up the f function and carry along the proper charges
!
        FF(M) = F1+F2+F3+F4+F5
!
!---    Calculate osmotic coefficient
!
        IF(SUM_M.NE.0)PHI = (2.D+0/SUM_M)*
     &                      (PHI1+PHI2+PHI3+PHI4+PHI5+PHI6+PHI7+PHI8)
        PHI = PHI+1.D+0
        
        DO NSP = 1,NSPL
          IF( NSP.EQ.IKH2O ) THEN
            IF( IKH2O.EQ.(M-1) ) THEN
              CLX = SP_CX(IKH2O) + DSP_CX(IKH2O)
            ELSE
              CLX = SP_CX(IKH2O)
            ENDIF
            CMX = CLX/(SLX*PORDX)
            CMX = CMX/(RHOLX*XLWX)
            LNA = -.0180153*SUM_M*PHI
            H2OACT = DEXP(LNA)
            LNG(M,IKH2O) = LNA - DLOG(CMX)
            ACTVX(M,IKH2O) = DEXP(LNG(M,IKH2O))
          ELSE
            LNG(M,NSP) = LNG(M,NSP)+SP_L(1,NSP)*SP_L(1,NSP)*FF(M)+
     &                 DABS(SP_L(1,NSP))*G4M(M)
            ACTVX(M,NSP) = DEXP(LNG(M,NSP))
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PITZER ---
!
      RETURN
      END 

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RATER( NEQX,M,NS,N,IRCX,NSPKX,NSPRX,NSPPX,SP_CX,ACTVX,
     &      VTOMX,VTOLX,VTOSX,EQKX,RRCX,RSBX )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
!
!     Kinetic reaction rates
!
!     YSPLX aqueous fraction of component species NEQ at node N
!     YSPGX gas fraction of component species NEQ at node N
!     YSPNX NAPL fraction of component species NEQ at node N
!     C(N,NSL) component species concentration (kmol/m^3 node)
!     SP_C(N,NSP) species concentration (kmol/m^3 node)
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 3 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
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
      REAL*8 ACTVX(LSPL+1,LSPL), SP_CX(LSPR)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RATER'
!
!---    TST kinetic reaction  ---
!
        IF( (IRCKT(IRCX).GE.10 .AND. IRCKT(IRCX).LE.12) .OR.
     &    (IRCKT(IRCX).GE.5 .AND. IRCKT(IRCX).LE.9) ) THEN
!
!---      Ion activity product mol/kg water, loop over species in
!         kinetic reaction  ---
!
          QX = 1.D+0
          DO L = 1,NSPKX
            NSPX = IRC_K(L+2,IRCX)
!
!---        Aqueous species,
!           concentration in molality, mol/kg H2O  ---
!
            IF( NSPX.LE.NSPL ) THEN
              CMX = SP_CX(NSPX)*VTOMX
              ACX = ACTVX(1,NSPX)
!
!---          Reactants  ---
!
              N1 = MAX(1,N*IRCKN(L))
              IF( L.LE.NSPRX ) THEN
                QX = QX*((CMX*ACX)**RC_K(L,N1,IRCX))
!
!---          Products  ---
!
              ELSE
                QX = QX/((CMX*ACX)**RC_K(L,N1,IRCX))
              ENDIF
!
!---        Solid species, skip  ---
!
            ELSEIF( NSPX.LE.NSPL+NSPS ) THEN
              CYCLE
            ENDIF
          ENDDO
!
!---      Initial reactive surface area, initial mineral volume 
!         fraction, current mineral volume fraction, minimum current 
!         mineral volumefraction allows re-precipitation of dissolved 
!         primary minerals NSP_M - mineral species number  ---
!
          NSPX = IRC_K(3+NSPKX,IRCX)
          NSP_M =  NSPX - NSPL
!
!---      Primary mineral  ---
!
          IF( RS_S(2,NSP_M,N).GT.EPSL ) THEN
            AOX = RS_S(1,NSP_M,N)*VOL(N)*RS_S(2,NSP_M,N)*SP_S(1,NSP_M)
            VFMOX = RS_S(2,NSP_M,N)
            IF( ISP_MN(NSPX).EQ.1 ) THEN
              VFMX = 1.D-3*(SP_CMN(N,NSP_M)+SP_CX(NSPX))
     &          *SP_S(2,NSP_M)/SP_S(1,NSP_M)
            ELSE
              VFMX = 1.D-3*SP_CX(NSPX)*SP_S(2,NSP_M)/SP_S(1,NSP_M)
            ENDIF
            VFMX = MAX( VFMX,1.D-2 )
!
!---      Secondary mineral, initial reactive surface area
!         for seconary minerals is set to 0.25 m^2/dm^3  ---
!
          ELSE
            AOX = 0.25D+3*VOL(N)
            VFMOX = 1.D-2
            VFMX = 1.D-3*SP_CX(NSPX)*SP_S(2,NSP_M)/SP_S(1,NSP_M)
            VFMX = MAX( VFMX,1.D-2 )
          ENDIF
          VFMX = MAX( VFMX,0.D+0 )
!
!---      Reactive surface area  ---
!
          AX = AOX*(((POR0(2,N)*VFMX)/
     &      (POR(2,N)*VFMOX))**(2.D+0/3.D+0))
          IF( ISLC(56).EQ.2 ) AX = AOX
!
!---      Reaction rate, mol/s  ---
!
          RRBX = -AX*RRCX*(1.D+0-(QX/EQKX))
!
!---      pH dependence  ---
!
          IF( IRCKT(IRCX).GE.5 .AND. IRCKT(IRCX).LE.9
     &      .AND. ISPLK(1).NE.0 ) THEN
            NSP_PHX = MOD(ISPLK(1),1000)
            PHX = -LOG10(1.D-3*SP_CX(NSP_PHX)*VTOLX)
            IF( IRCKT(IRCX).GE.8 .AND. IRCKT(IRCX).LE.9 ) THEN
              RRBX = RRBX*MAX( 0.D+0,
     &          (7.9201D-1 - 1.3479D-1*PHX + 5.2D-3*(PHX**2)))
            ELSE
              N9 = MAX(1,N*IRCKN(NSPKX+9))
              RRBX = RRBX*(1.D+1**(-RC_K(NSPKX+9,N9,IRCX)*PHX))
            ENDIF
          ENDIF
!
!---      Reaction rate, mol/m^3 aqu s  ---
!
          RRBX = RRBX*VTOLX/VOL(N)
!
!---      Direction limited  ---
!
          IF( IRCKT(IRCX).EQ.6 .OR. IRCKT(IRCX).EQ.8
     &      .OR. IRCKT(IRCX).EQ.11 ) THEN
            RRBX = MAX( RRBX,0.D+0 )
          ELSEIF( IRCKT(IRCX).EQ.7 .OR. IRCKT(IRCX).EQ.9
     &      .OR. IRCKT(IRCX).EQ.12 ) THEN
            RRBX = MIN( RRBX,0.D+0 )
          ENDIF
!
!---    Forward-backward kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.1 ) THEN
          N1 = MAX(1,N*IRCKN(NSPKX+1))
          N2 = MAX(1,N*IRCKN(NSPKX+2))
          FRRX = RC_K(NSPKX+1,N1,IRCX)
          BRRX = RC_K(NSPKX+2,N2,IRCX)
!
!---      Loop over reactants  ---
!
          DO L = 1,NSPRX
!
!---        Global species index  ---
!
            NSPX = IRC_K(L+2,IRCX)
            CMX = SP_CX(NSPX)*VTOLX
            N1 = MAX(1,N*IRCKN(L))
            FRRX = FRRX*(CMX**RC_K(L,N1,IRCX))
          ENDDO
!
!---      Loop over products  ---
!
          DO L = 1,NSPPX
!
!---        Global species index  ---
!
            NSPX = IRC_K(L+2+NSPRX,IRCX)
            CMX = SP_CX(NSPX)*VTOLX
            N1 = MAX(1,N*IRCKN(L+NSPRX))
            BRRX = BRRX*(CMX**RC_K(L+NSPRX,N1,IRCX))
          ENDDO
          RRBX = FRRX - BRRX
!
!---    Valocchi-Monod kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.2 ) THEN
!
!---      Concentration of biomass in mol/m^3 aqu  ---
!
          NSPX = IRC_K(5,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
!
!---      Partial rate of donor degradation  ---
!
          N6 = MAX(1,N*IRCKN(6))
          RRBX = RC_K(6,N6,IRCX)*CMX
!
!---      Concentration of donor in mol/kg water  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Partial rate of donor degradation  ---
!
          N4 = MAX(1,N*IRCKN(4))
          RRBX = RRBX*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---      Concentration of acceptor in mol/kg water  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Rate of donor degradation  ---
!
          N5 = MAX(1,N*IRCKN(5))
          RRBX = RRBX*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---    Schroth-Monod kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.17 ) THEN
!
!---      Concentration of biomass in mol/m^3 aqu  ---
!
          NSPX = IRC_K(5,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
!
!---      Partial rate of donor degradation  ---
!
          N6 = MAX(1,N*IRCKN(6))
          RRBX = RC_K(6,N6,IRCX)*CMX
!
!---      Concentration of donor in mol/kg water  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Partial rate of donor degradation  ---
!
          N4 = MAX(1,N*IRCKN(4))
          RRBX = RRBX*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---      Concentration of acceptor in mol/kg water  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Rate of donor degradation  ---
!
          N5 = MAX(1,N*IRCKN(5))
          RRBX = RRBX*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---      Summation of limiter concentrations  ---
!
          SCMX = 0.D+0
          DO L = 1,IRC_K(5,IRCX)
            NSPX = IRC_K(5+L,IRCX)
            SCMX = SCMX + SP_CX(NSPX)*VTOMX
          ENDDO
!
!---      Rate of total limiter control  ---
!
          N6 = MAX(1,N*IRCKN(6))
          RRBX = RRBX*(1.D+0-(SCMX/(RC_K(6,N6,IRCX)+SCMX)))
!
!---    Valocchi-Sorption kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.3 ) THEN
!
!---      Concentration of sorbed species in mol/gm soil  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOSX*1.D-3
          N4 = MAX(1,N*IRCKN(4))
          RRBX = CMX/RC_K(4,N4,IRCX)
!
!---      Concentration of aqueous species in mol/kg water  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Rate of sorption in mol/m^3 aqu s  ---
!
          N3 = MAX(1,N*IRCKN(3))
          RRBX = RC_K(3,N3,IRCX)*(CMX-RRBX)*1.D+3
!
!---    Valocchi-Biomass kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.4 ) THEN
!
!---      Concentration of biomass in mol/m^3 aqu  ---
!
          NSPX = IRC_K(5,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
!
!---      Partial rate of donor degradation  ---
!
          N6 = MAX(1,N*IRCKN(6))
          RRBX = RC_K(6,N6,IRCX)*CMX
!
!---      Concentration of donor in mol/m^3 aqu  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Partial rate of donor degradation  ---
!
          N4 = MAX(1,N*IRCKN(4))
          RRBX = RRBX*(CMX/(RC_K(4,N4,IRCX)+CMX))
!
!---      Concentration of acceptor in mol/m^3 aqu  ---
!
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOMX
!
!---      Rate of donor degradation  ---
!
          N5 = MAX(1,N*IRCKN(5))
          RRBX = RRBX*(CMX/(RC_K(5,N5,IRCX)+CMX))
!
!---      Concentration of biomass in mol/m^3 aqu  ---
!
          NSPX = IRC_K(5,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
!
!---      Rate of biomass production in mol/m^3 aqu s  ---
!
          N7 = MAX(1,N*IRCKN(7))
          N8 = MAX(1,N*IRCKN(8))
          RRBX = RC_K(7,N7,IRCX)*RRBX - RC_K(8,N8,IRCX)*CMX
!
!---    Liu's multi-rate kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.41 ) THEN
          N1 = MAX(1,N*IRCKN(NSPKX+1))
          RMX = RC_K(NSPKX+1,N1,IRCX)
!
!---      Loop over reactants  ---
!
          N2 = MAX(1,N*IRCKN(NSPKX+2))
          N3 = MAX(1,N*IRCKN(NSPKX+3))
          N4 = MAX(1,N*IRCKN(NSPKX+4))
          N5 = MAX(1,N*IRCKN(NSPKX+5))
          SDENX = RC_K(NSPKX+2,N2,IRCX)*VTOLX/VTOSX
          PFRCX = RC_K(NSPKX+3,N3,IRCX)
          XLGK1 = RC_K(NSPKX+4,N4,IRCX)
          XLGK2 = RC_K(NSPKX+5,N5,IRCX)
          FRRX = 1.D0
          BRRX = 1.D0
          DO L = 1,NSPRX
!
!---        Global species index  ---
!
            NSPX = IRC_K(L+2,IRCX)
            CMX = SP_CX(NSPX)*VTOMX
            IF(NSPX.LE.NSPL) THEN
              ACX = ACTVX(1,NSPX)
            ELSE
              ACX =1.D0
            ENDIF
            N1 = MAX(1,N*IRCKN(L))
            FRRX = FRRX*((CMX*ACX)**RC_K(L,N1,IRCX))
          ENDDO
!
!---      Loop over products  ---
!
          DO L = 1,NSPPX-1
!
!---        Global species index  ---
!
            NSPX = IRC_K(L+2+NSPRX,IRCX)
            CMX = SP_CX(NSPX)*VTOMX
            IF(NSPX.LE.NSPL) THEN
              ACX = ACTVX(1,NSPX)
            ELSE
              ACX = 1.D0
            ENDIF
            N1 = MAX(1,N*IRCKN(L+NSPRX))
            BRRX = BRRX*((CMX*ACX)**RC_K(L+NSPRX,N1,IRCX))
          ENDDO
          L=NSPPX
          NSPX = IRC_K(L+2+NSPRX,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
          FRRX = FRRX*10**XLGK1
          BRRX = BRRX*10**XLGK2
          RRBX = RMX*(SDENX*PFRCX*FRRX/(1.D0+FRRX+BRRX)-CMX)
!
!---    Liu's dual domain kinetic reaction  ---
!
        ELSEIF( IRCKT(IRCX).EQ.42 ) THEN
          N1 = MAX(1,N*IRCKN(NSPKX+1))
          RMX = RC_K(NSPKX+1,N1,IRCX)
!
!---      Loop over reactants  ---
!
          N3 = MAX(1,N*IRCKN(NSPKX+3))
          PFRCX = RC_K(NSPKX+3,N3,IRCX)
          FRRX = 0.D0
!
!---      Global species index  ---
!
          NSPX = IRC_K(3,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
          FRRX = FRRX+CMX**VTOLX
          NSPX = IRC_K(4,IRCX)
          CMX = SP_CX(NSPX)*VTOLX
          FRRX = FRRX-CMX**VTOLX
          RRBX = RMX*FRRX
!
!---    Multirate  ---
!
        ELSEIF( IRCKT(IRCX).EQ.20 ) THEN
!
!---      Neutral reaction rate, mol/m^2 s  ---
!
          TKX = T(2,N)+TABS
          N1 = MAX(1,N*IRCKN(1))
          N2 = MAX(1,N*IRCKN(2))
          N3 = MAX(1,N*IRCKN(3))
          TKRX = RC_K(3,N3,IRCX)+TABS
          RRCX = RC_K(1,N1,IRCX)*EXP( -RC_K(2,N2,IRCX)*
     &      ((1.D+0/TKX)-(1.D+0/TKRX))/(1.D-3*RCU) )
!
!---      Loop over mechanisms  ---
!
          DO NKRMX = 1,IRC_K(2,IRCX)
            IX = 3+((NKRMX-1)*6)
!
!---        Ion activity product mol/kg water, loop over species  ---
!
            QX = 1.D+0
            DO L = 1,IRC_K(IX,IRCX)
              IX = 3+((NKRMX-1)*6)+L
              NSPX = IRC_K(IX,IRCX)
!
!---          Aqueous species,
!             concentration in molality, mol/kg H2O  ---
!
              IF( NSPX.LE.NSPL ) THEN
                CMX = SP_CX(NSPX)*VTOMX
                ACX = ACTVX(1,NSPX)
                IX = 6+((NKRMX-1)*8)+L
                NX = MAX(1,N*IRCKN(IX))
                QX = QX*((CMX*ACX)**RC_K(IX,NX,IRCX))
              ENDIF
            ENDDO
            IX = 4+((NKRMX-1)*8)
            NX = MAX(1,N*IRCKN(IX))
            N1 = MAX(1,N*IRCKN(IX+1))
            N2 = MAX(1,N*IRCKN(IX+2))
            TKRX = RC_K(IX+2,N2,IRCX)+TABS
            RRCX = RRCX + RC_K(IX,NX,IRCX)*EXP( -RC_K(IX+1,N1,IRCX)*
     &        ((1.D+0/TKX)-(1.D+0/TKRX))/(1.D-3*RCU) )*QX
          ENDDO
!
!---      Initial reactive surface area, initial mineral volume 
!         fraction,current mineral volume fraction, minimum current 
!         mineral volumefraction allows re-precipitation of dissolved 
!         primary minerals NSP_M - mineral species number  ---
!
          NSPX = IRC_K(1,IRCX)
          NSP_M =  NSPX - NSPL
!
!---      Primary mineral  ---
!
          IF( RS_S(2,NSP_M,N).GT.EPSL ) THEN
            AOX = RS_S(1,NSP_M,N)*VOL(N)*RS_S(2,NSP_M,N)*SP_S(1,NSP_M)
            VFMOX = RS_S(2,NSP_M,N)
            IF( ISP_MN(NSPX).EQ.1 ) THEN
              VFMX = 1.D-3*(SP_CMN(N,NSP_M)+SP_CX(NSPX))
     &          *SP_S(2,NSP_M)/SP_S(1,NSP_M)
            ELSE
              VFMX = 1.D-3*SP_CX(NSPX)*SP_S(2,NSP_M)/SP_S(1,NSP_M)
            ENDIF
            VFMX = MAX( VFMX,1.D-2 )
!
!---      Secondary mineral, initial reactive surface area
!         for seconary minerals is set to 0.25 m^2/dm^3  ---
!
          ELSE
            AOX = 0.25D+3*VOL(N)
            VFMOX = 1.D-2
            VFMX = 1.D-3*SP_CX(NSPX)*SP_S(2,NSP_M)/SP_S(1,NSP_M)
            VFMX = MAX( VFMX,1.D-2 )
          ENDIF
          VFMX = MAX( VFMX,0.D+0 )
!
!---      Reactive surface area  ---
!
          AX = AOX*(((POR0(2,N)*VFMX)/
     &      (POR(2,N)*VFMOX))**(2.D+0/3.D+0))
          IF( ISLC(56).EQ.2 ) AX = AOX
!
!---      Reaction rate, mol/s  ---
!
          RRBX = -AX*RRCX
!
!---      Reaction rate, mol/m^3 aqu s  ---
!
          RRBX = RRBX*VTOLX/VOL(N)
        ENDIF
!
!---    Liu's dual domain model  ---
!
        IF( IRCKT(IRCX).EQ.42 ) THEN
          EQ_KX = EQ_K(NS+M,NEQX)
          IF( IMMB(NEQC+NEQX) == 1 ) THEN
!
!---        RC_K(2,IRCX) is pore fraction of immobile domain  ---
!
            N2 = MAX(1,N*IRCKN(NSPKX+2))
            EQ_KX = EQ_K(NS+M,NEQX)*RC_K(NSPKX+2,N2,IRCX)
          ENDIF
          RSBX = RSBX + RRBX*EQ_KX
        ELSE
          RSBX = RSBX + RRBX*EQ_K(NS+M,NEQX)
        ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RATER group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE READ_PTZRCOEF
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
!     ECKEChemX
!
!     Read binary ptcf.bin file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 February 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE TRNSPT
      USE SOLTN
      USE PTZRCOEF
      USE PROP
      USE HYST
      USE GRID
      USE GLB_PAR
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
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/READ_PTZRCOEF'
!
!---  Allocate and initialize memory, global and local, for 
!     Pitzer coefficient arrays  ---
!
      CALL ALLOC_PTZRCOEF
      CALL INTLZ_PTZRCOEF
!
!---  Open ptcf.bin file for reactive species modeling  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'ptcf.bin',MPI_MODE_RDONLY,
     &  MPI_INFO_NULL,IFILE,IERR )
!
!---  Initialize cummulative offset  ---
!
      IOFFSET = 0
!
!---  Read Pitzer coefficient integers
!     (duplicated across processors)  ---
!
      NVAR = 3
      ALLOCATE( IVARX(1:NVAR),STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,IVARX,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
      NCC_PZ = IVARX(1)
      NA_PZ = IVARX(2)
      NNN_PZ = IVARX(3)
!      PRINT *,'IVARX = ',(IVARX(M),M=1,NVAR),' ID = ',ID
      DEALLOCATE( IVARX,STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,JPA,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!      PRINT *,'JPA = ',(JPA(L),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,JPC,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!      PRINT *,'JPC = ',(JPC(L),L=1,LCAT),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNEU
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,JPN,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!      PRINT *,'JPN = ',(JPN(L),L=1,LNEU),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,B0,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'B0 = ',((B0(K,L),K=1,LCAT),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,B1,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'B1 = ',((B1(K,L),K=1,LCAT),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,B2,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'B2 = ',((B2(K,L),K=1,LCAT),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,CMXX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'CMXX = ',((CMXX(K,L),K=1,LCAT),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNCF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,TCC,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'TCC = ',(TCC(L),L=1,LNCF),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNAF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,TAA,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'TAA = ',(TAA(L),L=1,LNAF),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNCF*LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,PSIC,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'PSIC = ',((PSIC(K,L),K=1,LNCF),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNAF*LCAT
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,PSIA,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'PSIA = ',((PSIA(K,L),K=1,LNAF),L=1,LCAT),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNEU*LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ALAMB,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ALAMB = ',((ALAMB(K,L),K=1,LNEU),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNEU*LCAT
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,CLAMB,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'CLAMB = ',((CLAMB(K,L),K=1,LNEU),L=1,LCAT),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNEU*LNEU
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ELAMB,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ELAMB = ',((ELAMB(K,L),K=1,LNEU),L=1,LNEU),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNNF*LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,HOLAMB,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'HOLAMB = ',((HOLAMB(K,L),K=1,LNNF),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,BPPR,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'BPPR = ',((BPPR(K,L),K=1,LCAT),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,BPHI,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'BPHI = ',((BPHI(K,L),K=1,LCAT),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,BPR,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'BPR = ',((BPR(K,L),K=1,LCAT),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,BMMX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'BMMX = ',((BMMX(K,L),K=1,LCAT),L=1,LANI),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI*8
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATB0,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATB0 = ',(((ATB0(K,L,M),K=1,LCAT),L=1,LANI),M=1,8),
!     &  ' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI*8
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATB1,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATB1 = ',(((ATB1(K,L,M),K=1,LCAT),L=1,LANI),M=1,8),
!     &  ' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI*8
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATB2,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATB2 = ',(((ATB2(K,L,M),K=1,LCAT),L=1,LANI),M=1,8),
!     &  ' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LCAT*LANI*8
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATCMX,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATCMX = ',(((ATCMX(K,L,M),K=1,LCAT),L=1,LANI),M=1,8),
!     &  ' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNEU*LNEU*6
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATNLAM,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATNLAM = ',(((ATNLAM(K,L,M),K=1,LNEU),L=1,LNEU),M=1,6),
!     &  ' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNEU*LCAT*6
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATCLAM,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATCLAM = ',(((ATCLAM(K,L,M),K=1,LNEU),L=1,LCAT),M=1,6),
!     &  ' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNEU*LANI*6
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATALAM,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATALAM = ',(((ATALAM(K,L,M),K=1,LNEU),L=1,LANI),M=1,6),
!     &  ' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNNF*LANI*6
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATHLAM,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATHLAM = ',(((ATHLAM(K,L,M),K=1,LNNF),L=1,LANI),M=1,6),
!     &  ' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNCF*6
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATTC,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATTC = ',((ATTC(K,L),K=1,LNCF),L=1,6),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNCF*LANI*6
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATPC,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATPC = ',(((ATPC(K,L,M),K=1,LNCF),L=1,LANI),M=1,6),
!     &  ' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNAF*6
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATTA,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATTA = ',((ATTA(K,L),K=1,LNAF),L=1,6),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNAF*LCAT*6
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ATPA,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ATPA = ',(((ATPA(K,L,M),K=1,LNAF),L=1,LCAT),M=1,6),
!     &  ' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNCF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,CTCPH,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'CTCPH = ',(CTCPH(L),L=1,LNCF),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNCF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,CTC,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'CTC = ',(CTC(L),L=1,LNCF),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNCF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,CTCPR,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'CTCPR = ',(CTCPR(L),L=1,LNCF),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNCF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,CTCPPR,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'CTCPPR = ',(CTCPPR(L),L=1,LNCF),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNAF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,CTAPH,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'CTAPH = ',(CTAPH(L),L=1,LNAF),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNAF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,CTA,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'CTA = ',(CTA(L),L=1,LNAF),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNAF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,CTAPR,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'CTAPR = ',(CTAPR(L),L=1,LNAF),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LNAF
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,CTAPPR,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'CTAPPR = ',(CTAPPR(L),L=1,LNAF),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LMCG*LMCG
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ETH,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ETH = ',(ETH(L),L=1,LMCG*LMCG),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LMCG*LMCG
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ETHP,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ETHP = ',(ETHP(L),L=1,LMCG*LMCG),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LMCG*LMCG
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTR + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,ETHP2,NVAR,MPI_REAL8,
     &  STATUS,IERR)
!      PRINT *,'ETHP2 = ',(ETHP2(L),L=1,LMCG*LMCG),' ID = ',ID
!
!---  Read Pitzer coefficient array
!     (duplicated across processors)  ---
!
      NVAR = LSPL
      OFFSET = IOFFSET + NBYTB
      IOFFSET = IOFFSET + NVAR*NBYTI + 2*NBYTB
      CALL MPI_FILE_READ_AT( IFILE,OFFSET,IDD_PZ,NVAR,MPI_INTEGER,
     &  STATUS,IERR)
!      PRINT *,'IDD_PZ = ',(IDD_PZ(L),L=1,LSPL),' ID = ',ID
!
!---  Close ptcf.bin file  ---
!
      CALL MPI_FILE_CLOSE( IFILE,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of READ_PTZRCOEF group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RESET_SP
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
!
!     Reset reactive-species concentrations with old time-step
!     component-species concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 10 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE REACT
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
      SUB_LOG(ISUB_LOG) = '/RESET_SP'
!
!---  Loop over active nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        DO NSP = 1,NSPR
          SP_C(N,NSP) = SP_CO(N,NSP)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RESET_SP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RMNSP
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChem_X
!
!     Reconstitute mineral species concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 2 May 2006.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE REACT
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
      SUB_LOG(ISUB_LOG) = '/RMNSP'
!
!---  Loop over active nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    Loop over solid species  ---
!
        DO NSPX = 1,NSPS
          NSP = NSPL + NSPX
!
!---      Mineral species  ---
!
          IF( ISP_MN(NSP).EQ.1 ) THEN
            SP_C(N,NSP) = SP_C(N,NSP) + SP_CMN(N,NSPX)
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RMNSP group  ---
!
      RETURN
      END

!----------------------Function--------------------------------------!
!
      FUNCTION TFUNC(A,TX)
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
!     ECKEChemX
!
!     This subroutine calculates higher order electrostatic functions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by S. Yabusaki.
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
      REAL*8 A(8)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TFUNC'
      TFUNC = A(1) + A(2)*TX + A(3)/TX + A(4)*LOG(TX) + 
     &  A(5)/(TX-263.D+0) + A(6)*TX**2 + A(7)/(680.D+0-TX) +
     &  A(8)/(TX-227.D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TFUNC ---
!
      RETURN
      END 

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTCHEM
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     ECKEChemX
!
!     Load old reactive-species concentrations and
!     component-species concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 3 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
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
      SUB_LOG(ISUB_LOG) = '/UPDTCHEM'
!
!---  Loop over active nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        DO NEQ = 1,NEQC+NEQK
          NSL = NEQ + NSOLU
          CO(N,NSL) = C(N,NSL)
        ENDDO
        DO NSP = 1,NSPR
          SP_CO(N,NSP) = SP_C(N,NSP)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDTCHEM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ZLKSRC
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
!     ECKEChemX
!
!     Zero linked sources
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
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
      SUB_LOG(ISUB_LOG) = '/ZLKSRC'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    Zero linked CO2 source  ---
!
        SRCA(1,N) = 0.D+0

      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ZLKSRC group  ---
!
      RETURN
      END



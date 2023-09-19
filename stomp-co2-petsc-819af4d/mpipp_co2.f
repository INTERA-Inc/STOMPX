!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_BIN_CO2
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
!     Write a series of binary input files for execution with STOMPX-CO2
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
      USE CONST
      USE BCVP
      USE BCVS
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IDPX,JDPX,KDPX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IDX,JDX,KDX
      CHARACTER*6 FORM21
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM21 /'(I1)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRITE_BIN_CO2'
!
!---  Allocate temporary memory for mapping of nodes on processors  --
!
      ALLOCATE( IDPX(1:2,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IDPX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( JDPX(1:2,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: JDPX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( KDPX(1:2,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: KDPX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Create a mapping of nodes on processors, including
!     ghost cells  --
!
      DO KP = 1,KP_MPI
        DO JP = 1,JP_MPI
          DO IP = 1,IP_MPI
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            NC = 0
            DO M = 1,2
              IDPX(M,NP) = ID_MPI(M,IP)
              JDPX(M,NP) = JD_MPI(M,JP)
              KDPX(M,NP) = KD_MPI(M,KP)
            ENDDO
            DO K = KD_MPI(1,KP),KD_MPI(2,KP)
              DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                  NC = NC + 1
                  N = ND(I,J,K)
                  ND_MPI(NC,NP) = N
                ENDDO
              ENDDO
            ENDDO
            NC_MPI(NP) = NC
          ENDDO
        ENDDO
      ENDDO
!
!---  Create an array for looping without including ghost cells  ---
!
      ALLOCATE( IDX(1:2,1:LPX_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IDX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( JDX(1:2,1:LPY_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: JDX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( KDX(1:2,1:LPZ_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: KDX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Create a mapping of equation ordering, with nodes
!     counted by processor, then i,j,k indexing, not including
!     ghost cells  --
!
      IL = IFLD/IP_MPI
      JL = JFLD/JP_MPI
      KL = KFLD/KP_MPI
      DO KP = 1,KP_MPI
        IF( KP.EQ.1 ) THEN
          KDX(1,KP) = 1
          KDX(2,KP) = KL
        ELSEIF( KP.EQ.KP_MPI ) THEN
          KDX(1,KP) = (KP-1)*KL + 1
          KDX(2,KP) = KFLD
        ELSE
          KDX(1,KP) = (KP-1)*KL + 1
          KDX(2,KP) = KP*KL
        ENDIF
        DO JP = 1,JP_MPI
          IF( JP.EQ.1 ) THEN
            JDX(1,JP) = 1
            JDX(2,JP) = JL
          ELSEIF( JP.EQ.JP_MPI ) THEN
            JDX(1,JP) = (JP-1)*JL + 1
            JDX(2,JP) = JFLD
          ELSE
            JDX(1,JP) = (JP-1)*JL + 1
            JDX(2,JP) = JP*JL
          ENDIF
          DO IP = 1,IP_MPI
            IF( IP.EQ.1 ) THEN
              IDX(1,IP) = 1
              IDX(2,IP) = IL
            ELSEIF( IP.EQ.IP_MPI ) THEN
              IDX(1,IP) = (IP-1)*IL + 1
              IDX(2,IP) = IFLD
            ELSE
              IDX(1,IP) = (IP-1)*IL + 1
              IDX(2,IP) = IP*IL
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Create a mapping of nodes to processors, not including
!     ghost cells  ---
!
      DO KP = 1,KP_MPI
        DO JP = 1,JP_MPI
          DO IP = 1,IP_MPI
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            DO K = KDX(1,KP),KDX(2,KP)
              DO J = JDX(1,JP),JDX(2,JP)
                DO I = IDX(1,IP),IDX(2,IP)
                  N = ND(I,J,K)
                  NTP_MPI(N) = NP
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Write solution control and output data to solu.bin  ---
!
      CALL WRITE_SOLU_CO2( IDX,JDX,KDX )
!
!---  Write grid data to a series of binary grid files  ---
!
      CALL WRITE_GRID_CO2( IDPX,JDPX,KDPX )
!
!---   Write property data to mech.bin, hydr.bin satu.bin, perm.bin,
!     comp.bin, and tabl.bin
!
      CALL WRITE_PROP_CO2
!
!---  Write state condition data to a series of binary files  ---
!
      CALL WRITE_STATE_CO2
!
!---  Write boundary condition data to boco.bin  ---
!
      CALL WRITE_BOCO_CO2
!
!---  Write source data to sorc.bin  ---
!
      CALL WRITE_SORC_CO2
!
!---  Write coupled-well data to well.bin  ---
!
      CALL WRITE_COUP_WELL_CO2( IDX,JDX,KDX )
!
!---  Write solute/species transport data to tpor.bin  ---
!
      CALL WRITE_TPORT_CO2
!
!---  Write geomechanics property data to gmec.bin  ---
!
      IF( ISLC(50).NE.0 ) CALL WRITE_GMEC_CO2
!
!---  Write geomechanics grid data to gmgd.bin  ---
!
      IF( ISLC(50).NE.0 ) CALL WRITE_GMGD_CO2
!
!---  Write geomechanics boundary condition data to gmbc.bin  ---
!
      IF( ISLC(50).NE.0 ) CALL WRITE_GMBC_CO2( IDX,JDX,KDX )

!
!---  Reactive transport  ---
!
      IF( ISLC(40).NE.0 ) THEN
!
!---    Write ECKEChem inputs to 1 .bin file  ---
!
        CALL WRITE_REACT   
!
!---    Write ptcf.bin  ---
!
        IF( IACTV.EQ.2 ) CALL WRITE_PTZR
      ENDIF

!
!---  Write co2_prop.bin file  ---
!
      CALL WRITE_EOS_CO2
!
!---  Deallocate temporary memory for field node equation pointers  ---
!
      DEALLOCATE( IDX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IDX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( JDX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: JDX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( KDX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: KDX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Deallocate temporary memory for mapping of nodes on processors  --
!
      DEALLOCATE( IDPX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IDPX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( JDPX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: JDPX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( KDPX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: KDPX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_BIN_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_BOCO_CO2
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
!     Write boundary condition data to bcsc.bin
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
      USE CONST
      USE BCVP
      USE BCVS
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: PLBX,PGBX,YLSBX
      REAL*8, DIMENSION(:), ALLOCATABLE :: XPBCX,YPBCX,ZPBCX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IBCTX,IBCSPX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCCX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCNX,IBCDX,IBCMX,IBCINX
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBCX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRITE_BOCO_CO2'
!
!---  Open a new boco.bin file  --
!
      OPEN( UNIT=20,FILE='boco.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='boco.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Determine the number of boundary condition surfaces on each
!     processor, not including ghost cells  --
!
      ALLOCATE( NBCX(1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NBCX'
        CALL WRMSGS( INDX )
      ENDIF
      DO NP = 1,NP_MPI
        NBCX(NP) = 0
      ENDDO
      DO NB = 1,NBC
        N = IBCN(NB)
        NP = NTP_MPI(N)
        NBCX(NP) = NBCX(NP) + 1
      ENDDO
      NBCSX = 0
      DO NP = 1,NP_MPI
        NBCSX = NBCSX + NBCX(NP)
      ENDDO
      WRITE(20) (NBCX(NP),NP=1,NP_MPI)
!
!---  Write boundary condition input variables to bcsc.bin  --
!
      WRITE(20) (((BC(L,M,N),L=1,LBCV),M=1,LBTM),N=1,LBCIN)
!
!---  Allocate processor local temporary memory for boundary
!     condition indices  --
!
      ALLOCATE( IBCNX(1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCNX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IBCDX(1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCDX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IBCMX(1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCMX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IBCINX(1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCINX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IBCCX(1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCCX'
        CALL WRMSGS( INDX )
      ENDIF
      LX = LUK+LSOLU*LC
      ALLOCATE( IBCTX(1:LX,1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCTX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IBCSPX(1:(LSPBC+1),1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCSPX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( XPBCX(1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: XPBCX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( YPBCX(1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: YPBCX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( ZPBCX(1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ZPBCX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Create a processor-local map of boundary condition field nodes,
!     directions, number of time points, input links, cycling options,
!     and types, not including ghost cells  ---
!
      NC = 0
      DO KP = 1,KP_MPI
        DO JP = 1,JP_MPI
          DO IP = 1,IP_MPI
            KFLDX = KD_MPI(2,KP)-KD_MPI(1,KP)+1
            JFLDX = JD_MPI(2,JP)-JD_MPI(1,JP)+1
            IFLDX = ID_MPI(2,IP)-ID_MPI(1,IP)+1
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            DO NB = 1,NBC
              N = IBCN(NB)
              IF( NP.EQ.NTP_MPI(N) ) THEN
                NC = NC + 1
                I = ID(N)
                J = JD(N)
                K = KD(N)
                IBCNX(NC) = (K-KD_MPI(1,KP))*IFLDX*JFLDX +
     &            (J-JD_MPI(1,JP))*IFLDX + (I-ID_MPI(1,IP)) + 1
                IBCDX(NC) = IBCD(NB)
                IBCMX(NC) = IBCM(NB)
                IBCINX(NC) = IBCIN(NB)
                IBCCX(NC) = IBCC(NB)
                XPBCX(NC) = XPBC(NB)
                YPBCX(NC) = YPBC(NB)
                ZPBCX(NC) = ZPBC(NB)
                DO M = 1,LX
                  IBCTX(M,NC) = IBCT(M,NB)
                ENDDO
                DO M = 1,(LSPBC+1)
                  IBCSPX(M,NC) = IBCSP(M,NB)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      WRITE(20) (IBCNX(M),M=1,NBCSX)
      WRITE(20) (IBCDX(M),M=1,NBCSX)
      WRITE(20) (IBCMX(M),M=1,NBCSX)
      WRITE(20) (IBCINX(M),M=1,NBCSX)
      WRITE(20) (IBCCX(M),M=1,NBCSX)
      WRITE(20) ((IBCTX(L,M),L=1,LX),M=1,NBCSX)
      WRITE(20) ((IBCSPX(L,M),L=1,(LSPBC+1)),M=1,NBCSX)
      WRITE(20) (XPBCX(M),M=1,NBCSX)
      WRITE(20) (YPBCX(M),M=1,NBCSX)
      WRITE(20) (ZPBCX(M),M=1,NBCSX)
!
!---  Deallocate processor local temporary memory for boundary
!     condition indices  --
!
      DEALLOCATE( IBCNX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCNX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IBCDX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCDX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IBCMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCMX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IBCINX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCINX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IBCCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCCX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IBCTX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCTX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IBCSPX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCSPX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( XPBCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: XPBCX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( YPBCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: YPBCX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( ZPBCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: ZPBCX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Allocate processor local temporary memory for initial
!     condition boundary conditions  --
!
      ALLOCATE( PLBX(1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: PLBX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( PGBX(1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: PGBX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( YLSBX(1:NBCSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: YLSBX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Create a processor-local map of boundary condition field nodes,
!     directions, number of time points, input links, cycling options,
!     and types, not including ghost cells  ---
!
      NC = 0
      DO KP = 1,KP_MPI
        DO JP = 1,JP_MPI
          DO IP = 1,IP_MPI
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            DO NB = 1,NBC
              N = IBCN(NB)
              IF( NP.EQ.NTP_MPI(N) ) THEN
                NC = NC + 1
                PLBX(NC) = PLB(1,NB)
                PGBX(NC) = PGB(1,NB)
                YLSBX(NC) = YLSB(1,NB)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      WRITE(20) (PLBX(M),M=1,NBCSX)
      WRITE(20) (PGBX(M),M=1,NBCSX)
      WRITE(20) (YLSBX(M),M=1,NBCSX)
!
!---  Deallocate processor local temporary memory for initial
!     condition boundary conditions  --
!
      DEALLOCATE( PLBX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: PLBX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( PGBX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: PGBX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( YLSBX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: YLSBX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Deallocate temporary boundary condition memory  --
!
      DEALLOCATE( NBCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NBCX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Close the boco.bin file  ---
!
      CLOSE(20)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_BOCO_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_COUP_WELL_CO2( IDX,JDX,KDX )
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
!     Write state condition data to wltp.bin
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
      USE CONST
      USE BCVP
      USE BCVS
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER, DIMENSION(1:2,1:LPX_MPI) :: IDX
      INTEGER, DIMENSION(1:2,1:LPY_MPI) :: JDX
      INTEGER, DIMENSION(1:2,1:LPZ_MPI) :: KDX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IXP_CWX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IXPX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRITE_COUP_WELL_CO2'
!
!---  Open a new well.bin file  --
!
      OPEN( UNIT=20,FILE='well.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='well.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Allocate temporary memory  ---
!
      ALLOCATE( IXPX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IXPX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Write number of coupled wells  ---
!
      WRITE(20) N_CW
!
!---  Write coupled well parameters, indices and variables
!     if coupled wells are modeled  ---
!
      IF( N_CW.GT.0 ) THEN
!
!---    Write name of coupled wells  ---
!
        WRITE(20) (WNM_CW(N),N=1,LN_CW)
!
!---    Write type index of coupled wells  ---
!
        WRITE(20) (IT_CW(N),N=1,LN_CW)
!
!---    Write injection well state index of coupled wells  ---
!
        WRITE(20) ((ITS_CW(M,N),M=1,LWTP_CW),N=1,LN_CW)
!
!---    Write x-, y-, and z-direction coupled well fraction factors  ---
!
        WRITE(20) ((FF_CW(M,N),M=1,3),N=1,LN_CW)
!
!---    Write coupled well indices  ---
!
        WRITE(20) ((ID_CW(M,N),M=1,10),N=1,LN_CW)
!
!---    Write coupled well total mass limit  ---
!
        WRITE(20) (TML_CW(N),N=1,LN_CW)
!
!---    Write cyclic index  ---
!
        WRITE(20) (ICC_CW(N),N=1,LN_CW)
!
!---    Write number of time points  ---
!
        WRITE(20) (IM_CW(N),N=1,LN_CW)
!
!---    Write number of time points per time period  ---
!
        WRITE(20) ((IMP_CW(M,N),M=1,LWTP_CW),N=1,LN_CW)
!
!---    Write coupled well reference nodes  ---
!
        WRITE(20) (IREF_CW(N),N=1,LVREF)
!
!---    Allocate temporary memory to create an array of well node to 
!       field node equation pointers  ---
!
        ALLOCATE( IXP_CWX(1:LWF_CW),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IXP_CWX'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Coupled equations using I,J,K ordering on parallel 
!       processors, not including ghost cells  ---
!
        NC = 0
        DO KP = 1,KP_MPI
          DO JP = 1,JP_MPI
            DO IP = 1,IP_MPI
              DO K = KDX(1,KP),KDX(2,KP)
                DO J = JDX(1,JP),JDX(2,JP)
                  DO I = IDX(1,IP),IDX(2,IP)
                    N = ND(I,J,K)
                    IF( IXP(N).EQ.0 ) CYCLE
                    NC = NC + 1
                    IXPX(N) = NC
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
!---    Loop over coupled wells  ---
!
        DO NCW = 1,N_CW
!
!---      Loop over field nodes with coupled-well nodes ---
!
          DO NWF = ID_CW(5,NCW),ID_CW(6,NCW)
            N = IWF_CW(NWF)
!
!---        Assign the field-node equation number, according to 
!           coupled equations, using I,J,K ordering on parallel
!           processors ---
!
            IXP_CWX(NWF) = IXPX(N)
          ENDDO
        ENDDO
!
!---    Write field to well node pointers  ---
!
        WRITE(20) (IWF_CW(N),N=1,LWF_CW)
        WRITE(20) (IXP_CWX(N),N=1,LWF_CW)
        WRITE(20) (IWN_CW(N),N=1,LWN_CW)
        WRITE(20) (IWP_CW(N),N=1,LWN_CW)
!
!---    Write coupled well arrays  ---
!
        WRITE(20) (IS_CW(M),M=1,LWI_CW)
        WRITE(20) ((XTP_CW(L,M),L=1,2),M=1,LWI_CW)
        WRITE(20) ((YTP_CW(L,M),L=1,2),M=1,LWI_CW)
        WRITE(20) ((ZTP_CW(L,M),L=1,2),M=1,LWI_CW)
        WRITE(20) ((PAR_CW(L,M),L=1,5),M=1,LWI_CW)
        WRITE(20) (INV_CW(M),M=1,LWN_CW)
!
!---    Write time dependent driving parameters for the 
!       coupled well  ---
!
        WRITE(20) (((VAR_CW(L,M,N),L=1,8),M=1,LWT_CW),N=1,LN_CW)
!
!---    Write time dependent solute driving parameters for the coupled
!       well  ---
!
        WRITE(20) (((VARC_CW(L,M,N),L=1,LSOLU_CW),M=1,LWT_CW),N=1,LN_CW)
!
!---    Write time dependent species driving parameters for the coupled
!       well  ---
!
        WRITE(20) (((VARSP_CW(L,M,N),L=1,LSPC_CW),M=1,LWT_CW),N=1,LN_CW)
!
!---    Write x-direction well projections  ---
!
        WRITE(20) (PLX_CW(N),N=1,LWN_CW)
!
!---    Write y-direction well projections  ---
!
        WRITE(20) (PLY_CW(N),N=1,LWN_CW)
!
!---    Write z-direction well projections  ---
!
        WRITE(20) (PLZ_CW(N),N=1,LWN_CW)
!
!---    Write x-direction well node points  ---
!
        WRITE(20) ((XP_CW(M,N),M=1,2),N=1,LWN_CW)
!
!---    Write y-direction well node points  ---
!
        WRITE(20) ((YP_CW(M,N),M=1,2),N=1,LWN_CW)
!
!---    Write z-direction well node points  ---
!
        WRITE(20) ((ZP_CW(M,N),M=1,2),N=1,LWN_CW)
!
!---    Write active solutes in coupled well  ---
!
        WRITE(20) ((ISOLU_CW(M,N),M=1,LSOLU_CW),N=1,LN_CW)
!
!---    Write active species in coupled well  ---
!
        WRITE(20) ((ISPC_CW(M,N),M=1,LSPC_CW),N=1,LN_CW)
!
!---    Write active component species in coupled well  ---
!
        WRITE(20) ((ISOLC_CW(M,N),M=1,LEQC),N=1,LN_CW)
!
!---    Write active kinetic species in coupled well  ---
!
        WRITE(20) ((ISOLK_CW(M,N),M=1,LEQK),N=1,LN_CW)
!
!---    Write number of species in coupled well  ---
!
        WRITE(20) (NSP_CW(N),N=1,LN_CW)
!
!---    Deallocate temporary coupled-well index  ---
!
        DEALLOCATE( IXP_CWX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: IXP_CWX'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Deallocate temporary memory  ---
!
      DEALLOCATE( IXPX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IXPX'
        CALL WRMSGS( INDX )
      ENDIF
!!
!!---  Write number of leaky wells  ---
!!
!      WRITE(20) N_LW
!!
!!---  Write name of leaky wells  ---
!!
!      WRITE(20) (WNM_LW(N),N=1,N_LW)
!!
!!---  Write type index of leaky wells  ---
!!
!      WRITE(20) (IT_LW(N),N=1,N_LW)
!!
!!---  Write x-, y-, and z-direction leaky well fraction factors  ---
!!
!      WRITE(20) ((FF_LW(M,N),M=1,3),N=1,N_LW)
!!
!!---  Write leaky well indices  ---
!!
!      WRITE(20) ((ID_LW(M,N),M=1,5),N=1,N_LW)
!!
!!---  Loop over leaky wells  ---
!!
!      DO N = 1,N_LW
!!
!!---    Write type index for leaky well interval  ---
!!
!        WRITE(20) (IS_LW(M),M=ID_LW(1,N),ID_LW(2,N))
!!
!!---    Write rock/soil/grout/fill type for leaky well interval  ---
!!
!        WRITE(20) (IZ_LW(M),M=ID_LW(1,N),ID_LW(2,N))
!!
!!---    Write x-transition points for leaky well interval  ---
!!
!        WRITE(20) ((XTP_LW(L,M),L=1,2),M=ID_LW(1,N),ID_LW(2,N))
!!
!!---    Write y-transition points for leaky well interval  ---
!!
!        WRITE(20) ((YTP_LW(L,M),L=1,2),M=ID_LW(1,N),ID_LW(2,N))
!!
!!---    Write z-transition points for leaky well interval  ---
!!
!        WRITE(20) ((ZTP_LW(L,M),L=1,2),M=ID_LW(1,N),ID_LW(2,N))
!!
!!---    Write parameters for leaky well interval  ---
!!
!        WRITE(20) ((PAR_LW(L,M),L=1,5),M=ID_LW(1,N),ID_LW(2,N))
!      ENDDO
!
!---  Close the well.bin file  ---
!
      CLOSE(20)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_COUP_WELL_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_EOS_CO2
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
!     Write co2_prop.bin file
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
      USE CONST
      USE BCVP
      USE BCVS
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/WRITE_EOS_CO2'
!
!---  Open a new co2_prop.bin file  --
!
      OPEN( UNIT=20,FILE='co2_prop.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='co2_prop.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
      INCG = 1
!
!---  Write number of pressure points  ---
!
      WRITE(20) IP_TA(INCG)
!
!---  Write pressure points  ---
!
      WRITE(20) (P_TA(I,INCG),I=1,IP_TA(INCG))
!
!---  Write number of temperature points for each pressure point  ---
!
      WRITE(20) (IT_TA(I,INCG),I=1,IP_TA(INCG))
!
!---  Loop over the number of pressure points, writing temperatures  ---
!
      DO IPX = 1,IP_TA(INCG)
        WRITE(20) (T_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
      ENDDO
!
!---  Loop over the number of pressure points, writing density  ---
!
      DO IPX = 1,IP_TA(INCG)
        WRITE(20) (RHO_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
      ENDDO
!
!---  Loop over the number of pressure points, writing enthalpy  ---
!
      DO IPX = 1,IP_TA(INCG)
        WRITE(20) (H_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
      ENDDO
!
!---  Loop over the number of pressure points, writing internal 
!     energy  ---
!
      DO IPX = 1,IP_TA(INCG)
        WRITE(20) (U_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
      ENDDO
!
!---  Loop over the number of pressure points, writing fugacity  ---
!
      DO IPX = 1,IP_TA(INCG)
        WRITE(20) (FUG_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
      ENDDO
!
!---  Loop over the number of pressure points, writing entropy  ---
!
      DO IPX = 1,IP_TA(INCG)
        WRITE(20) (S_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
      ENDDO
!
!---  Write index of vapor-liquid point  ---
!
      WRITE(20) (IV_TA(I,INCG),I=1,IP_TA(INCG))
!
!---  Write number of vapor-liquid data points  ---
!
      WRITE(20) I_LV(INCG)
!
!---  Loop over the number of vapor-liquid data points, writing 
!     temperature, pressure, liquid density, liquid enthalpy,
!     liquid internal energy, liquid entropy, vapor density,
!     vapor enthalpy, vapor internal energy, internal entropy,
!     vapor fugacity  ---
!
      DO IPX = 1,IP_TA(INCG)
        WRITE(20) T_LV(IPX,INCG),P_LV(IPX,INCG),RHOL_LV(IPX,INCG),
     &    HL_LV(IPX,INCG),UL_LV(IPX,INCG),SL_LV(IPX,INCG),
     &    RHOV_LV(IPX,INCG),HV_LV(IPX,INCG),UV_LV(IPX,INCG),
     &    SV_LV(IPX,INCG),FUG_LV(IPX,INCG)
      ENDDO
!
!---  Close the co2_prop.bin file  ---
!
      CLOSE(20)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_EOS_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_GMBC_CO2( IDX,JDX,KDX )
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
!     Write geomechanical boundary condition data to gmbc.bin
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 September 2022 (orange shirt day).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
      USE CONST
      USE BCVP
      USE BCVS
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBC_GMX
      INTEGER, DIMENSION(:), ALLOCATABLE :: NTL_MPI,NTR_MPI
      INTEGER, DIMENSION(1:2,1:LPX_MPI) :: IDX
      INTEGER, DIMENSION(1:2,1:LPY_MPI) :: JDX
      INTEGER, DIMENSION(1:2,1:LPZ_MPI) :: KDX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IBCT_GMX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCC_GMX,IBCN_GMX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCD_GMX,IBCM_GMX,IBCIN_GMX
      INTEGER, DIMENSION(8) :: IPX,IP1X,IP2X,JPX,JP1X,JP2X,KPX,KP1X,KP2X
!
!----------------------Data Statements---------------------------------!
!
      DATA IP1X / 1,2,3,4,5,6,7,8 /
      DATA IP2X / 0,2,0,4,0,6,0,8 /
      DATA JP1X / 1,2,3,4,5,6,7,8 /
      DATA JP2X / 0,0,3,4,0,0,7,8 /
      DATA KP1X / 1,2,3,4,5,6,7,8 /
      DATA KP2X / 0,0,0,0,5,6,7,8 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRITE_GMBC_CO2'
!
!---  Allocate memory for mapping of finite element nodes to 
!     processors (NTR_MPI) and global finite element nodes to local
!     finite element nodes (NTL_MPI)  --
!
      ALLOCATE( NTR_MPI(1:LFEN),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NTR_MPI'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( NTL_MPI(1:LFEN),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NTL_MPI'
        CALL WRMSGS( INDX )
      ENDIF
      DO L = 1,LFEN
        NTR_MPI(L) = 0
        NTL_MPI(L) = 0
      ENDDO
!
!---  Create a mapping of finite element nodes to processors, 
!     not including ghost cells  ---
!
      DO KP = 1,KP_MPI
        DO JP = 1,JP_MPI
          DO IP = 1,IP_MPI
            IF( IP.EQ.1 ) THEN
              DO M = 1,8
                IPX(M) = IP1X(M)
              ENDDO
            ELSE
              DO M = 1,8
                IPX(M) = IP2X(M)
              ENDDO
            ENDIF
            IF( JP.EQ.1 ) THEN
              DO M = 1,8
                JPX(M) = JP1X(M)
              ENDDO
            ELSE
              DO M = 1,8
                JPX(M) = JP2X(M)
              ENDDO
            ENDIF
            IF( KP.EQ.1 ) THEN
              DO M = 1,8
                KPX(M) = KP1X(M)
              ENDDO
            ELSE
              DO M = 1,8
                KPX(M) = KP2X(M)
              ENDDO
            ENDIF
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            DO K = KDX(1,KP),KDX(2,KP)
              DO J = JDX(1,JP),JDX(2,JP)
                DO I = IDX(1,IP),IDX(2,IP)
                  N = ND(I,J,K)
                  DO M = 1,8
                    MX = MAX( IPX(M),JPX(M),KPX(M) )
                    IF( MX.EQ.0 ) CYCLE
                    NE = ND_GM(M,N)
                    NTR_MPI(NE) = NP
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Open a new gmbc.bin file  --
!
      OPEN( UNIT=20,FILE='gmbc.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='gmbc.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Determine the number of geomechanical boundary conditions
!     on each processor, not including ghost cells  --
!
      ALLOCATE( NBC_GMX(1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NBC_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      DO NP = 1,NP_MPI
        NBC_GMX(NP) = 0
      ENDDO
      DO NB = 1,NBC_GM
        N = IBCN_GM(NB)
        IFENX = 1
        DO M = 1,3
          IF( IBCT_GM(M,NB).EQ.1 .OR. IBCT_GM(M,NB).EQ.2 ) IFENX = 0
        ENDDO
!
!---    Boundary node number refers to a finite-element node number  --
!
        IF( IFENX.EQ.1 ) THEN
          NP = NTR_MPI(N)
!
!---    Boundary node number refers to a field node number  --
!
        ELSE
          NP = NTP_MPI(N)
        ENDIF
        NBC_GMX(NP) = NBC_GMX(NP) + 1
      ENDDO
      NBCS_GMX = 0
      DO NP = 1,NP_MPI
        NBCS_GMX = NBCS_GMX + NBC_GMX(NP)
      ENDDO
!
!---  Write number of geomechanical boundary conditions
!     on each processor, not including ghost cells to gmbc.bin  --
!
      WRITE(20) (NBC_GMX(NP),NP=1,NP_MPI)
!
!---  Write geomechanical boundary condition input variables 
!     to gmbc.bin  --
!
      WRITE(20) (((BC_GM(L,M,N),L=1,LBCV_GM),M=1,LBTM_GM),N=1,LBCIN_GM)
!
!---  Write reference geomechanical boundary condition index
!     to gmbc.bin  --
!
      WRITE(20) (IBCR_GM(N),N=1,LBCIN_GM)
!
!---  Allocate processor local temporary memory for boundary
!     condition indices  --
!
      ALLOCATE( IBCN_GMX(1:NBCS_GMX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCN_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IBCD_GMX(1:NBCS_GMX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCD_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IBCM_GMX(1:NBCS_GMX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCM_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IBCIN_GMX(1:NBCS_GMX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCIN_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IBCC_GMX(1:NBCS_GMX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCC_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IBCT_GMX(1:3,1:NBCS_GMX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IBCT_GMX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Create a processor-local map of geomechanical boundary condition,
!     nodes, directions, number of time points, input links, 
!     cycling options, and types, not including ghost cells  ---
!
      NC = 0
      DO KP = 1,KP_MPI
        DO JP = 1,JP_MPI
          DO IP = 1,IP_MPI
            KFLDX = KD_MPI(2,KP)-KD_MPI(1,KP)+1
            JFLDX = JD_MPI(2,JP)-JD_MPI(1,JP)+1
            IFLDX = ID_MPI(2,IP)-ID_MPI(1,IP)+1
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            NCX = 0
!
!---       For processor NP create a map of global to local
!          finite element node numbers  ---
!
            DO K = KD_MPI_GM(1,KP),KD_MPI_GM(2,KP)
              DO J = JD_MPI_GM(1,JP),JD_MPI_GM(2,JP)
                DO I = ID_MPI_GM(1,IP),ID_MPI_GM(2,IP)
                  NCX = NCX + 1
                  NFEN = (K-1)*(IFLD+1)*(JFLD+1) + (J-1)*(IFLD+1) + I
                  NTL_MPI(NFEN) = NCX
                ENDDO
              ENDDO
            ENDDO
            DO NB = 1,NBC_GM
              IFENX = 1
              DO M = 1,3
                IF( IBCT_GM(M,NB).EQ.1 .OR. 
     &            IBCT_GM(M,NB).EQ.2 ) IFENX = 0
              ENDDO
              N = IBCN_GM(NB)
!
!---          Boundary node number refers to a finite-element 
!             node number  --
!
              IF( IFENX.EQ.1 ) THEN
                IF( NP.EQ.NTR_MPI(N) ) THEN
                  NC = NC + 1
                  MB = IBCIN_GM(NB)
                  IBCN_GMX(NC) = NTL_MPI(N)
                  IBCD_GMX(NC) = IBCD_GM(NB)
                  IBCM_GMX(NC) = IBCM_GM(MB)
                  IBCIN_GMX(NC) = IBCIN_GM(NB)
                  IBCC_GMX(NC) = IBCC_GM(NB)
                  DO M = 1,3
                    IBCT_GMX(M,NC) = IBCT_GM(M,NB)
                  ENDDO
                ENDIF
!
!---          Boundary node number refers to a field node number  --
!
              ELSE
                IF( NP.EQ.NTP_MPI(N) ) THEN
                  NC = NC + 1
                  I = ID(N)
                  J = JD(N)
                  K = KD(N)
                  MB = IBCIN_GM(NB)
                  IBCN_GMX(NC) = (K-KD_MPI(1,KP))*IFLDX*JFLDX +
     &              (J-JD_MPI(1,JP))*IFLDX + (I-ID_MPI(1,IP)) + 1
                  IBCD_GMX(NC) = IBCD_GM(NB)
                  IBCM_GMX(NC) = IBCM_GM(MB)
                  IBCIN_GMX(NC) = IBCIN_GM(NB)
                  IBCC_GMX(NC) = IBCC_GM(NB)
                  DO M = 1,3
                    IBCT_GMX(M,NC) = IBCT_GM(M,NB)
                  ENDDO
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      WRITE(20) (IBCN_GMX(M),M=1,NBCS_GMX)
      WRITE(20) (IBCD_GMX(M),M=1,NBCS_GMX)
      WRITE(20) (IBCM_GMX(M),M=1,NBCS_GMX)
      WRITE(20) (IBCIN_GMX(M),M=1,NBCS_GMX)
      WRITE(20) (IBCC_GMX(M),M=1,NBCS_GMX)
      WRITE(20) ((IBCT_GMX(L,M),L=1,3),M=1,NBCS_GMX)
!
!---  Deallocate processor local temporary memory for boundary
!     condition indices  --
!
      DEALLOCATE( IBCN_GMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCN_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IBCD_GMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCD_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IBCM_GMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCM_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IBCIN_GMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCIN_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IBCC_GMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCC_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IBCT_GMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IBCT_GMX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Deallocate temporary boundary condition memory  --
!
      DEALLOCATE( NTR_MPI,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NTR_MPI'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( NTL_MPI,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NTL_MPI'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( NBC_GMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NBC_GMX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Close the gmbc.bin file  ---
!
      CLOSE(20)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_GMBC_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_GMEC_CO2
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
!     Write geomechanics property data to gmec.bin
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 September 2022 (orange shirt day).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
      USE CONST
      USE BCVP
      USE BCVS
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PROP_GMX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IPROP_GMX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IHCM_GMX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRITE_GMEC_CO2'
!
!---  Allocate temporary memory for geomechanical property arrays  --
!
      ALLOCATE( PROP_GMX(1:7,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: PROP_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IHCM_GMX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IHCM_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IPROP_GMX(1:3,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IPROP_GMX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Fill temporary geomechanical property arrays  --
!
      DO N = 1,NFLD
        IF( IZ(N).LE.0 ) THEN
          DO M = 1,7
            PROP_GMX(M,N) = 0.D+0
          ENDDO
          DO M = 1,3
            IPROP_GMX(M,N) = 0
          ENDDO
          IHCM_GMX(N) = 0
        ELSE
          DO M = 1,7
            PROP_GMX(M,N) = PROP_GM(M,IZ(N))
          ENDDO
          DO M = 1,3
            IPROP_GMX(M,N) = IPROP_GM(M,IZ(N))
          ENDDO
          IHCM_GMX(N) = IHCM_GM(IZ(N))
        ENDIF
      ENDDO
!
!---  Open a new gmec.bin file  --
!
      OPEN( UNIT=20,FILE='gmec.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='gmec.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write geomechanical residual limit array to gmec.bin  --
!
      WRITE(20) (RSDM_GM(L),L=1,LEPD)
!
!---  Write geomechanical property arrays to gmec.bin  --
!
      WRITE(20) (((PROP_GMX(M,ND_MPI(K,L)),M=1,7),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      WRITE(20) (((IPROP_GMX(M,ND_MPI(K,L)),M=1,3),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      WRITE(20) ((IHCM_GMX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Close the gmec.bin file  ---
!
      CLOSE(20)
!
!---  Deallocate temporary memory for geomechanical property arrays  --
!
      DEALLOCATE( PROP_GMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: PROP_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IHCM_GMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IHCM_GMX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IPROP_GMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IPROP_GMX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_GMEC_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_GMGD_CO2
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
!     Write geomechanics grid data to gmgd.bin
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 04 October 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
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
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IDPX_GM,JDPX_GM,KDPX_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: NF_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: NGHNX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NCSGNX,NPSGNX
      INTEGER, DIMENSION(:), ALLOCATABLE :: NLSGNX,NLRGNX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRITE_GMGD_CO2'
!
!---  Allocate memory for global finite element nodes  --
!
      ALLOCATE( NF_GM(1:LFEN),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NF_GM'
        CALL WRMSGS( INDX )
      ENDIF
      DO NFEN = 1,NFEN_GM
        NF_GM(NFEN) = NFEN
      ENDDO
!
!---  Allocate temporary memory for mapping of nodes on processors  --
!
      ALLOCATE( IDPX_GM(1:2,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IDPX_GM'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( JDPX_GM(1:2,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: JDPX_GM'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( KDPX_GM(1:2,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: KDPX_GM'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Set processor indice limits for finite-element nodes  ---
!
      IL = (IFLD+1)/IP_MPI
      JL = (JFLD+1)/JP_MPI
      KL = (KFLD+1)/KP_MPI
      DO IP = 1,IP_MPI
        IF( IP.EQ.1 ) THEN
          ID_MPI_GM(1,IP) = 1
          ID_MPI_GM(2,IP) = IL + 1
        ELSEIF( IP.EQ.IP_MPI ) THEN
          ID_MPI_GM(1,IP) = (IP-1)*IL
          ID_MPI_GM(2,IP) = (IFLD+1)
        ELSE
          ID_MPI_GM(1,IP) = (IP-1)*IL
          ID_MPI_GM(2,IP) = IP*IL + 1
        ENDIF
      ENDDO
      DO JP = 1,JP_MPI
        IF( JP.EQ.1 ) THEN
          JD_MPI_GM(1,JP) = 1
          JD_MPI_GM(2,JP) = JL + 1
        ELSEIF( JP.EQ.JP_MPI ) THEN
          JD_MPI_GM(1,JP) = (JP-1)*JL
          JD_MPI_GM(2,JP) = (JFLD+1)
        ELSE
          JD_MPI_GM(1,JP) = (JP-1)*JL
          JD_MPI_GM(2,JP) = JP*JL + 1
        ENDIF
      ENDDO
      DO KP = 1,KP_MPI
        IF( KP.EQ.1 ) THEN
          KD_MPI_GM(1,KP) = 1
          KD_MPI_GM(2,KP) = KL + 1
        ELSEIF( KP.EQ.KP_MPI ) THEN
          KD_MPI_GM(1,KP) = (KP-1)*KL
          KD_MPI_GM(2,KP) = (KFLD+1)
        ELSE
          KD_MPI_GM(1,KP) = (KP-1)*KL
          KD_MPI_GM(2,KP) = KP*KL + 1
        ENDIF
      ENDDO
      IF( IP_MPI.EQ.1 ) THEN
        ID_MPI_GM(1,IP_MPI) = 1
        ID_MPI_GM(2,IP_MPI) = (IFLD+1)
      ENDIF
      IF( JP_MPI.EQ.1 ) THEN
        JD_MPI_GM(1,JP_MPI) = 1
        JD_MPI_GM(2,JP_MPI) = (JFLD+1)
      ENDIF
      IF( KP_MPI.EQ.1 ) THEN
        KD_MPI_GM(1,KP_MPI) = 1
        KD_MPI_GM(2,KP_MPI) = (KFLD+1)
      ENDIF
!
!---  Create a mapping of finite-element nodes on processors, including
!     FE ghost nodes  --
!
      DO KP = 1,KP_MPI
        DO JP = 1,JP_MPI
          DO IP = 1,IP_MPI
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            DO M = 1,2
              IDPX_GM(M,NP) = ID_MPI_GM(M,IP)
              JDPX_GM(M,NP) = JD_MPI_GM(M,JP)
              KDPX_GM(M,NP) = KD_MPI_GM(M,KP)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Open a new gmgd.bin file  --
!
      OPEN( UNIT=20,FILE='gmgd.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='gmgd.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write processor indice limits for finite-element nodes  ---
!
      WRITE(20) ((IDPX_GM(M,NP),M=1,2),NP=1,NP_MPI)
      WRITE(20) ((JDPX_GM(M,NP),M=1,2),NP=1,NP_MPI)
      WRITE(20) ((KDPX_GM(M,NP),M=1,2),NP=1,NP_MPI)
!
!---  Write global finite-element node pointer data to gmgd.bin  ---
!
      WRITE(20) (((ND_GM(M,ND_MPI(K,L)),M=1,8),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
!
!---  Deallocate memory for field node mapping on processors  ---
!
      DEALLOCATE( ND_MPI,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: ND_MPI'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Allocate memory for finite-element node mapping on processors  ---
!
      MFX_MPI = LFX_MPI+3
      IF( LFX_MPI.EQ.1 ) MFX_MPI = 2
      MFY_MPI = LFY_MPI+3
      IF( LFY_MPI.EQ.1 ) MFY_MPI = 2
      MFZ_MPI = LFZ_MPI+3
      IF( LFZ_MPI.EQ.1 ) MFZ_MPI = 2
      ALLOCATE( ND_MPI(MFX_MPI*MFY_MPI*MFZ_MPI,LP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ND_MPI'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Create a mapping of finite-element nodes on processors  ---
!
      DO KP = 1,KP_MPI
        DO JP = 1,JP_MPI
          DO IP = 1,IP_MPI
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            NC = 0
            DO K = KD_MPI_GM(1,KP),KD_MPI_GM(2,KP)
              DO J = JD_MPI_GM(1,JP),JD_MPI_GM(2,JP)
                DO I = ID_MPI_GM(1,IP),ID_MPI_GM(2,IP)
                  NC = NC + 1
                  NFEN = (K-1)*(IFLD+1)*(JFLD+1) + (J-1)*(IFLD+1) + I
                  ND_MPI(NC,NP) = NFEN
                ENDDO
              ENDDO
            ENDDO
            NC_MPI(NP) = NC
          ENDDO
        ENDDO
      ENDDO
!
!---  Write number of FE nodes on processors, including FE ghost
!     nodes to gmgd.bin  ---
!
      WRITE(20) (NC_MPI(NP),NP=1,NP_MPI)
!
!---  Write global field node pointer grid data to gmgd.bin  ---
!
      WRITE(20) (((NE_GM(M,ND_MPI(K,L)),M=1,8),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
!
!---  Write global finite-element node data to gmgd.bin  ---
!
      WRITE(20) ((NF_GM(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Write global Jacobian matrix pointer array to gmgd.bin  ---
!
      WRITE(20) ((IM_GM(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Deallocate memory for global finite-element node data  ---
!
      DEALLOCATE( NF_GM,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NF_GM'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Total number of ghost finite-element nodes per processor  --- 
!
      ALLOCATE( NGHNX(1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NGHNX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Number of finite-element ghost nodes to be sent, indexed by 
!     sending direction and sending processor number, where processor  
!     number is one plus the processor rank  --- 
!
      ALLOCATE( NCSGNX(1:6,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NCSGNX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Process receiving finite-element ghost nodes to be sent, 
!     indexed by sending direction and sending processor number, 
!     where processor number is one plus the processor rank  --- 
!
      ALLOCATE( NPSGNX(1:6,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NPSGNX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Count the number finite-element ghost nodes on each processor
!     and determine the receiving processors  ---
!
      DO NP = 1,NP_MPI
        DO M = 1,6
          NCSGNX(M,NP) = 0
          NPSGNX(M,NP) = 0
        ENDDO
        NGHNX(NP) = 0
      ENDDO
      NC = 0
      NGN = 0
      IL = (IFLD+1)/IP_MPI
      JL = (JFLD+1)/JP_MPI
      KL = (KFLD+1)/KP_MPI
      DO KP = 1,KP_MPI
        IF( KP.EQ.1 ) THEN
          KD1X = 1
          KD2X = KL
        ELSEIF( KP.EQ.KP_MPI ) THEN
          KD1X = (KP-1)*KL + 1
          KD2X = KFLD + 1
        ELSE
          KD1X = (KP-1)*KL + 1
          KD2X = KP*KL
        ENDIF
        IF( KP_MPI.EQ.1 ) THEN
          KD1X = 1
          KD2X = KFLD + 1
        ENDIF
        DO JP = 1,JP_MPI
          IF( JP.EQ.1 ) THEN
            JD1X = 1
            JD2X = JL
          ELSEIF( JP.EQ.JP_MPI ) THEN
            JD1X = (JP-1)*JL + 1
            JD2X = JFLD + 1
          ELSE
            JD1X = (JP-1)*JL + 1
            JD2X = JP*JL
          ENDIF
          IF( JP_MPI.EQ.1 ) THEN
            JD1X = 1
            JD2X = JFLD + 1
          ENDIF
          DO IP = 1,IP_MPI
            IF( IP.EQ.1 ) THEN
              ID1X = 1
              ID2X = IL
            ELSEIF( IP.EQ.IP_MPI ) THEN
              ID1X = (IP-1)*IL + 1
              ID2X = IFLD + 1
            ELSE
              ID1X = (IP-1)*IL + 1
              ID2X = IP*IL
            ENDIF
            IF( IP_MPI.EQ.1 ) THEN
              ID1X = 1
              ID2X = IFLD + 1
            ENDIF
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            NCX = 0
            DO K = KD1X,KD2X
              DO J = JD1X,JD2X
                DO I = ID1X,ID2X
                  NFEN = (K-1)*(IFLD+1)*(JFLD+1) + (J-1)*(IFLD+1) + I
                  NC = NC + 1
                  NCX = NCX + 1
!
!---              Bottom FE ghost node  ---
!
                  IF( K.EQ.KD1X .AND. KP.GT.1 ) THEN
                    NCSGNX(1,NP) = NCSGNX(1,NP) + 1
                    KPB = KP - 1
                    NPB = (KPB-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
                    NPSGNX(1,NP) = NPB
                    NGN = NGN + 1
                  ENDIF
!
!---              South FE ghost node  ---
!
                  IF( J.EQ.JD1X .AND. JP.GT.1 ) THEN
                    NCSGNX(2,NP) = NCSGNX(2,NP) + 1
                    JPS = JP - 1
                    NPS = (KP-1)*IP_MPI*JP_MPI + (JPS-1)*IP_MPI + IP
                    NPSGNX(2,NP) = NPS
                    NGN = NGN + 1
                  ENDIF
!
!---              West FE ghost node  ---
!
                  IF( I.EQ.ID1X .AND. IP.GT.1 ) THEN
                    NCSGNX(3,NP) = NCSGNX(3,NP) + 1
                    IPW = IP - 1
                    NPW = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IPW
                    NPSGNX(3,NP) = NPW
                    NGN = NGN + 1
                  ENDIF
!
!---              East FE ghost node  ---
!
                  IF( I.EQ.ID2X .AND. IP.LT.IP_MPI ) THEN
                    NCSGNX(4,NP) = NCSGNX(4,NP) + 1
                    IPE = IP + 1
                    NPE = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IPE
                    NPSGNX(4,NP) = NPE
                    NGN = NGN + 1
                  ENDIF
!
!---              North FE ghost node  ---
!
                  IF( J.EQ.JD2X .AND. JP.LT.JP_MPI ) THEN
                    NCSGNX(5,NP) = NCSGNX(5,NP) + 1
                    JPN = JP + 1
                    NPN = (KP-1)*IP_MPI*JP_MPI + (JPN-1)*IP_MPI + IP
                    NPSGNX(5,NP) = NPN
                    NGN = NGN + 1
                  ENDIF
!
!---              Top FE ghost node  ---
!
                  IF( K.EQ.KD2X .AND. KP.LT.KP_MPI ) THEN
                    NCSGNX(6,NP) = NCSGNX(6,NP) + 1
                    KPT = KP + 1
                    NPT = (KPT-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
                    NPSGNX(6,NP) = NPT
                    NGN = NGN + 1
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            DO M = 1,6
              NGHNX(NP) = NGHNX(NP) + NCSGNX(M,NP)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Number of FE ghost nodes to be sent, indexed by sending direction
!     and sending processor number, where processor number is one 
!     plus the processor rank  --- 
!
      WRITE(20) ((NCSGNX(K,L),K=1,6),L=1,NP_MPI)
!
!---  Process receiving FE ghost nodes to be sent, indexed by sending 
!     direction and sending processor number, where processor number 
!     is one plus the processor rank  --- 
!
      WRITE(20) ((NPSGNX(K,L),K=1,6),L=1,NP_MPI)
!
!---  Allocate memory for local FE node numbers on 
!     processors sending FE node numbers. Ordering is ijk for bottom, 
!     south, west, east, north, top FE ghost nodes on each 
!     processor  ---
!
      ALLOCATE( NLSGNX(1:NGN),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NLSGNX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Allocate memory for local FE node numbers on processors
!     receiving FE ghost node numbers. Ordering is ijk for bottom,
!     south, west, east, north, top FE ghost nodes on each 
!     processor  ---
!
      ALLOCATE( NLRGNX(1:NGN),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NLRGNX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Initialize sending and receiving FE node number indices  ---
!
      DO N = 1,NGN
        NLSGNX(N) = 0
        NLRGNX(N) = 0
      ENDDO
!
!---  Create a listing of the local FE node numbers on processors
!     sending host FE nodes. Ordering is ijk for bottom, 
!     south, west, east, north, top FE ghost nodes on each 
!     processor  ---
!
      NGN = 0
      IL = (IFLD+1)/IP_MPI
      JL = (JFLD+1)/JP_MPI
      KL = (KFLD+1)/KP_MPI
      DO KP = 1,KP_MPI
!
!---    Starting and ending global k indices without FE ghost nodes  ---
!
        IF( KP.EQ.1 ) THEN
          K1X = 1
          K2X = KL
        ELSEIF( KP.EQ.KP_MPI ) THEN
          K1X = (KP-1)*KL + 1
          K2X = KFLD + 1
        ELSE
          K1X = (KP-1)*KL + 1
          K2X = KP*KL
        ENDIF
        IF( KP_MPI.EQ.1 ) THEN
          K1X = 1
          K2X = KFLD + 1
        ENDIF
!
!---    Starting and ending global k indices with FE ghost nodes  ---
!
        IF( KP.EQ.1 ) THEN
          K3X = 1
          K4X = KL + 1
        ELSEIF( KP.EQ.KP_MPI ) THEN
          K3X = (KP-1)*KL
          K4X = KFLD + 1
        ELSE
          K3X = (KP-1)*KL
          K4X = KP*KL + 1
        ENDIF
        IF( KP_MPI.EQ.1 ) THEN
          K3X = 1
          K4X = KFLD + 1
        ENDIF
        DO JP = 1,JP_MPI
!
!---      Starting and ending global j indices without FE 
!         ghost nodes  ---
!
          IF( JP.EQ.1 ) THEN
            J1X = 1
            J2X = JL
          ELSEIF( JP.EQ.JP_MPI ) THEN
            J1X = (JP-1)*JL + 1
            J2X = JFLD + 1
          ELSE
            J1X = (JP-1)*JL + 1
            J2X = JP*JL
          ENDIF
          IF( JP_MPI.EQ.1 ) THEN
            J1X = 1
            J2X = JFLD + 1
          ENDIF
!
!---      Starting and ending global j indices with FE ghost nodes  ---
!
          IF( JP.EQ.1 ) THEN
            J3X = 1
            J4X = JL + 1
          ELSEIF( JP.EQ.JP_MPI ) THEN
            J3X = (JP-1)*JL
            J4X = JFLD + 1
          ELSE
            J3X = (JP-1)*JL
            J4X = JP*JL + 1
          ENDIF
          IF( JP_MPI.EQ.1 ) THEN
            J3X = 1
            J4X = JFLD + 1
          ENDIF
          DO IP = 1,IP_MPI
!
!---        Starting and ending global i indices without FE ghost 
!           nodes  ---
!
            IF( IP.EQ.1 ) THEN
              I1X = 1
              I2X = IL
            ELSEIF( IP.EQ.IP_MPI ) THEN
              I1X = (IP-1)*IL + 1
              I2X = IFLD + 1
            ELSE
              I1X = (IP-1)*IL + 1
              I2X = IP*IL
            ENDIF
            IF( IP_MPI.EQ.1 ) THEN
              I1X = 1
              I2X = IFLD + 1
            ENDIF
!
!---        Starting and ending global i indices with FE 
!           ghost nodes  ---
!
            IF( IP.EQ.1 ) THEN
              I3X = 1
              I4X = IL + 1
            ELSEIF( IP.EQ.IP_MPI ) THEN
              I3X = (IP-1)*IL
              I4X = IFLD + 1
            ELSE
              I3X = (IP-1)*IL
              I4X = IP*IL + 1
            ENDIF
            IF( IP_MPI.EQ.1 ) THEN
              I3X = 1
              I4X = IFLD + 1
            ENDIF
!
!---        Sending bottom tier of FE nodes on top processor
!           to top tier of FE ghost nodes on bottom processor, with
!           NLSGNX holding sending FE nodes  ---
!
            IF( KP.GT.1 ) THEN
              KGX = K3X + 1
              DO JGX = J1X,J2X
                DO IGX = I1X,I2X
                  NGN = NGN + 1
                  NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                  NLSGNX(NGN) = NLSX
                ENDDO
              ENDDO
            ENDIF
!
!---        Sending south tier of FE nodes on north processor
!           to north tier of FE ghost nodes on south processor, with
!           NLSGNX holding sending FE nodes  ---
!
            IF( JP.GT.1 ) THEN
              JGX = J3X + 1
              DO KGX = K1X,K2X
                DO IGX = I1X,I2X
                  NGN = NGN + 1
                  NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                  NLSGNX(NGN) = NLSX
                ENDDO
              ENDDO
            ENDIF
!
!---        Sending west tier of FE nodes on east processor
!           to east tier of FE ghost nodes on west processor, with
!           NLSGNX holding sending FE nodes  ---
!
            IF( IP.GT.1 ) THEN
              IGX = I3X + 1
              DO KGX = K1X,K2X
                DO JGX = J1X,J2X
                  NGN = NGN + 1
                  NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                  NLSGNX(NGN) = NLSX
                ENDDO
              ENDDO
            ENDIF
!
!---        Sending east tier of FE nodes on west processor
!           to west tier of FE ghost nodes on east processor, with
!           NLSGNX holding sending FE nodes  ---
!
            IF( IP.LT.IP_MPI ) THEN
!
!---          Starting and ending global i indices with FE 
!             ghost nodes  ---
!
              IGX = I4X - 1
              DO KGX = K1X,K2X
                DO JGX = J1X,J2X
                  NGN = NGN + 1
                  NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                   NLSGNX(NGN) = NLSX
                ENDDO
              ENDDO
            ENDIF
!
!---        Sending north tier of FE nodes on south processor
!           to south tier of FE ghost nodes on north processor, with
!           NLSGNX holding sending FE nodes  ---
!
            IF( JP.LT.JP_MPI ) THEN
!
!---          Starting and ending global j indices with FE ghost 
!             nodes  ---
!
              JGX = J4X - 1
              DO KGX = K1X,K2X
                DO IGX = I1X,I2X
                  NGN = NGN + 1
                  NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                  NLSGNX(NGN) = NLSX
                ENDDO
              ENDDO
            ENDIF
!
!---        Sending top tier of FE nodes on bottom processor
!           to bottom tier of FE ghost nodes on top processor, with
!           NLSGNX holding sending FE nodes  ---
!
            IF( KP.LT.KP_MPI ) THEN
!
!---          Starting and ending global k indices with FE 
!             ghost nodes  ---
!
              KGX = K4X - 1
              DO JGX = J1X,J2X
                DO IGX = I1X,I2X
                  NGN = NGN + 1
                  NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                  NLSGNX(NGN) = NLSX
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Create a listing of the local FE node numbers on processors
!     receiving FE ghost nodes. Ordering is ijk for bottom, 
!     south, west, east, north, top FE ghost nodes on each 
!     processor  ---
!
      NGN = 0
      IL = (IFLD+1)/IP_MPI
      JL = (JFLD+1)/JP_MPI
      KL = (KFLD+1)/KP_MPI
      DO KP = 1,KP_MPI
!
!---    Starting and ending global k indices without FE ghost nodes  ---
!
        IF( KP.EQ.1 ) THEN
          K1X = 1
          K2X = KL
        ELSEIF( KP.EQ.KP_MPI ) THEN
          K1X = (KP-1)*KL + 1
          K2X = KFLD + 1
        ELSE
          K1X = (KP-1)*KL + 1
          K2X = KP*KL
        ENDIF
        IF( KP_MPI.EQ.1 ) THEN
          K1X = 1
          K2X = KFLD + 1
        ENDIF
!
!---    Starting and ending global k indices with FE ghost nodes  ---
!
        IF( KP.EQ.1 ) THEN
          K3X = 1
          K4X = KL + 1
        ELSEIF( KP.EQ.KP_MPI ) THEN
          K3X = (KP-1)*KL
          K4X = KFLD + 1
        ELSE
          K3X = (KP-1)*KL
          K4X = KP*KL + 1
        ENDIF
        IF( KP_MPI.EQ.1 ) THEN
          K3X = 1
          K4X = KFLD + 1
        ENDIF
        DO JP = 1,JP_MPI
!
!---      Starting and ending global j indices without FE ghost
!         nodes  ---
!
          IF( JP.EQ.1 ) THEN
            J1X = 1
            J2X = JL
          ELSEIF( JP.EQ.JP_MPI ) THEN
            J1X = (JP-1)*JL + 1
            J2X = JFLD + 1
          ELSE
            J1X = (JP-1)*JL + 1
            J2X = JP*JL
          ENDIF
          IF( JP_MPI.EQ.1 ) THEN
            J1X = 1
            J2X = JFLD + 1
          ENDIF
!
!---      Starting and ending global j indices with FE ghost nodes  ---
!
          IF( JP.EQ.1 ) THEN
            J3X = 1
            J4X = JL + 1
          ELSEIF( JP.EQ.JP_MPI ) THEN
            J3X = (JP-1)*JL
            J4X = JFLD + 1
          ELSE
            J3X = (JP-1)*JL
            J4X = JP*JL + 1
          ENDIF
          IF( JP_MPI.EQ.1 ) THEN
            J3X = 1
            J4X = JFLD + 1
          ENDIF
          DO IP = 1,IP_MPI
!
!---        Starting and ending global i indices without FE ghost 
!           nodes  ---
!
            IF( IP.EQ.1 ) THEN
              I1X = 1
              I2X = IL
            ELSEIF( IP.EQ.IP_MPI ) THEN
              I1X = (IP-1)*IL + 1
              I2X = IFLD + 1
            ELSE
              I1X = (IP-1)*IL + 1
              I2X = IP*IL
            ENDIF
            IF( IP_MPI.EQ.1 ) THEN
              I1X = 1
              I2X = IFLD + 1
            ENDIF
!
!---        Starting and ending global i indices with FE 
!           ghost nodes  ---
!
            IF( IP.EQ.1 ) THEN
              I3X = 1
              I4X = IL + 1
            ELSEIF( IP.EQ.IP_MPI ) THEN
              I3X = (IP-1)*IL
              I4X = IFLD + 1
            ELSE
              I3X = (IP-1)*IL
              I4X = IP*IL + 1
            ENDIF
            IF( IP_MPI.EQ.1 ) THEN
              I3X = 1
              I4X = IFLD + 1
            ENDIF
!
!---        Receiving bottom tier of FE nodes on top processor
!           from top tier of FE ghost nodes on bottom processor, with
!           NLRGNX holding receiving FE ghost nodes  ---
!
            IF( KP.GT.1 ) THEN
              KGX = K3X
              DO JGX = J1X,J2X
                DO IGX = I1X,I2X
                  NGN = NGN + 1
                  NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                  NLRGNX(NGN) = NLRX
                ENDDO
              ENDDO
            ENDIF
!
!---        Receiving south tier of FE nodes on north processor
!           from north tier of FE ghost nodes on south processor, with
!           NLRGNX holding receiving FE ghost nodes  ---
!
            IF( JP.GT.1 ) THEN
              JGX = J3X
              DO KGX = K1X,K2X
                DO IGX = I1X,I2X
                  NGN = NGN + 1
                  NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                  NLRGNX(NGN) = NLRX
                ENDDO
              ENDDO
            ENDIF
!
!---        Receiving west tier of FE nodes on east processor
!           from east tier of FE ghost nodes on west processor, with
!           NLRGNX holding receiving FE ghost nodes  ---
!
            IF( IP.GT.1 ) THEN
              IGX = I3X
              DO KGX = K1X,K2X
                DO JGX = J1X,J2X
                  NGN = NGN + 1
                  NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                  NLRGNX(NGN) = NLRX
                ENDDO
              ENDDO
            ENDIF
!
!---        Receiving east tier of FE nodes on west processor
!           from west tier of FE ghost nodes on east processor, with
!           NLRGNX holding receiving FE ghost nodes  ---
!
            IF( IP.LT.IP_MPI ) THEN
!
!---          Starting and ending global i indices with FE 
!             ghost nodes  ---
!
              IGX = I4X
              DO KGX = K1X,K2X
                DO JGX = J1X,J2X
                  NGN = NGN + 1
                  NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                  NLRGNX(NGN) = NLRX
                ENDDO
              ENDDO
            ENDIF
!
!---        Receiving north tier of FE nodes on south processor
!           from south tier of FE ghost nodes on north processor, with
!           NLRGNX holding receiving FE ghost nodes  ---
!
            IF( JP.LT.JP_MPI ) THEN
!
!---          Starting and ending global j indices with FE 
!             ghost nodes  ---
!
              JGX = J4X
              DO KGX = K1X,K2X
                DO IGX = I1X,I2X
                  NGN = NGN + 1
                  NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                  NLRGNX(NGN) = NLRX
                ENDDO
              ENDDO
            ENDIF
!
!---        Receiving top tier of FE nodes on bottom processor
!           from bottom tier of FE ghost nodes on bottom processor, with
!           NLRGNX holding receiving FE ghost nodes  ---
!
            IF( KP.LT.KP_MPI ) THEN
!
!---          Starting and ending global k indices with FE 
!             ghost nodes  ---
!
              KGX = K4X
              DO JGX = J1X,J2X
                DO IGX = I1X,I2X
                  NGN = NGN + 1
                  NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &              (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                  NLRGNX(NGN) = NLRX
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Create a listing of the local FE node numbers on processors
!     sending and receiving FE ghost nodes. Ordering is ijk for bottom, 
!     south, west, east, north, top FE ghost nodes on each 
!     processor  ---
!
      WRITE(20) (NLSGNX(L),L=1,NGN)
      WRITE(20) (NLRGNX(L),L=1,NGN)
!
!---  Deallocate NCSGNX  --- 
!
      DEALLOCATE( NCSGNX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NCSGNX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Deallocate NPSGNX  --- 
!
      DEALLOCATE( NPSGNX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NPSGNX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Deallocate NLRGNX  --- 
!
      DEALLOCATE( NLRGNX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NLRGNX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Deallocate NLSGNX  --- 
!
      DEALLOCATE( NLSGNX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NLSGNX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Close the gmgd.bin file  ---
!
      CLOSE(20)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_GMGD_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_GRID_CO2( IDPX,JDPX,KDPX )
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
!     Write grid data to a series of binary grid files.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
      USE CONST
      USE BCVP
      USE BCVS
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: GRVXX,GRVYX,GRVZX
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: AFXX,AFYX,AFZX
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: DXGPX,DYGPX,DZGPX
      REAL*8, DIMENSION(:), ALLOCATABLE :: RPX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IXPX,NDX
      INTEGER, DIMENSION(1:2,1:NP_MPI) :: IDPX,JDPX,KDPX
      INTEGER, DIMENSION(:), ALLOCATABLE :: NUKFLX,NUKFOX
      INTEGER, DIMENSION(:), ALLOCATABLE :: NUKTLX,NUKTOX
      INTEGER, DIMENSION(:), ALLOCATABLE :: NGHCX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NCSGCX,NPSGCX
      INTEGER, DIMENSION(:), ALLOCATABLE :: NLSGCX,NLRGCX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRITE_GRID_CO2'
!
!---  Open a new grid1.bin file  --
!
      OPEN( UNIT=20,FILE='grid1.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='grid1.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write grid data to the grid1.bin file--
!
      WRITE(20) (NC_MPI(NP),NP=1,NP_MPI)
      WRITE(20) ((IDPX(M,NP),M=1,2),NP=1,NP_MPI)
      WRITE(20) ((JDPX(M,NP),M=1,2),NP=1,NP_MPI)
      WRITE(20) ((KDPX(M,NP),M=1,2),NP=1,NP_MPI)
!
!---  Write grid data to the grid1.bin file--
!
      WRITE(20) (((XE(M,ND_MPI(K,L)),M=1,8),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Close the grid1.bin file  ---
!
      CLOSE(20)
!
!---  Open a new grid2.bin file  --
!
      OPEN( UNIT=20,FILE='grid2.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='grid2.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
      WRITE(20) (((YE(M,ND_MPI(K,L)),M=1,8),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Close the grid2.bin file  ---
!
      CLOSE(20)
!
!---  Open a new grid3.bin file  --
!
      OPEN( UNIT=20,FILE='grid3.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='grid3.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
      WRITE(20) (((ZE(M,ND_MPI(K,L)),M=1,8),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Close the grid3.bin file  ---
!
      CLOSE(20)
!
!---  Open a new grid4.bin file  --
!
      OPEN( UNIT=20,FILE='grid4.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='grid4.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write grid data to the grid4.bin file--
!
      WRITE(20) ((XP(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((YP(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((ZP(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Allocate temporary memory  ---
!
      ALLOCATE( RPX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: RPX'
        CALL WRMSGS( INDX )
      ENDIF
      DO N = 1,NFLD
        RPX(N) = RP(ID(N))
      ENDDO
      WRITE(20) ((RPX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      DEALLOCATE( RPX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: RPX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Close the grid4.bin file  ---
!
      CLOSE(20)
!
!---  Open a new grid5.bin file  --
!
      OPEN( UNIT=20,FILE='grid5.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='grid5.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write grid data to the grid5.bin file--
!
      WRITE(20) ((DXGF(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((DYGF(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((DZGF(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Close the grid5.bin file  ---
!
      CLOSE(20)
!
!---  Open a new grid6.bin file  --
!
      OPEN( UNIT=20,FILE='grid6.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='grid6.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Allocate temporary memory and write grid data to grid6.bin  ---
!
      ALLOCATE( DXGPX(1:2,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: DXGPX'
        CALL WRMSGS( INDX )
      ENDIF
      DO N = 1,NFLD
        NPX = NSX(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        DXGPX(1,N) = DXGP(NPX)
        DXGPX(2,N) = DXGP(NQX)
      ENDDO
      WRITE(20) (((DXGPX(M,ND_MPI(K,L)),M=1,2),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      DEALLOCATE( DXGPX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: DXGPX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( DYGPX(1:2,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: DYGPX'
        CALL WRMSGS( INDX )
      ENDIF
      DO N = 1,NFLD
        NPY = NSY(N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        DYGPX(1,N) = DYGP(NPY)
        DYGPX(2,N) = DYGP(NQY)
      ENDDO
      WRITE(20) (((DYGPX(M,ND_MPI(K,L)),M=1,2),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      DEALLOCATE( DYGPX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: DYGPX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( DZGPX(1:2,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: DZGPX'
        CALL WRMSGS( INDX )
      ENDIF
      DO N = 1,NFLD
        NPZ = NSZ(N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        DZGPX(1,N) = DZGP(NPZ)
        DZGPX(2,N) = DZGP(NQZ)
      ENDDO
      WRITE(20) (((DZGPX(M,ND_MPI(K,L)),M=1,2),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      DEALLOCATE( DZGPX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: DZGPX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Close the grid6.bin file  ---
!
      CLOSE(20)
!
!---  Open a new grid7.bin file  --
!
      OPEN( UNIT=20,FILE='grid7.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='grid7.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Allocate temporary memory and write grid data to grid7.bin  ---
!
      ALLOCATE( AFXX(1:2,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: AFXX'
        CALL WRMSGS( INDX )
      ENDIF
      DO N = 1,NFLD
        NPX = NSX(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        AFXX(1,N) = AFX(NPX)
        AFXX(2,N) = AFX(NQX)
      ENDDO
      WRITE(20) (((AFXX(M,ND_MPI(K,L)),M=1,2),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      DEALLOCATE( AFXX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: AFXX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( AFYX(1:2,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: AFYX'
        CALL WRMSGS( INDX )
      ENDIF
      DO N = 1,NFLD
        NPY = NSY(N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        AFYX(1,N) = AFY(NPY)
        AFYX(2,N) = AFY(NQY)
      ENDDO
      WRITE(20) (((AFYX(M,ND_MPI(K,L)),M=1,2),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      DEALLOCATE( AFYX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: AFYX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( AFZX(1:2,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: AFZX'
        CALL WRMSGS( INDX )
      ENDIF
      DO N = 1,NFLD
        NPZ = NSZ(N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        AFZX(1,N) = AFZ(NPZ)
        AFZX(2,N) = AFZ(NQZ)
      ENDDO
      WRITE(20) (((AFZX(M,ND_MPI(K,L)),M=1,2),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      DEALLOCATE( AFZX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: AFZX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Close the grid7.bin file  ---
!
      CLOSE(20)
!
!---  Open a new grid8.bin file  --
!
      OPEN( UNIT=20,FILE='grid8.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='grid8.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Allocate temporary memory and write grid data to grid8.bin  ---
!
      WRITE(20) ((VOL(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      ALLOCATE( GRVXX(1:2,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: GRVXX'
        CALL WRMSGS( INDX )
      ENDIF
      DO N = 1,NFLD
        NPX = NSX(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        GRVXX(1,N) = GRVX(NPX)
        GRVXX(2,N) = GRVX(NQX)
      ENDDO
      WRITE(20) (((GRVXX(M,ND_MPI(K,L)),M=1,2),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      DEALLOCATE( GRVXX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: GRVXX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( GRVYX(1:2,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: GRVYX'
        CALL WRMSGS( INDX )
      ENDIF
      DO N = 1,NFLD
        NPY = NSY(N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        GRVYX(1,N) = GRVY(NPY)
        GRVYX(2,N) = GRVY(NQY)
      ENDDO
      WRITE(20) (((GRVYX(M,ND_MPI(K,L)),M=1,2),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      DEALLOCATE( GRVYX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: GRVYX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( GRVZX(1:2,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: GRVZX'
        CALL WRMSGS( INDX )
      ENDIF
      DO N = 1,NFLD
        NPZ = NSZ(N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        GRVZX(1,N) = GRVZ(NPZ)
        GRVZX(2,N) = GRVZ(NQZ)
      ENDDO
      WRITE(20) (((GRVZX(M,ND_MPI(K,L)),M=1,2),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      DEALLOCATE( GRVZX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: GRVZX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Close the grid8.bin file  ---
!
      CLOSE(20)
!
!---  Open a new grid9.bin file  --
!
      OPEN( UNIT=20,FILE='grid9.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='grid9.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Allocate temporary memory and write grid data to grid9.bin  ---
!
      ALLOCATE( IXPX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IXPX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( NDX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NDX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( NUKFLX(1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NUKFLX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( NUKFOX(1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NUKFOX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( NUKTLX(1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NUKTLX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( NUKTOX(1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NUKTOX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Total number of ghost cells per processor  --- 
!
      ALLOCATE( NGHCX(1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NGHCX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Number of ghost cells to be sent, indexed by sending direction
!     and sending processor number, where processor number is one 
!     plus the processor rank  --- 
!
      ALLOCATE( NCSGCX(1:6,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NCSGCX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Process receiving ghost cells to be sent, indexed by sending 
!     direction and sending processor number, where processor number 
!     is one plus the processor rank  --- 
!
      ALLOCATE( NPSGCX(1:6,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NPSGCX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Create a mapping of equation ordering, with nodes
!     counted by processor, then i,j,k indexing, not including
!     ghost cells and count the number ghost cells on each 
!     processor  ---
!
      DO NP = 1,NP_MPI
        DO M = 1,6
          NCSGCX(M,NP) = 0
          NPSGCX(M,NP) = 0
        ENDDO
        NGHCX(NP) = 0
      ENDDO
      NC = 0
      NGC = 0
      NUKFGX = 0
      NUKTGX = 0
      IL = IFLD/IP_MPI
      JL = JFLD/JP_MPI
      KL = KFLD/KP_MPI
      DO KP = 1,KP_MPI
        IF( KP.EQ.1 ) THEN
          KD1X = 1
          KD2X = KL
        ELSEIF( KP.EQ.KP_MPI ) THEN
          KD1X = (KP-1)*KL + 1
          KD2X = KFLD
        ELSE
          KD1X = (KP-1)*KL + 1
          KD2X = KP*KL
        ENDIF
        IF( KP_MPI.EQ.1 ) THEN
          KD1X = 1
          KD2X = KFLD
        ENDIF
        DO JP = 1,JP_MPI
          IF( JP.EQ.1 ) THEN
            JD1X = 1
            JD2X = JL
          ELSEIF( JP.EQ.JP_MPI ) THEN
            JD1X = (JP-1)*JL + 1
            JD2X = JFLD
          ELSE
            JD1X = (JP-1)*JL + 1
            JD2X = JP*JL
          ENDIF
          IF( JP_MPI.EQ.1 ) THEN
            JD1X = 1
            JD2X = JFLD
          ENDIF
          DO IP = 1,IP_MPI
            IF( IP.EQ.1 ) THEN
              ID1X = 1
              ID2X = IL
            ELSEIF( IP.EQ.IP_MPI ) THEN
              ID1X = (IP-1)*IL + 1
              ID2X = IFLD
            ELSE
              ID1X = (IP-1)*IL + 1
              ID2X = IP*IL
            ENDIF
            IF( IP_MPI.EQ.1 ) THEN
              ID1X = 1
              ID2X = IFLD
            ENDIF
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            NCX = 0
            DO K = KD1X,KD2X
              DO J = JD1X,JD2X
                DO I = ID1X,ID2X
                  N = ND(I,J,K)
                  IF( IXP(N).EQ.0 ) THEN
                    IXPX(N) = 0
                  ELSE
                    NC = NC + 1
                    NCX = NCX + 1
                    IXPX(N) = NC
!
!---                Bottom FE ghost node  ---
!
                    IF( K.EQ.KD1X .AND. KP.GT.1 ) THEN
                      NCSGCX(1,NP) = NCSGCX(1,NP) + 1
                      KPB = KP - 1
                      NPB = (KPB-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
                      NPSGCX(1,NP) = NPB
                      NGC = NGC + 1
                    ENDIF
!
!---                South FE ghost node  ---
!
                    IF( J.EQ.JD1X .AND. JP.GT.1 ) THEN
                      NCSGCX(2,NP) = NCSGCX(2,NP) + 1
                      JPS = JP - 1
                      NPS = (KP-1)*IP_MPI*JP_MPI + (JPS-1)*IP_MPI + IP
                      NPSGCX(2,NP) = NPS
                      NGC = NGC + 1
                    ENDIF
!
!---                West FE ghost node  ---
!
                    IF( I.EQ.ID1X .AND. IP.GT.1 ) THEN
                      NCSGCX(3,NP) = NCSGCX(3,NP) + 1
                      IPW = IP - 1
                      NPW = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IPW
                      NPSGCX(3,NP) = NPW
                      NGC = NGC + 1
                    ENDIF
!
!---                East FE ghost node  ---
!
                    IF( I.EQ.ID2X .AND. IP.LT.IP_MPI ) THEN
                      NCSGCX(4,NP) = NCSGCX(4,NP) + 1
                      IPE = IP + 1
                      NPE = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IPE
                      NPSGCX(4,NP) = NPE
                      NGC = NGC + 1
                    ENDIF
!
!---                North FE ghost node  ---
!
                    IF( J.EQ.JD2X .AND. JP.LT.JP_MPI ) THEN
                      NCSGCX(5,NP) = NCSGCX(5,NP) + 1
                      JPN = JP + 1
                      NPN = (KP-1)*IP_MPI*JP_MPI + (JPN-1)*IP_MPI + IP
                      NPSGCX(5,NP) = NPN
                      NGC = NGC + 1
                    ENDIF
!
!---                Top FE ghost node  ---
!
                    IF( K.EQ.KD2X .AND. KP.LT.KP_MPI ) THEN
                      NCSGCX(6,NP) = NCSGCX(6,NP) + 1
                      KPT = KP + 1
                      NPT = (KPT-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
                      NPSGCX(6,NP) = NPT
                      NGC = NGC + 1
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            NUKFLX(NP) = NCX*ISVC
            NUKFOX(NP) = (NC-NCX)*ISVC
            NUKFGX = NUKFGX + NCX*ISVC
            NUKTLX(NP) = NCX
            NUKTOX(NP) = (NC-NCX)
            NUKTGX = NUKTGX + NCX
            DO M = 1,6
              NGHCX(NP) = NGHCX(NP) + NCSGCX(M,NP)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Locate coupled-well equations on the last processor  ---
!
      NUKFGX = NUKFGX + N_CW
      NUKFLX(NP_MPI) = NUKFLX(NP_MPI) + N_CW
!
!---  Allocate memory for local node numbers on processors
!     sending ghost cells. Ordering is ijk for bottom, south,
!     west, east, north, top ghost cells on each processor  ---
!
      ALLOCATE( NLSGCX(1:NGC),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NLSGCX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Allocate memory for local node numbers on processors
!     receiving ghost cells. Ordering is ijk for bottom, south,
!     west, east, north, top ghost cells on each processor  ---
!
      ALLOCATE( NLRGCX(1:NGC),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NLRGCX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Create a listing of the local node numbers on processors
!     sending host cells. Ordering is ijk for bottom, 
!     south, west, east, north, top ghost cells on each processor  ---
!
      NGC = 0
      IL = IFLD/IP_MPI
      JL = JFLD/JP_MPI
      KL = KFLD/KP_MPI
      DO KP = 1,KP_MPI
!
!---    Starting and ending global k indices without ghost cells  ---
!
        IF( KP.EQ.1 ) THEN
          K1X = 1
          K2X = KL
        ELSEIF( KP.EQ.KP_MPI ) THEN
          K1X = (KP-1)*KL + 1
          K2X = KFLD
        ELSE
          K1X = (KP-1)*KL + 1
          K2X = KP*KL
        ENDIF
        IF( KP_MPI.EQ.1 ) THEN
          K1X = 1
          K2X = KFLD
        ENDIF
!
!---    Starting and ending global k indices with ghost cells  ---
!
        IF( KP.EQ.1 ) THEN
          K3X = 1
          K4X = KL + 1
        ELSEIF( KP.EQ.KP_MPI ) THEN
          K3X = (KP-1)*KL
          K4X = KFLD
        ELSE
          K3X = (KP-1)*KL
          K4X = KP*KL + 1
        ENDIF
        IF( KP_MPI.EQ.1 ) THEN
          K3X = 1
          K4X = KFLD
        ENDIF
        DO JP = 1,JP_MPI
!
!---      Starting and ending global j indices without ghost cells  ---
!
          IF( JP.EQ.1 ) THEN
            J1X = 1
            J2X = JL
          ELSEIF( JP.EQ.JP_MPI ) THEN
            J1X = (JP-1)*JL + 1
            J2X = JFLD
          ELSE
            J1X = (JP-1)*JL + 1
            J2X = JP*JL
          ENDIF
          IF( JP_MPI.EQ.1 ) THEN
            J1X = 1
            J2X = JFLD
          ENDIF
!
!---      Starting and ending global j indices with ghost cells  ---
!
          IF( JP.EQ.1 ) THEN
            J3X = 1
            J4X = JL + 1
          ELSEIF( JP.EQ.JP_MPI ) THEN
            J3X = (JP-1)*JL
            J4X = JFLD
          ELSE
            J3X = (JP-1)*JL
            J4X = JP*JL + 1
          ENDIF
          IF( JP_MPI.EQ.1 ) THEN
            J3X = 1
            J4X = JFLD
          ENDIF
          DO IP = 1,IP_MPI
!
!---        Starting and ending global i indices without ghost 
!           cells  ---
!
            IF( IP.EQ.1 ) THEN
              I1X = 1
              I2X = IL
            ELSEIF( IP.EQ.IP_MPI ) THEN
              I1X = (IP-1)*IL + 1
              I2X = IFLD
            ELSE
              I1X = (IP-1)*IL + 1
              I2X = IP*IL
            ENDIF
            IF( IP_MPI.EQ.1 ) THEN
              I1X = 1
              I2X = IFLD
            ENDIF
!
!---        Starting and ending global i indices with ghost cells  ---
!
            IF( IP.EQ.1 ) THEN
              I3X = 1
              I4X = IL + 1
            ELSEIF( IP.EQ.IP_MPI ) THEN
              I3X = (IP-1)*IL
              I4X = IFLD
            ELSE
              I3X = (IP-1)*IL
              I4X = IP*IL + 1
            ENDIF
            IF( IP_MPI.EQ.1 ) THEN
              I3X = 1
              I4X = IFLD
            ENDIF
!
!---        Sending bottom tier of field cells on top processor
!           to top tier of ghost cells on bottom processor, with
!           NLSGCX holding sending field cells  ---
!
            IF( KP.GT.1 ) THEN
              KGX = K3X + 1
              DO JGX = J1X,J2X
                DO IGX = I1X,I2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLSGCX(NGC) = NLSX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
!---        Sending south tier of field cells on north processor
!           to north tier of ghost cells on south processor, with
!           NLSGCX holding sending field cells  ---
!
            IF( JP.GT.1 ) THEN
              JGX = J3X + 1
              DO KGX = K1X,K2X
                DO IGX = I1X,I2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLSGCX(NGC) = NLSX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
!---        Sending west tier of field cells on east processor
!           to east tier of ghost cells on west processor, with
!           NLSGCX holding sending field cells  ---
!
            IF( IP.GT.1 ) THEN
              IGX = I3X + 1
              DO KGX = K1X,K2X
                DO JGX = J1X,J2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLSGCX(NGC) = NLSX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
!---        Sending east tier of field cells on west processor
!           to west tier of ghost cells on east processor, with
!           NLSGCX holding sending field cells  ---
!
            IF( IP.LT.IP_MPI ) THEN
!
!---          Starting and ending global i indices with ghost cells  ---
!
              IGX = I4X - 1
              DO KGX = K1X,K2X
                DO JGX = J1X,J2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLSGCX(NGC) = NLSX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
!---        Sending north tier of field cells on south processor
!           to south tier of ghost cells on north processor, with
!           NLSGCX holding sending field cells  ---
!
            IF( JP.LT.JP_MPI ) THEN
!
!---          Starting and ending global j indices with ghost cells  ---
!
              JGX = J4X - 1
              DO KGX = K1X,K2X
                DO IGX = I1X,I2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLSGCX(NGC) = NLSX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
!---        Sending top tier of field cells on bottom processor
!           to bottom tier of ghost cells on top processor, with
!           NLSGCX holding sending field cells  ---
!
            IF( KP.LT.KP_MPI ) THEN
!
!---          Starting and ending global k indices with ghost cells  ---
!
              KGX = K4X - 1
              DO JGX = J1X,J2X
                DO IGX = I1X,I2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLSX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLSGCX(NGC) = NLSX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Create a listing of the local node numbers on processors
!     receiving ghost cells. Ordering is ijk for bottom, 
!     south, west, east, north, top ghost cells on each processor  ---
!
      NGC = 0
      IL = IFLD/IP_MPI
      JL = JFLD/JP_MPI
      KL = KFLD/KP_MPI
      DO KP = 1,KP_MPI
!
!---    Starting and ending global k indices without ghost cells  ---
!
        IF( KP.EQ.1 ) THEN
          K1X = 1
          K2X = KL
        ELSEIF( KP.EQ.KP_MPI ) THEN
          K1X = (KP-1)*KL + 1
          K2X = KFLD
        ELSE
          K1X = (KP-1)*KL + 1
          K2X = KP*KL
        ENDIF
        IF( KP_MPI.EQ.1 ) THEN
          K1X = 1
          K2X = KFLD
        ENDIF
!
!---    Starting and ending global k indices with ghost cells  ---
!
        IF( KP.EQ.1 ) THEN
          K3X = 1
          K4X = KL + 1
        ELSEIF( KP.EQ.KP_MPI ) THEN
          K3X = (KP-1)*KL
          K4X = KFLD
        ELSE
          K3X = (KP-1)*KL
          K4X = KP*KL + 1
        ENDIF
        IF( KP_MPI.EQ.1 ) THEN
          K3X = 1
          K4X = KFLD
        ENDIF
        DO JP = 1,JP_MPI
!
!---      Starting and ending global j indices without ghost cells  ---
!
          IF( JP.EQ.1 ) THEN
            J1X = 1
            J2X = JL
          ELSEIF( JP.EQ.JP_MPI ) THEN
            J1X = (JP-1)*JL + 1
            J2X = JFLD
          ELSE
            J1X = (JP-1)*JL + 1
            J2X = JP*JL
          ENDIF
          IF( JP_MPI.EQ.1 ) THEN
            J1X = 1
            J2X = JFLD
          ENDIF
!
!---      Starting and ending global j indices with ghost cells  ---
!
          IF( JP.EQ.1 ) THEN
            J3X = 1
            J4X = JL + 1
          ELSEIF( JP.EQ.JP_MPI ) THEN
            J3X = (JP-1)*JL
            J4X = JFLD
          ELSE
            J3X = (JP-1)*JL
            J4X = JP*JL + 1
          ENDIF
          IF( JP_MPI.EQ.1 ) THEN
            J3X = 1
            J4X = JFLD
          ENDIF
          DO IP = 1,IP_MPI
!
!---        Starting and ending global i indices without ghost 
!           cells  ---
!
            IF( IP.EQ.1 ) THEN
              I1X = 1
              I2X = IL
            ELSEIF( IP.EQ.IP_MPI ) THEN
              I1X = (IP-1)*IL + 1
              I2X = IFLD
            ELSE
              I1X = (IP-1)*IL + 1
              I2X = IP*IL
            ENDIF
            IF( IP_MPI.EQ.1 ) THEN
              I1X = 1
              I2X = IFLD
            ENDIF
!
!---        Starting and ending global i indices with ghost cells  ---
!
            IF( IP.EQ.1 ) THEN
              I3X = 1
              I4X = IL + 1
            ELSEIF( IP.EQ.IP_MPI ) THEN
              I3X = (IP-1)*IL
              I4X = IFLD
            ELSE
              I3X = (IP-1)*IL
              I4X = IP*IL + 1
            ENDIF
            IF( IP_MPI.EQ.1 ) THEN
              I3X = 1
              I4X = IFLD
            ENDIF
!
!---        Receiving bottom tier of field cells on top processor
!           from top tier of ghost cells on bottom processor, with
!           NLRGCX holding receiving ghost cells  ---
!
            IF( KP.GT.1 ) THEN
              KGX = K3X
              DO JGX = J1X,J2X
                DO IGX = I1X,I2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLRGCX(NGC) = NLRX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
!---        Receiving south tier of field cells on north processor
!           from north tier of ghost cells on south processor, with
!           NLRGCX holding receiving ghost cells  ---
!
            IF( JP.GT.1 ) THEN
              JGX = J3X
              DO KGX = K1X,K2X
                DO IGX = I1X,I2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLRGCX(NGC) = NLRX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
!---        Receiving west tier of field cells on east processor
!           from east tier of ghost cells on west processor, with
!           NLRGCX holding receiving ghost cells  ---
!
            IF( IP.GT.1 ) THEN
              IGX = I3X
              DO KGX = K1X,K2X
                DO JGX = J1X,J2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLRGCX(NGC) = NLRX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
!---        Receiving east tier of field cells on west processor
!           from west tier of ghost cells on east processor, with
!           NLRGCX holding receiving ghost cells  ---
!
            IF( IP.LT.IP_MPI ) THEN
!
!---          Starting and ending global i indices with ghost cells  ---
!
              IGX = I4X
              DO KGX = K1X,K2X
                DO JGX = J1X,J2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLRGCX(NGC) = NLRX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
!---        Receiving north tier of field cells on south processor
!           from south tier of ghost cells on north processor, with
!           NLRGCX holding receiving ghost cells  ---
!
            IF( JP.LT.JP_MPI ) THEN
!
!---          Starting and ending global j indices with ghost cells  ---
!
              JGX = J4X
              DO KGX = K1X,K2X
                DO IGX = I1X,I2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLRGCX(NGC) = NLRX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
!---        Receiving top tier of field cells on bottom processor
!           from bottom tier of ghost cells on bottom processor, with
!           NLRGCX holding receiving ghost cells  ---
!
            IF( KP.LT.KP_MPI ) THEN
!
!---          Starting and ending global k indices with ghost cells  ---
!
              KGX = K4X
              DO JGX = J1X,J2X
                DO IGX = I1X,I2X
                  N = ND(IGX,JGX,KGX)
                  IF( IXP(N).GT.0 ) THEN
                    NGC = NGC + 1
                    NLRX = (KGX-K3X)*(J4X-J3X+1)*(I4X-I3X+1) +
     &                (JGX-J3X)*(I4X-I3X+1) + (IGX-I3X+1)
                    NLRGCX(NGC) = NLRX
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Create a mapping of the sequential node numbering  --
!
      DO N = 1,NFLD
        NDX(N) = N
      ENDDO
      WRITE(20) ((IXPX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) (((INBS(J,ND_MPI(K,L)),J=1,6),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((NDX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) (NUKFLX(L),L=1,NP_MPI)
      WRITE(20) (NUKFOX(L),L=1,NP_MPI)
      WRITE(20) NUKFGX
      WRITE(20) (NUKTLX(L),L=1,NP_MPI)
      WRITE(20) (NUKTOX(L),L=1,NP_MPI)
      WRITE(20) NUKTGX
      WRITE(20) IFLD,JFLD,KFLD,NFLD,NXP,ICS,IP_MPI,JP_MPI,KP_MPI,NP_MPI
!
!---  Total number of ghost cells per processor  --- 
!
      WRITE(20) (NGHCX(L),L=1,NP_MPI)
!
!---  Number of ghost cells to be sent, indexed by sending direction
!     and sending processor number, where processor number is one 
!     plus the processor rank  --- 
!
      WRITE(20) ((NCSGCX(K,L),K=1,6),L=1,NP_MPI)
!
!---  Process receiving ghost cells to be sent, indexed by sending 
!     direction and sending processor number, where processor number 
!     is one plus the processor rank  --- 
!
      WRITE(20) ((NPSGCX(K,L),K=1,6),L=1,NP_MPI)
!
!---  Create a listing of the local node numbers on processors
!     sending and receiving ghost cells. Ordering is ijk for bottom, 
!     south, west, east, north, top ghost cells on each processor  ---
!
      WRITE(20) (NLSGCX(L),L=1,NGC)
      WRITE(20) (NLRGCX(L),L=1,NGC)
      DEALLOCATE( NDX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NDX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( NUKFLX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NUKFLX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( NUKFOX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NUKFOX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( NUKTLX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NUKTLX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( NUKTOX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NUKTOX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( NGHCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NGHCX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( NCSGCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NCSGCX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( NPSGCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NPSGCX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( NLSGCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NLSGCX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( NLRGCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NLRGCX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Close the grid9.bin file  ---
!
      CLOSE(20)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_GRID_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_PROP_CO2
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
!     Write property data to mech.bin, hydr.bin satu1.bin, satu2.bin,
!     perm.bin, comp.bin, and tabl.bin
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
      USE CONST
      USE BCVP
      USE BCVS
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: CPSX,RHOSX
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CMPX,PORX,TORX,PERMX
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SCHRX
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RPGCX,RPLCX
      INTEGER, DIMENSION(:), ALLOCATABLE :: ITORX,IPRFX,ISMX,ISCHRX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IRPGX,IRPLX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISLTBLX,IRGTBLX,IRLTBLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRITE_PROP_CO2'
!
!---  Allocate temporary memory for mechanical property arrays  --
!
      ALLOCATE( CPSX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: CPSX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( RHOSX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: RHOSX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( CMPX(1:4,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: CMPX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( PORX(1:6,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: PORX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( TORX(1:6,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: TORX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( ITORX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ITORX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Fill temporary mechanical property arrays  --
!
      DO N = 1,NFLD
        IF( IZ(N).LE.0 ) THEN
          RHOSX(N) = 2.65D+3
          CPSX(N) = 9.D+2
          DO M = 1,4
            CMPX(M,N) = 0.D+0
          ENDDO
          DO M = 1,6
            PORX(M,N) = 3.D-1
            TORX(M,N) = 0.D+0
          ENDDO
          ITORX(N) = 0
        ELSE
          RHOSX(N) = RHOS(IZ(N))
          CPSX(N) = CPS(IZ(N))
          DO M = 1,4
            CMPX(M,N) = CMP(M,IZ(N))
          ENDDO
          DO M = 1,6
            PORX(M,N) = POR(M,IZ(N))
            TORX(M,N) = TOR(M,IZ(N))
          ENDDO
          ITORX(N) = ITOR(IZ(N))
        ENDIF
      ENDDO
!
!---  Open a new mech.bin file  --
!
      OPEN( UNIT=20,FILE='mech.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='mech.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write rock/soil array to mech.bin  --
!
      WRITE(20) ((IZ(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Write mechanical property arrays to mech.bin  --
!
      WRITE(20) ((RHOSX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((CPSX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) (((CMPX(M,ND_MPI(K,L)),M=1,4),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      WRITE(20) (((PORX(M,ND_MPI(K,L)),M=1,6),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      WRITE(20) (((TORX(M,ND_MPI(K,L)),M=1,6),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      WRITE(20) ((ITORX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Close the mech.bin file  ---
!
      CLOSE(20)
!
!---  Deallocate temporary memory for mechanical property arrays  --
!
      DEALLOCATE( CPSX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: CPSX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( RHOSX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: RHOSX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( CMPX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: CMPX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( PORX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: PORX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( TORX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: TORX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( ITORX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: ITORX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Allocate temporary memory for hydraulic property arrays  --
!
      ALLOCATE( PERMX(1:9,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: PERMX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IPRFX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IPRFX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Fill temporary hydraulic property arrays  --
!
      DO N = 1,NFLD
        IF( IZ(N).LE.0 ) THEN
          DO M = 1,9
            PERMX(M,N) = 0.D+0
          ENDDO
          IPRFX(N) = 0
        ELSE
          DO M = 1,9
            PERMX(M,N) = PERM(M,IZ(N))
          ENDDO
          IPRFX(N) = IPRF(IZ(N))
        ENDIF
      ENDDO
!
!---  Open a new hydr.bin file  --
!
      OPEN( UNIT=20,FILE='hydr.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='hydr.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write hydraulic property arrays to hydr.bin  --
!
      WRITE(20) (((PERMX(M,ND_MPI(K,L)),M=1,9),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      WRITE(20) ((IPRFX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Close the hydr.bin file  ---
!
      CLOSE(20)
!
!---  Deallocate temporary memory for hydraulic property arrays  --
!
      DEALLOCATE( PERMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: PERMX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IPRFX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IPRFX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Allocate temporary memory for saturation function arrays  --
!
      ALLOCATE( SCHRX(1:LSCHR,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: SCHRX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( ISCHRX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ISCHRX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( ISMX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ISMX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( ISLTBLX(1:4,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ISLTBLX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IRGTBLX(1:4,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IRGTBLX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IRLTBLX(1:4,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IRLTBLX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Fill temporary saturation function arrays  --
!
      DO N = 1,NFLD
        ISCHRX(N) = 0
        ISMX(N) = 0
        DO M = 1,LSCHR
          SCHRX(M,N) = 0.D+0
        ENDDO
        DO M = 1,4
          ISLTBLX(M,N) = 0
          IRGTBLX(M,N) = 0
          IRLTBLX(M,N) = 0
        ENDDO
        IF( IZ(N).LE.0 ) THEN
          ISCHRX(N) = 0
          ISMX(N) = 0
          DO M = 1,LSCHR
            SCHRX(M,N) = 0.D+0
          ENDDO
          DO M = 1,4
            ISLTBLX(M,N) = 0
            IRGTBLX(M,N) = 0
            IRLTBLX(M,N) = 0
          ENDDO
        ELSE
          ISCHRX(N) = ISCHR(IZ(N))
          ISMX(N) = ISM(IZ(N))
          DO M = 1,LSCHR
            SCHRX(M,N) = SCHR(M,IZ(N))
          ENDDO
          DO M = 1,4
            ISLTBLX(M,N) = ISLTBL(M,IZ(N))
            IRGTBLX(M,N) = IRGTBL(M,IZ(N))
            IRLTBLX(M,N) = IRLTBL(M,IZ(N))
          ENDDO
        ENDIF
      ENDDO
!
!---  Open a new satu1.bin file  --
!
      OPEN( UNIT=20,FILE='satu1.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='satu1.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write saturation function arrays to satu1.bin  --
!
      WRITE(20) (((SCHRX(M,ND_MPI(K,L)),M=1,LSCHR),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
!
!---  Close the satu1.bin file  ---
!
      CLOSE(20)
!
!---  Open a new satu2.bin file  --
!
      OPEN( UNIT=20,FILE='satu2.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='satu2.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write saturation function arrays to satu2.bin  --
!
      WRITE(20) ((ISCHRX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((ISMX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) (((ISLTBLX(M,ND_MPI(K,L)),M=1,4),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      WRITE(20) (((IRGTBLX(M,ND_MPI(K,L)),M=1,4),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      WRITE(20) (((IRLTBLX(M,ND_MPI(K,L)),M=1,4),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
!
!---  Close the satu2.bin file  ---
!
      CLOSE(20)
!
!---  Deallocate temporary memory for saturation function arrays  --
!
      DEALLOCATE( SCHRX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: SCHRX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( ISCHRX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: ISCHRX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( ISMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: ISMX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( ISLTBLX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: ISLTBLX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IRGTBLX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IRGTBLX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IRLTBLX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IRLTBLX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Allocate temporary memory for relative permeability
!     function arrays  --
!
      ALLOCATE( RPGCX(1:LRPGC,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: RPGCX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IRPGX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IRPGX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( RPLCX(1:LRPLC,1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: RPLCX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( IRPLX(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IRPLX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Fill temporary relative permeability function arrays  --
!
      DO N = 1,NFLD
        IRPGX(N) = 0
        IRPLX(N) = 0
        DO M = 1,LRPGC
          RPGCX(M,N) = 0.D+0
        ENDDO
        DO M = 1,LRPLC
          RPLCX(M,N) = 0.D+0
        ENDDO
        IF( IZ(N).LE.0 ) THEN
          IRPGX(N) = 0
          IRPLX(N) = 0
          DO M = 1,LRPGC
            RPGCX(M,N) = 0.D+0
          ENDDO
          DO M = 1,LRPLC
            RPLCX(M,N) = 0.D+0
          ENDDO
        ELSE
          IRPGX(N) = IRPG(IZ(N))
          IRPLX(N) = IRPL(IZ(N))
          DO M = 1,LRPGC
            RPGCX(M,N) = RPGC(M,IZ(N))
          ENDDO
          DO M = 1,LRPLC
            RPLCX(M,N) = RPLC(M,IZ(N))
          ENDDO
        ENDIF
      ENDDO
!
!---  Open a new perm.bin file  --
!
      OPEN( UNIT=20,FILE='perm.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='perm.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write relative permeability function arrays to perm.bin  --
!
      WRITE(20) (((RPGCX(M,ND_MPI(K,L)),M=1,LRPGC),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      WRITE(20) ((IRPGX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) (((RPLCX(M,ND_MPI(K,L)),M=1,LRPLC),K=1,NC_MPI(L)),
     &  L=1,NP_MPI)
      WRITE(20) ((IRPLX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Close the perm.bin file  ---
!
      CLOSE(20)
!
!---  Open a new comp.bin file  --
!
      OPEN( UNIT=20,FILE='comp.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='comp.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write initial compressibility state arrays to comp.bin  --
!
      WRITE(20) ((PCMP(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((TCMP(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Close the comp.bin file  ---
!
      CLOSE(20)
!
!---  Open a new tabl.bin file  --
!
      OPEN( UNIT=20,FILE='tabl.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='tabl.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write table lookup arrays to tabl.bin  --
!
      WRITE(20) NTBL
      WRITE(20) (TBLX(L),L=1,LTBL)
      WRITE(20) (TBLY(L),L=1,LTBL)
      WRITE(20) (TBLDDX(L),L=1,LTBL)
      WRITE(20) (TBLDDY(L),L=1,LTBL)
!
!---  Close the tabl.bin file  ---
!
      CLOSE(20)
!
!---  Deallocate temporary memory for relative permeability
!     function arrays  --
!
      DEALLOCATE( RPGCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: RPGCX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IRPGX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IRPGX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( RPLCX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: RPLCX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( IRPLX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IRPLX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_PROP_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_SOLU_CO2( IDX,JDX,KDX )
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
!     Write solution control and output data to solu.bin
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
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
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CNVSF
      REAL*8, DIMENSION(:), ALLOCATABLE :: CNVREF,CNVPLOT
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARX
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSFX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVARX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISFNX,ISFSX
      INTEGER, DIMENSION(1:2,1:LPX_MPI) :: IDX
      INTEGER, DIMENSION(1:2,1:LPY_MPI) :: JDX
      INTEGER, DIMENSION(1:2,1:LPZ_MPI) :: KDX
      INTEGER, DIMENSION(6) :: ISFCX
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: CHREFX
      CHARACTER(64), DIMENSION(:,:), ALLOCATABLE :: CHSFX
      CHARACTER*6 FORM21
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM21 /'(I1)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRITE_SOLU_CO2'
!
!---  Open a new solu.bin file  --
!
      OPEN( UNIT=20,FILE='solu.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='solu.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Load parameters into a temporary integer array and then
!     write the array to the solu.bin file  --
!
      ALLOCATE( IVARX(1:341),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IVARX'
        CALL WRMSGS( INDX )
      ENDIF
      IVARX(1) = L_BH
      IVARX(2) = L_CW
      IVARX(3) = L_DP
      IVARX(4) = L_EC
      IVARX(5) = L_FRC
      IVARX(6) = L_LV
      IVARX(7) = L_LW
      IVARX(8) = L_SFC
      IVARX(9) = LAD
      IVARX(10) = LALC
      IVARX(11) = LAN
      IVARX(12) = LANI
      IVARX(13) = LANW
      IVARX(14) = LATM
      IVARX(15) = LBAL
      IVARX(16) = LBC
      IVARX(17) = LBC_BH
      IVARX(18) = LBC_EC
      IVARX(19) = LBC_FRC
      IVARX(20) = LBC_GM
      IVARX(21) = LBCA
      IVARX(22) = LBCC
      IVARX(23) = LBCG
      IVARX(24) = LBCGC
      IVARX(25) = LBCH
      IVARX(26) = LBCI
      IVARX(27) = LBCIN
      IVARX(28) = LBCIN_GM
      IVARX(29) = LBCL
      IVARX(30) = LBCN
      IVARX(31) = LBCN2
      IVARX(32) = LBCS
      IVARX(33) = LBCT
      IVARX(34) = LBCU
      IVARX(35) = LBCV
      IVARX(36) = LBCV_GM
      IVARX(37) = LBD
      IVARX(38) = LBN_BH
      IVARX(39) = LBN_BHC
      IVARX(40) = LBR
      IVARX(41) = LBTM
      IVARX(42) = LBTM_GM
      IVARX(43) = LC
      IVARX(44) = LCAT
      IVARX(45) = LCDC
      IVARX(46) = LCDP
      IVARX(47) = LCDS
      IVARX(48) = LCH_HT
      IVARX(49) = LCHEM
      IVARX(50) = LCKN
      IVARX(51) = LCMP
      IVARX(52) = LCN
      IVARX(53) = LCN_HT
      IVARX(54) = LCOAX_BH
      IVARX(55) = LCP_HT
      IVARX(56) = LD
      IVARX(57) = LDCO2
      IVARX(58) = LEPD
      IVARX(59) = LEQC
      IVARX(60) = LEQE
      IVARX(61) = LEQK
      IVARX(62) = LESITE
      IVARX(63) = LF_FRC
      IVARX(64) = LF_FRCC
      IVARX(65) = LFC_BH
      IVARX(66) = LFC_FRC
      IVARX(67) = LFD
      IVARX(68) = LFD_DP
      IVARX(69) = LFD_EC
      IVARX(70) = LFDA
      IVARX(71) = LFDC
      IVARX(72) = LFDCR
      IVARX(73) = LFDD
      IVARX(74) = LFDG
      IVARX(75) = LFDGC
      IVARX(76) = LFDH
      IVARX(77) = LFDI
      IVARX(78) = LFDL
      IVARX(79) = LFDM
      IVARX(80) = LFDN
      IVARX(81) = LFDN2
      IVARX(82) = LFDNH
      IVARX(83) = LFDR
      IVARX(84) = LFDRG
      IVARX(85) = LFDRL
      IVARX(86) = LFDRN
      IVARX(87) = LFDS
      IVARX(88) = LFDT
      IVARX(89) = LFEN
      IVARX(90) = LFILES
      IVARX(91) = LFW
      IVARX(92) = LFX
      IVARX(93) = LFX_MPI
      IVARX(94) = LFXY
      IVARX(95) = LFY
      IVARX(96) = LFY_MPI
      IVARX(97) = LFYZ
      IVARX(98) = LFZ
      IVARX(99) = LFZ_MPI
      IVARX(100) = LFZX
      IVARX(101) = LG
      IVARX(102) = LGC
      IVARX(103) = LHBW
      IVARX(104) = LHE_HT
      IVARX(105) = LHF_HT
      IVARX(106) = LHYD
      IVARX(107) = LI_BH
      IVARX(108) = LINC
      IVARX(109) = LINH
      IVARX(110) = LIS
      IVARX(111) = LJA
      IVARX(112) = LJB
      IVARX(113) = LJC
      IVARX(114) = LJC_BH
      IVARX(115) = LJC_GM
      IVARX(116) = LJD
      IVARX(117) = LJE
      IVARX(118) = LJF
      IVARX(119) = LJG
      IVARX(120) = LJG_BCF
      IVARX(121) = LJG_BCM
      IVARX(122) = LJG_BH
      IVARX(123) = LJG_FCB
      IVARX(124) = LJG_FCM
      IVARX(125) = LJG_FRC
      IVARX(126) = LJG_GM
      IVARX(127) = LJG_MCB
      IVARX(128) = LJG_MCF
      IVARX(129) = LJH
      IVARX(130) = LJH_BCF
      IVARX(131) = LJH_BCM
      IVARX(132) = LJH_BH
      IVARX(133) = LJH_FCB
      IVARX(134) = LJH_FCM
      IVARX(135) = LJH_FRC
      IVARX(136) = LJH_GM
      IVARX(137) = LJH_MCB
      IVARX(138) = LJH_MCF
      IVARX(139) = LJI
      IVARX(140) = LJJ
      IVARX(141) = LJK
      IVARX(142) = LJK_BCF
      IVARX(143) = LJK_BCM
      IVARX(144) = LJK_BH
      IVARX(145) = LJK_FCB
      IVARX(146) = LJK_FCM
      IVARX(147) = LJK_FRC
      IVARX(148) = LJK_MCB
      IVARX(149) = LJK_MCF
      IVARX(150) = LJL
      IVARX(151) = LJL_BCF
      IVARX(152) = LJL_BCM
      IVARX(153) = LJL_BH
      IVARX(154) = LJL_FCB
      IVARX(155) = LJL_FCM
      IVARX(156) = LJL_FRC
      IVARX(157) = LJL_MCB
      IVARX(158) = LJL_MCF
      IVARX(159) = LJM
      IVARX(160) = LJN
      IVARX(161) = LJN_BH
      IVARX(162) = LJO
      IVARX(163) = LJO_GM
      IVARX(164) = LL
      IVARX(165) = LM
      IVARX(166) = LMC
      IVARX(167) = LMCG
      IVARX(168) = LMNP
      IVARX(169) = LMPH
      IVARX(170) = LN
      IVARX(171) = LN_BH
      IVARX(172) = LN_BHC
      IVARX(173) = LN_CW
      IVARX(174) = LN_LW
      IVARX(175) = LN2
      IVARX(176) = LNAF
      IVARX(177) = LNC_FRC
      IVARX(178) = LNCF
      IVARX(179) = LNEU
      IVARX(180) = LNGC
      IVARX(181) = LNHC
      IVARX(182) = LNNF
      IVARX(183) = LNNGC
      IVARX(184) = LNOTES
      IVARX(185) = LNW
      IVARX(186) = LNWN
      IVARX(187) = LNWS
      IVARX(188) = LNWT
      IVARX(189) = LNWV
      IVARX(190) = LO_PH
      IVARX(191) = LO_TH
      IVARX(192) = LOBDS
      IVARX(193) = LOBDT
      IVARX(194) = LOUPV
      IVARX(195) = LP_MPI
      IVARX(196) = LP_TA
      IVARX(197) = LPC
      IVARX(198) = LPCF
      IVARX(199) = LPE_HT
      IVARX(200) = LPF_EOR
      IVARX(201) = LPH
      IVARX(202) = LPLANT
      IVARX(203) = LPOLYC
      IVARX(204) = LPOLYN
      IVARX(205) = LPP_HT
      IVARX(206) = LPT
      IVARX(207) = LPTA
      IVARX(208) = LPTM
      IVARX(209) = LPX_MPI
      IVARX(210) = LPY_MPI
      IVARX(211) = LPZ_MPI
      IVARX(212) = LR
      IVARX(213) = LRC
      IVARX(214) = LRCE
      IVARX(215) = LRCG
      IVARX(216) = LRCK
      IVARX(217) = LRCL
      IVARX(218) = LRCN
      IVARX(219) = LRCS
      IVARX(220) = LRCT
      IVARX(221) = LREF
      IVARX(222) = LREK
      IVARX(223) = LREL
      IVARX(224) = LREM
      IVARX(225) = LRFN
      IVARX(226) = LRK
      IVARX(227) = LRPGC
      IVARX(228) = LRPL
      IVARX(229) = LRPLC
      IVARX(230) = LRPNC
      IVARX(231) = LS
      IVARX(232) = LSALC
      IVARX(233) = LSCHR
      IVARX(234) = LSEC
      IVARX(235) = LSEE
      IVARX(236) = LSEK
      IVARX(237) = LSF
      IVARX(238) = LSFCA
      IVARX(239) = LSFCC
      IVARX(240) = LSFCN
      IVARX(241) = LSFCP
      IVARX(242) = LSFCT
      IVARX(243) = LSFDOM
      IVARX(244) = LSFV
      IVARX(245) = LSFVGC
      IVARX(246) = LSOLSR
      IVARX(247) = LSOLU
      IVARX(248) = LSOLU_BH
      IVARX(249) = LSOLU_CW
      IVARX(250) = LSP
      IVARX(251) = LSPBC
      IVARX(252) = LSPC_CW
      IVARX(253) = LSPE
      IVARX(254) = LSPG
      IVARX(255) = LSPILL
      IVARX(256) = LSPK
      IVARX(257) = LSPL
      IVARX(258) = LSPLK
      IVARX(259) = LSPN
      IVARX(260) = LSPR
      IVARX(261) = LSPS
      IVARX(262) = LSPT
      IVARX(263) = LSR
      IVARX(264) = LSR_BH
      IVARX(265) = LSR_FRC
      IVARX(266) = LSRX
      IVARX(267) = LSRY
      IVARX(268) = LSRZ
      IVARX(269) = LSTC
      IVARX(270) = LSTM
      IVARX(271) = LSTM_BH
      IVARX(272) = LSTM_FRC
      IVARX(273) = LSU
      IVARX(274) = LSV
      IVARX(275) = LSW
      IVARX(276) = LSX
      IVARX(277) = LSXC
      IVARX(278) = LSXG
      IVARX(279) = LSXGC
      IVARX(280) = LSXL
      IVARX(281) = LSXLC
      IVARX(282) = LSXN
      IVARX(283) = LSXN2
      IVARX(284) = LSXNC
      IVARX(285) = LSXS
      IVARX(286) = LSXT
      IVARX(287) = LSY
      IVARX(288) = LSYC
      IVARX(289) = LSYG
      IVARX(290) = LSYGC
      IVARX(291) = LSYL
      IVARX(292) = LSYLC
      IVARX(293) = LSYN
      IVARX(294) = LSYN2
      IVARX(295) = LSYNC
      IVARX(296) = LSYS
      IVARX(297) = LSYT
      IVARX(298) = LSZ
      IVARX(299) = LSZC
      IVARX(300) = LSZG
      IVARX(301) = LSZGC
      IVARX(302) = LSZL
      IVARX(303) = LSZLC
      IVARX(304) = LSZN
      IVARX(305) = LSZN2
      IVARX(306) = LSZNC
      IVARX(307) = LSZS
      IVARX(308) = LSZT
      IVARX(309) = LSZW
      IVARX(310) = LT
      IVARX(311) = LT_BH
      IVARX(312) = LT_FRC
      IVARX(313) = LT_FRCC
      IVARX(314) = LT_PH
      IVARX(315) = LT_TA
      IVARX(316) = LT_TH
      IVARX(317) = LTBL
      IVARX(318) = LTC_FRC
      IVARX(319) = LTP_HT
      IVARX(320) = LUGR
      IVARX(321) = LUK
      IVARX(322) = LUK_BH
      IVARX(323) = LUK_CW
      IVARX(324) = LUK_SFC
      IVARX(325) = LUKW
      IVARX(326) = LVIC_FRC
      IVARX(327) = LVPLOT
      IVARX(328) = LVREF
      IVARX(329) = LWELL
      IVARX(330) = LWF_CW
      IVARX(331) = LWF_LW
      IVARX(332) = LWI_CW
      IVARX(333) = LWI_LW
      IVARX(334) = LWN_CW
      IVARX(335) = LWN_LW
      IVARX(336) = LWSI
      IVARX(337) = LWT_CW
      IVARX(338) = LWTI
      IVARX(339) = LWTP_CW
      IVARX(340) = LXP_FRC
      IVARX(341) = LXYZG
      WRITE(20) (IVARX(M),M=1,341)
      DEALLOCATE( IVARX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IVARX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Load time-step variables and PETSc convergence parameters
!     into a temporary real array and then
!     write the array to the solu.bin file  --
!
      ALLOCATE( VARX(1:18),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: VARX'
        CALL WRMSGS( INDX )
      ENDIF
      VARX(1) = TM
      VARX(2) = TMMX
      VARX(3) = TMPR
      VARX(4) = DT 
      VARX(5) = DTI 
      VARX(6) = DTMX 
      VARX(7) = DTMN 
      VARX(8) = DTAF 
      VARX(9) = DTCF 
      VARX(10) = DTO 
      VARX(11) = DTSO
      VARX(12) = RSDMX
      VARX(13) = RLXF
      VARX(14) = CRNTMXC
      VARX(15) = RTOL_PETSC
      VARX(16) = ATOL_PETSC
      VARX(17) = DTOL_PETSC
      VARX(18) = REAL(MAXITS_PETSC)
      WRITE(20) (VARX(M),M=1,18)
      DEALLOCATE( VARX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: VARX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Write execution period values to the solu.bin file  --
!
      WRITE(20) (TMPS(M),M=1,LEPD)
      WRITE(20) (TMPE(M),M=1,LEPD)
      WRITE(20) (TMPD(M),M=1,LEPD)
      WRITE(20) (TMPX(M),M=1,LEPD)
      WRITE(20) (TMPN(M),M=1,LEPD)
      WRITE(20) (TMPA(M),M=1,LEPD)
      WRITE(20) (TMPC(M),M=1,LEPD)
      WRITE(20) (RSDM(M),M=1,LEPD)
      WRITE(20) (RSD(M),M=1,(LUK*(1+LWELL+LSPILL)))
      WRITE(20) (WFMN(M),M=1,20)
!
!---  Load solution control integers into a temporary array and then
!     write the array to the solu.bin file  --
!
      ALLOCATE( IVARX(1:39),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IVARX'
        CALL WRMSGS( INDX )
      ENDIF
      IVARX(1) = IVRSN
      IVARX(2) = ISIC
      IVARX(3) = ICNV
      IVARX(4) = IEO
      IVARX(5) = ILES
      IVARX(6) = IOM
      IVARX(7) = ICODE
      IVARX(8) = IEQT
      IVARX(9) = IEQW
      IVARX(10) = IEQA
      IVARX(11) = IEQN
      IVARX(12) = IEQO
      IVARX(13) = IEQC
      IVARX(14) = IEQS
      IVARX(15) = IEQD
      IVARX(16) = IEQDO
      IVARX(17) = IEQHA
      IVARX(18) = IEQHN
      IVARX(19) = IEQHO
      IVARX(20) = IEQDA
      IVARX(21) = IAQU
      IVARX(22) = IGAS
      IVARX(23) = INAPL
      IVARX(24) = NEPD
      IVARX(25) = MEPD
      IVARX(26) = IEPD
      IVARX(27) = NRIMX
      IVARX(28) = NSTEP
      IVARX(29) = NRST
      IVARX(30) = NITER
      IVARX(31) = NTSR
      IVARX(32) = NGC
      IVARX(33) = MXSTEP
      IVARX(34) = IUNM
      IVARX(35) = IUNKG
      IVARX(36) = IUNS
      IVARX(37) = IUNK
      IVARX(38) = IUNMOL
      IVARX(39) = ISVC
      WRITE(20) (IVARX(M),M=1,39)
      DEALLOCATE( IVARX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IVARX'
        CALL WRMSGS( INDX )
      ENDIF
      WRITE(20) (ISLC(M),M=1,100)
      WRITE(20) (IDMN(M),M=1,20)
!
!---  Write unit conversions for outputs to the solu.bin file  ---
!
      ALLOCATE( CNVPLOT(1:LVPLOT),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: CNVPLOT'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( CNVREF(1:LVREF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: CNVREF'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( CNVSF(1:2,1:LSF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: CNVSF'
        CALL WRMSGS( INDX )
      ENDIF
      CNVTM = 1.D+0
      INDX = 5
      CALL RDUNIT( UNTM,CNVTM,INDX)
      CNVLN = 1.D+0
      INDX = 5
      CALL RDUNIT( UNLN,CNVLN,INDX)
      DO M = 1,LVPLOT
        IF( IPLOT(M).EQ.4 ) THEN
          IF( UNPLOT(4).EQ.'c' ) THEN
            CNVPLOT(M) = -1.D+0
          ELSEIF( UNPLOT(4).EQ.'k' ) THEN
            CNVPLOT(M) = -2.D+0
          ELSEIF( UNPLOT(4).EQ.'f' ) THEN
            CNVPLOT(M) = -3.D+0
          ELSEIF( UNPLOT(4).EQ.'r' ) THEN
            CNVPLOT(M) = -4.D+0
          ENDIF          
        ELSE
          INDX = 5
          CNVPLOT(M) = 1.D+0
          CALL RDUNIT( UNPLOT(IPLOT(M)),CNVPLOT(M),INDX )
        ENDIF
      ENDDO
      DO M = 1,LVREF
        IF( IREF(M).EQ.4 ) THEN
          IF( UNREF(4).EQ.'c' ) THEN
            CNVREF(M) = -1.D+0
          ELSEIF( UNREF(4).EQ.'k' ) THEN
            CNVREF(M) = -2.D+0
          ELSEIF( UNREF(4).EQ.'f' ) THEN
            CNVREF(M) = -3.D+0
          ELSEIF( UNREF(4).EQ.'r' ) THEN
            CNVREF(M) = -4.D+0
          ENDIF          
        ELSE
          INDX = 5
          CNVREF(M) = 1.D+0
          CALL RDUNIT( UNREF(IREF(M)),CNVREF(M),INDX )
        ENDIF
      ENDDO
      DO M = 1,LSF
        DO L = 1,2
          INDX = 5
          CNVSF(L,M) = 1.D+0
          CALL RDUNIT( UNSF(L,M),CNVSF(L,M),INDX )
        ENDDO
      ENDDO
!
!---  Write output unit conversions to the solu.bin file  ---
!
      WRITE(20) CNVTM
      WRITE(20) CNVLN
      WRITE(20) (CNVPLOT(M),M=1,LVPLOT)
      WRITE(20) (CNVREF(M),M=1,LVREF)
      WRITE(20) ((CNVSF(L,M),L=1,2),M=1,LSF)
      DEALLOCATE( CNVPLOT,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: CNVPLOT'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( CNVREF,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: CNVREF'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( CNVSF,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: CNVSF'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Write output unit to the solu.bin file  ---
!
      WRITE(20) UNTM
      WRITE(20) UNLN
      WRITE(20) (UNPLOT(IPLOT(M)),M=1,LVPLOT)
      WRITE(20) (UNREF(IREF(M)),M=1,LVREF)
      WRITE(20) ((UNSF(L,M),L=1,2),M=1,LSF)
!
!---  Write output control variables to the solu.bin file  ---
!
      WRITE(20) (PRTM(M),M=1,LPTM)
!
!---  Write output control integer arrays to the solu.bin file  ---
!
      WRITE(20) (IPLOT(M),M=1,LVPLOT)
      WRITE(20) (IREF(M),M=1,LVREF)
      WRITE(20) (NDREF(M),M=1,LREF)
      WRITE(20) (ISFT(M),M=1,LSF)
      WRITE(20) (ISFF(M),M=1,LSF)
      WRITE(20) (ISFD(M),M=1,LSF)
      WRITE(20) (ISFGP(M),M=1,LSF)
!
!---  Load output control integers into a temporary array and then
!     write the array to the solu.bin file  --
!
      ALLOCATE( IVARX(1:14),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IVARX'
        CALL WRMSGS( INDX )
      ENDIF
      IVARX(1) = NPRTM
      IVARX(2) = NVPLOT
      IVARX(3) = NREF
      IVARX(4) = NVREF
      IVARX(5) = ICNO
      IVARX(6) = ICNS
      IVARX(7) = NSF
      IVARX(8) = NSFGP
      IVARX(9) = IHSF
      IVARX(10) = IFQS
      IVARX(11) = IFQO
      IVARX(12) = ISGNS
      IVARX(13) = ISGNO
      IVARX(14) = ISGNP
      WRITE(20) (IVARX(M),M=1,14)
      DEALLOCATE( IVARX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IVARX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Load output character strings into a temporary array and then
!     write the array to the solu.bin file  --
!
      ALLOCATE( CHREFX(1:LVREF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: CHREFX'
        CALL WRMSGS( INDX )
      ENDIF
      DO M = 1,LVREF
        IRNV = IREF(M)
        IRNV_CW = IREF_CW(M)
        CHREFX(M) = '            '
        CHREFX(M)(1:6) = CHREF(IRNV)
!
!---    Insert well number after the letter 'W'  ---
!
        IF( IRNV_CW.GT.0 ) THEN
          ICX = ICOUNT(IRNV_CW)
          WRITE(FORM21(3:3),'(I1)') ICX
          ICH = INDEX(CHREFX(M)(1:),' ') - 1
          IWX = INDEX(CHREFX(M)(1:),'W') + 1
          CHREFX(M)(IWX+ICX:ICH+ICX) = CHREFX(M)(IWX:ICH)
          WRITE(CHREFX(M)(IWX:IWX+ICX-1),FORM21) IRNV_CW
        ENDIF
      ENDDO
      WRITE(20) (CHREFX(M),M=1,LVREF)
      DEALLOCATE( CHREFX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: CHREFX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Append a ".bin" to the surface flux file names, and then
!     write the surface flux file names to the solu.bin file  ---
!
      DO M = 1,LSF
        NCH = INDEX(FNSF(M)(1:),'  ')
        FNSF(M)(NCH:NCH+4) = '.bin'
      ENDDO
      WRITE(20) (FNSF(M),M=1,LSF)
!
!---  Load output character strings into a temporary array and then
!     write the array to the solu.bin file  --
!
      ALLOCATE( CHSFX(1:2,1:LSF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: CHSFX'
        CALL WRMSGS( INDX )
      ENDIF
      DO NS = 1,NSF
        DO M = 1,2
          CHSFX(M,NS) = CHSF(ISFT(NS),ABS(ISFD(NS)))
          NCH = INDEX( CHSFX(M,NS),'  ' )
          IF( M.EQ.1 ) CHSFX(M,NS)(NCH:NCH+1) = 'R('
          IF( M.EQ.2 ) CHSFX(M,NS)(NCH:NCH+1) = 'I('
          WRITE(CHSFX(M,NS)(NCH+2:NCH+4),'(I3)') NS
          CHSFX(M,NS)(NCH+5:NCH+5) = ')'
        ENDDO
      ENDDO
!
!---  Write the surface flux header character strings
!     to the solu.bin file  ---
!
      WRITE(20) ((CHSFX(L,M),L=1,2),M=1,LSF)
      DEALLOCATE( CHSFX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: CHSFX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Determine the total number of surface-flux nodes  ---
!
      NSFNX = 0
      DO NS = 1,NSF
        IF( ISFD(NS).EQ.4 ) THEN
          NSFDOMX = NSFDOM(NS)
        ELSE
          NSFDOMX = 1
        ENDIF
        DO NC = 1,NSFDOMX
          IF( ISFD(NS).EQ.4 ) THEN
            ISFDX = ISFDOM(4,NC,NS)
            DO I = 1,3
              ISFCX((I-1)*2+1) = ISFDOM(I,NC,NS)
              ISFCX((I-1)*2+2) = ISFDOM(I,NC,NS)
            ENDDO
          ELSE
            ISFDX = ISFD(NS)
            DO I = 1,6
              ISFCX(I) = ABS(ISFC(I,NS))
            ENDDO
          ENDIF
        ENDDO
        NSFNX = NSFNX + (ISFCX(2)-ISFCX(1)+1)*(ISFCX(4)-ISFCX(3)+1)*
     &    (ISFCX(6)-ISFCX(5)+1)
      ENDDO
!
!---  Determine the number of surface-flux surfaces on each
!     processor, not including ghost cells  --
!
      ALLOCATE( NSFX(1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NSFX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( ISFNX(1:NSFNX,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ISFNX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( ISFSX(1:NSFNX,1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ISFSX'
        CALL WRMSGS( INDX )
      ENDIF
      DO NP = 1,NP_MPI
        NSFX(NP) = 0
        DO MP = 1,NSFNX
          ISFNX(MP,NP) = 0
          ISFSX(MP,NP) = 0
        ENDDO
      ENDDO
      DO KP = 1,KP_MPI
        DO JP = 1,JP_MPI
          DO IP = 1,IP_MPI
            KFLDX = KD_MPI(2,KP)-KD_MPI(1,KP)+1
            JFLDX = JD_MPI(2,JP)-JD_MPI(1,JP)+1
            IFLDX = ID_MPI(2,IP)-ID_MPI(1,IP)+1
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            DO K = KDX(1,KP),KDX(2,KP)
              DO J = JDX(1,JP),JDX(2,JP)
                DO I = IDX(1,IP),IDX(2,IP)
                  N = ND(I,J,K)
                  DO NS = 1,NSF
                    IF( ISFD(NS).EQ.4 ) THEN
                      NSFDOMX = NSFDOM(NS)
                    ELSE
                      NSFDOMX = 1
                    ENDIF
                    DO NC = 1,NSFDOMX
                      IF( ISFD(NS).EQ.4 ) THEN
                        ISFDX = ISFDOM(4,NC,NS)
                        DO L = 1,3
                          ISFCX((L-1)*2+1) = ISFDOM(L,NC,NS)
                          ISFCX((L-1)*2+2) = ISFDOM(L,NC,NS)
                        ENDDO
                      ELSE
                        ISFDX = ISFD(NS)
                        DO L = 1,6
                          ISFCX(L) = ABS(ISFC(L,NS))
                        ENDDO
                      ENDIF
                      DO IX = ISFCX(1),ISFCX(2)
                      DO JX = ISFCX(3),ISFCX(4)
                      DO KX = ISFCX(5),ISFCX(6)
                        IF( ND(IX,JX,KX).EQ.N ) THEN
!
!---                      Local count of surface nodes  --
!
                          NSFX(NP) = NSFX(NP) + 1
!
!---                      Local node number  --
!
                          ISFNX(NSFX(NP),NP) = 
     &                      (K-KD_MPI(1,KP))*IFLDX*JFLDX +
     &                      (J-JD_MPI(1,JP))*IFLDX + 
     &                      (I-ID_MPI(1,IP)) + 1
!
!---                      Surface index  --
!
                          ISFSX(NSFX(NP),NP) = NS
                        ENDIF
                      ENDDO
                      ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      WRITE(20) (NSFX(NP),NP=1,NP_MPI)
      WRITE(20) ((ISFNX(L,NP),L=1,NSFX(NP)),NP=1,NP_MPI)
      WRITE(20) ((ISFSX(L,NP),L=1,NSFX(NP)),NP=1,NP_MPI)
      DEALLOCATE( NSFX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NSFX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( ISFNX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: ISFNX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( ISFSX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Dellocation Error: ISFSX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Close the solu.bin file  ---
!
      CLOSE(20)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_SOLU_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_SORC_CO2
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
!     Write source data to bcsc.bin
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
      USE CONST
      USE BCVP
      USE BCVS
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISRNX,ISRMX,ISRTX,ISRINX
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSRX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRITE_SORC_CO2'
!
!---  Open a new sorc.bin file  --
!
      OPEN( UNIT=20,FILE='sorc.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='sorc.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Determine the number of source nodes on each
!     processor, not including ghost cells  ---
!
      ALLOCATE( NSRX(1:NP_MPI),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NSRX'
        CALL WRMSGS( INDX )
      ENDIF
      DO NP = 1,NP_MPI
        NSRX(NP) = 0
      ENDDO
      DO NS = 1,NSR
        DO ISDX = ISRDM(1,NS),ISRDM(2,NS)
          DO JSDX = ISRDM(3,NS),ISRDM(4,NS)
            DO KSDX = ISRDM(5,NS),ISRDM(6,NS)
              NSDX = ND(ISDX,JSDX,KSDX)
              NP = NTP_MPI(NSDX)
              NSRX(NP) = NSRX(NP) + 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      NSRSX = 0
      DO NP = 1,NP_MPI
        NSRSX = NSRSX + NSRX(NP)
      ENDDO
      WRITE(20) (NSRX(NP),NP=1,NP_MPI)
!
!---  Write source input variables to bcsc.bin  --
!
      LX = 8+LSOLU+LSPT+LNGC
      WRITE(20) (((SRC(L,M,N),L=1,LX),M=1,LSTM),N=1,LSR)
!
!---  Allocate processor local temporary memory for source indices  --
!
      ALLOCATE( ISRNX(1:NSRSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ISRNX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( ISRMX(1:NSRSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ISRMX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( ISRINX(1:NSRSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ISRINX'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( ISRTX(1:NSRSX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: ISRTX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Create a processor-local map of source time points and types
!     not including ghost cells  ---
!
      NC = 0
      DO KP = 1,KP_MPI
        DO JP = 1,JP_MPI
          DO IP = 1,IP_MPI
            KFLDX = KD_MPI(2,KP)-KD_MPI(1,KP)+1
            JFLDX = JD_MPI(2,JP)-JD_MPI(1,JP)+1
            IFLDX = ID_MPI(2,IP)-ID_MPI(1,IP)+1
            NP = (KP-1)*IP_MPI*JP_MPI + (JP-1)*IP_MPI + IP
            DO NS = 1,NSR
              DO ISDX = ISRDM(1,NS),ISRDM(2,NS)
                DO JSDX = ISRDM(3,NS),ISRDM(4,NS)
                  DO KSDX = ISRDM(5,NS),ISRDM(6,NS)
                    NSDX = ND(ISDX,JSDX,KSDX)
                    IF( NP.EQ.NTP_MPI(NSDX) ) THEN
                      NC = NC + 1
                      ISRNX(NC) = (KSDX-KD_MPI(1,KP))*IFLDX*JFLDX +
     &                  (JSDX-JD_MPI(1,JP))*IFLDX + 
     &                  (ISDX-ID_MPI(1,IP)) + 1
                      ISRTX(NC) = ISRT(NS)
                      ISRMX(NC) = ISRM(NS)
                      ISRINX(NC) = NS
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      WRITE(20) (ISRNX(M),M=1,NSRSX)
      WRITE(20) (ISRMX(M),M=1,NSRSX)
      WRITE(20) (ISRINX(M),M=1,NSRSX)
      WRITE(20) (ISRTX(M),M=1,NSRSX)
!
!---  Deallocate processor local temporary memory for source indices  --
!
      DEALLOCATE( ISRNX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: ISRNX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( ISRMX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: ISRMX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( ISRINX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: ISRINX'
        CALL WRMSGS( INDX )
      ENDIF
      DEALLOCATE( ISRTX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: ISRTX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Deallocate temporary source memory  --
!
      DEALLOCATE( NSRX,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NSRX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Close the sorc.bin file  ---
!
      CLOSE(20)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_SORC_CO2 group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_STATE_CO2
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
!     Write state condition data to state.bin
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
      USE CONST
      USE BCVP
      USE BCVS
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/WRITE_STATE_CO2'
!
!---  Open a new state1.bin file  --
!
      OPEN( UNIT=20,FILE='state1.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='state1.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write state condition arrays to state1.bin  --
!
      WRITE(20) ((T(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((PL(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((PG(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((SG(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((SL(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((SGT(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((ASLMIN(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((TMS(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Close the state1.bin file  ---
!
      CLOSE(20)
!
!---  Open a new state2.bin file  --
!
      OPEN( UNIT=20,FILE='state2.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='state2.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write state condition arrays to state2.bin  --
!
      WRITE(20) ((YLS(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((PVA(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((PVW(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((XLA(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((PORD(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((PORT(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((NPHAZ(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Write Eclipse gas saturation and pressure arrays to state2.bin  ---
!
      WRITE(20) ((SI(1,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
      WRITE(20) ((SI(2,ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
!
!---  Close the state2.bin file  ---
!
      CLOSE(20)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_STATE_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRITE_TPORT_CO2
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
!     Write solute/species transport data to wltp.bin
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE TABL
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE NCG_PT
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE COUP_WELL
      USE CONST
      USE BCVP
      USE BCVS
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: PCSLX
      REAL*8, DIMENSION(:), ALLOCATABLE :: DISPLX,DISPTX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRITE_TPORT_CO2'
!
!---  Open a new tpor.bin file  --
!
      OPEN( UNIT=20,FILE='tpor.bin',STATUS='UNKNOWN',
     &  FORM='UNFORMATTED' )
      CLOSE(UNIT=20,STATUS='DELETE')
      OPEN(UNIT=20, FILE='tpor.bin', STATUS='NEW', 
     &  FORM='UNFORMATTED')
!
!---  Write solute parameters, indices, and variables
!     if solutes are modeled  ---
!
      IF( IEQC.NE.0 ) THEN
!
!---    Write number of solutes  ---
!
        WRITE(20) NSOLU
!
!---    Write solute names  ---
!
        WRITE(20) (SOLUT(NSL),NSL=1,NSOLU)
!
!---    Write initial solute concentrations  ---
!
        DO NSL = 1,NSOLU
          WRITE(20) ((C(ND_MPI(K,L),NSL),K=1,NC_MPI(L)),L=1,NP_MPI)
        ENDDO
!
!---    Write solute aqueous diffusion coefficient  ---
!
        WRITE(20) (SMDL(NSL),NSL=1,NSOLU)
!
!---    Write solute gas diffusion coefficient  ---
!
        WRITE(20) (SMDG(NSL),NSL=1,NSOLU)
!
!---    Write index for solute gas-aqueous partition coefficient  ---
!
        WRITE(20) (IPCGL(NSL),NSL=1,NSOLU)
!
!---    Write solute gas-aqueous partition coefficient parameters  ---
!
        WRITE(20) ((PCGL(M,NSL),M=1,5),NSL=1,NSOLU)
!
!---    Write solute half-life, s  ---
!
        WRITE(20) (HLF(NSL),NSL=1,NSOLU)
!
!---    Write solute chain-decay fraction  ---
!
        WRITE(20) ((CHDF(MSL,NSL),MSL=1,NSOLU),NSL=1,NSOLU)
!
!---    Write number of Bateman chain decay series  ---
!
        WRITE(20) NBCDS
!
!---    Write Bateman chain decay series index  ---
!
        WRITE(20) (IBCDS(L),L=1,LCDC+LSOLU)
!
!---    Write Bateman number of solutes in chain decay path  ---
!
        WRITE(20) (NBCDP(L),L=1,LCDC)
!
!---    Write Bateman chain decay path indices  ---
!
        WRITE(20) (((IBCDP(K,L,M),K=1,LCDS),L=1,LCDP),M=1,LCDC)
!
!---    Write Courant number calculation index  ---
!
        WRITE(20) ICRNT
!
!---    Write maximum Courant number  ---
!
        WRITE(20) CRNTMXT
!
!---    Allocate temporary memory for solid-aqueous partition array  --
!
        ALLOCATE( PCSLX(1:LFD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: PCSLX'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Write solute solid-aqueous partition coefficient  ---
!
        DO NSL = 1,NSOLU
          DO M = 1,5
!
!---        Fill temporary solid-aqueous partition array  --
!
            DO N = 1,NFLD
              IF( IZ(N).LE.0 ) THEN
                IF( M.EQ.1 ) THEN
                  PCSLX(N) = 1.D-20
                ELSE
                  PCSLX(N) = 0.D+0
                ENDIF
              ELSE
                PCSLX(N) = PCSL(M,IZ(N),NSL)
              ENDIF
            ENDDO
            WRITE(20) ((PCSLX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
          ENDDO
        ENDDO
!
!---    Deallocate temporary memory for solid-aqueous partition array  --
!
        DEALLOCATE( PCSLX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: PCSLX'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Allocate temporary memory for longitudinal and transverse
!       dispersivity arrays  --
!
        ALLOCATE( DISPLX(1:LFD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: DISPLX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( DISPTX(1:LFD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: DISPTX'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Fill temporary dispersivity arrays  --
!
        DO N = 1,NFLD
          IF( IZ(N).LE.0 ) THEN
            DISPLX(N) = 0.D+0
            DISPTX(N) = 0.D+0
          ELSE
            DISPLX(N) = DISPL(IZ(N))
            DISPTX(N) = DISPT(IZ(N))
          ENDIF
        ENDDO
!
!---    Write dispersivity arrays  ---
!
        WRITE(20) ((DISPLX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
        WRITE(20) ((DISPTX(ND_MPI(K,L)),K=1,NC_MPI(L)),L=1,NP_MPI)
        IDISP = 1
!
!---    Deallocate temporary memory for dispersivity arrays  --
!
        DEALLOCATE( DISPLX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: DISPLX'
          CALL WRMSGS( INDX )
        ENDIF
        DEALLOCATE( DISPTX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: DISPTX'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Write initial condition type array for both solutes
!     and reactive species, if solutes or reactive species
!     are modeled  ---
!
      IF( NSOLU.GT.0 .OR. ISLC(40).NE.0 ) THEN
        DO M = 1,LSOLU+LSPT
          WRITE(20) ((ICT(ND_MPI(K,L),M),K=1,NC_MPI(L)),L=1,NP_MPI)
        ENDDO
      ENDIF
!
!---  Close the tpor.bin file  ---
!
      CLOSE(20)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRITE_TPORT_CO2 group
!
      RETURN
      END



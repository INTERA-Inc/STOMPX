!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXN( NSL )
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
!     Compute solute transport flux NAPL-phase, excluding boundaries.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July, 1993.
!     Last Modified by MD White, Battelle, PNL, October 13, 1995.
!     sfxn.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE NAPL
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
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
      SUB_LOG(ISUB_LOG) = '/SFXN'
      M = 1
!
!---  X-direction solute flux NAPL-phase, excluding boundaries
!
      DO N = 1,NFLD
        I = ID(N)
        IF( I.EQ.1 ) CYCLE
        J = JD(N)
        K = KD(N)
        NW = N-1
        IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &    INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
        NPX = NSX(N)
        TCOR = (T(2,N)+TABS)/TSPRF
        SDFNP = SMDN(NSL)*TCOR*(VISRN/VISN(2,N))
        VMCP = SN(2,N)*PORD(2,N)
        DFP = TORN(2,N)*VMCP*SDFNP
        FCNP = YN(N,NSL)/(VMCP+SMALL)
        TCOR = (T(2,NW)+TABS)/TSPRF
        SDFNW = SMDN(NSL)*TCOR*(VISRN/VISN(2,NW))
        VMCW = SN(2,NW)*PORD(2,NW)
        DFW = TORN(2,NW)*VMCW*SDFNW
        FCNW = YN(NW,NSL)/(VMCW+SMALL)
        IF( IDISP .EQ. 1 ) THEN
          CALL ADVW( PORD,SN,UN,VN,WN,UNWX,VNWX,WNWX,N,M )
          UNWSQ = UNWX*UNWX
          VNWSQ = VNWX*VNWX
          WNWSQ = WNWX*WNWX
          ZVW = SQRT(UNWSQ+VNWSQ+WNWSQ)
          INDX = 17
          DPLW = DIFMN(DISPL(IZ(NW)),DISPL(IZ(N)),DXGF(NW),DXGF(N),
     &      UNWX,INDX)
          DPTW = DIFMN(DISPT(IZ(NW)),DISPT(IZ(N)),DXGF(NW),DXGF(N),
     &      UNWX,INDX)
          DPW = (DPLW*UNWSQ + DPTW*(VNWSQ+WNWSQ))/(ZVW+SMALL)
        ELSE
          DPW = 0.D+0
        ENDIF
        INDX = 16
        DFW = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),UN(1,NPX),INDX)
        DDW = (DFW+DPW)/DXGP(NPX)
        IF( ISLC(1).GE.1 )  THEN
          UC(NPX,NSL) = UC(NPX,NSL) + DDW*(C(NW,NSL)*FCNW-C(N,NSL)*FCNP)
        ELSE
          AN = MAX( UN(1,NPX),ZERO ) +
     &      DDW*MAX( (ONE-(TENTH*ABS(UN(1,NPX))/(DDW+SMALL)))**5,ZERO )
          ANP = MAX( -UN(1,NPX),ZERO ) +
     &      DDW*MAX( (ONE-(TENTH*ABS(UN(1,NPX))/(DDW+SMALL)))**5,ZERO )
         UC(NPX,NSL) = UC(NPX,NSL)+(C(NW,NSL)*AN*FCNW-C(N,NSL)*ANP*FCNP)
        ENDIF
      ENDDO
!
!---  Y-direction solute flux NAPL-phase, excluding boundaries
!
      DO N = 1,NFLD
        J = JD(N)
        IF( J.EQ.1 ) CYCLE
        I = ID(N)
        K = KD(N)
        NS = N-IFLD
        IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &    INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
        NPY = NSY(N)
        TCOR = (T(2,N)+TABS)/TSPRF
        SDFNP = SMDN(NSL)*TCOR*(VISRN/VISN(2,N))
        VMCP = SN(2,N)*PORD(2,N)
        DFP = TORN(2,N)*VMCP*SDFNP
        FCNP = YN(N,NSL)/(VMCP+SMALL)
        TCOR = (T(2,NS)+TABS)/TSPRF
        SDFNS = SMDN(NSL)*TCOR*(VISRN/VISN(2,NS))
        VMCS = SN(2,NS)*PORD(2,NS)
        DFS = TORN(2,NS)*VMCS*SDFNS
        FCNS = YN(NS,NSL)/(VMCS+SMALL)
        IF( IDISP .EQ. 1 ) THEN
          CALL ADVS( PORD,SN,UN,VN,WN,UNSX,VNSX,WNSX,N,M )
          UNSSQ = UNSX*UNSX
          VNSSQ = VNSX*VNSX
          WNSSQ = WNSX*WNSX
          ZVS = SQRT(UNSSQ+VNSSQ+WNSSQ)
          INDX = 17
          DPLS = DIFMN(DISPL(IZ(NS)),DISPL(IZ(N)),DYGF(NS),DYGF(N),
     &      VNSX,INDX)
          DPTS = DIFMN(DISPT(IZ(NS)),DISPT(IZ(N)),DYGF(NS),DYGF(N),
     &      VNSX,INDX)
          DPS = (DPLS*VNSSQ + DPTS*(UNSSQ+WNSSQ))/(ZVS+SMALL)
        ELSE
          DPS = 0.D+0
        ENDIF
        INDX = 16
        DFS = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),VN(1,NPY),INDX)
        DDS = (DFS+DPS)/DYGP(NPY)/RP(I)
        IF( ISLC(1).GE.1 )  THEN
          VC(NPY,NSL) = VC(NPY,NSL) + DDS*(C(NS,NSL)*FCNS-C(N,NSL)*FCNP)
        ELSE
          AN = MAX( VN(1,NPY),ZERO ) +
     &      DDS*MAX( (ONE-(TENTH*ABS(VN(1,NPY))/(DDS+SMALL)))**5,ZERO )
          ANP = MAX( -VN(1,NPY),ZERO ) +
     &      DDS*MAX( (ONE-(TENTH*ABS(VN(1,NPY))/(DDS+SMALL)))**5,ZERO )
         VC(NPY,NSL) = VC(NPY,NSL)+(C(NS,NSL)*AN*FCNS-C(N,NSL)*ANP*FCNP)
        ENDIF
      ENDDO
!
!---  Z-direction solute flux NAPL-phase, excludinN boundaries
!
      DO N = 1,NFLD
        K = KD(N)
        IF( K.EQ.1 ) CYCLE
        J = JD(N)
        I = ID(N)
        NB = N-IJFLD
        IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &    INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
        NPZ = NSZ(N)
        TCOR = (T(2,N)+TABS)/TSPRF
        SDFNP = SMDN(NSL)*TCOR*(VISRN/VISN(2,N))
        VMCP = SN(2,N)*PORD(2,N)
        DFP = TORN(2,N)*VMCP*SDFNP
        FCNP = YN(N,NSL)/(VMCP+SMALL)
        TCOR = (T(2,NB)+TABS)/TSPRF
        SDFNB = SMDN(NSL)*TCOR*(VISRN/VISN(2,NB))
        VMCB = SN(2,NB)*PORD(2,NB)
        DFB = TORN(2,NB)*VMCB*SDFNB
        FCNB = YN(NB,NSL)/(VMCB+SMALL)
        IF( IDISP .EQ. 1 ) THEN
          CALL ADVB( PORD,SN,UN,VN,WN,UNBX,VNBX,WNBX,N,M )
          UNBSQ = UNBX*UNBX
          VNBSQ = VNBX*VNBX
          WNBSQ = WNBX*WNBX
          ZVB = SQRT(UNBSQ+VNBSQ+WNBSQ)
          INDX = 17
          DPLB = DIFMN(DISPL(IZ(NB)),DISPL(IZ(N)),DZGF(NB),DZGF(N),
     &      WNBX,INDX)
          DPTB = DIFMN(DISPT(IZ(NB)),DISPT(IZ(N)),DZGF(NB),DZGF(N),
     &      WNBX,INDX)
          DPB = (DPLB*WNBSQ + DPTB*(VNBSQ+UNBSQ))/(ZVB+SMALL)
        ELSE
          DPB = 0.D+0
        ENDIF
        INDX = 16
        DFB = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX)
        DDB = (DFB+DPB)/DZGP(NPZ)
        IF( ISLC(1).GE.1 ) THEN
          WC(NPZ,NSL) = WC(NPZ,NSL) + DDB*(C(NB,NSL)*FCNB-C(N,NSL)*FCNP)
        ELSE
          AN = MAX( WN(1,NPZ),ZERO ) +
     &      DDB*MAX( (ONE-(TENTH*ABS(WN(1,NPZ))/(DDB+SMALL)))**5,ZERO )
          ANP = MAX( -WN(1,NPZ),ZERO ) +
     &      DDB*MAX( (ONE-(TENTH*ABS(WN(1,NPZ))/(DDB+SMALL)))**5,ZERO )
         WC(NPZ,NSL) = WC(NPZ,NSL)+(C(NB,NSL)*AN*FCNB-C(N,NSL)*ANP*FCNP)
        ENDIF
      ENDDO
!
!---  End of SFXN group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END


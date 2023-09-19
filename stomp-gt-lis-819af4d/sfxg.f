!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXG( NSL )
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
!     Compute solute transport flux gas-phase, excluding boundaries.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July, 1993.
!     Last Modified by MD White, Battelle, PNL, October 13, 1995.
!     sfxg.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE HYST
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/SFXG'
      M = 1
!
!---  X-direction solute flux gas-phase, excluding boundaries
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
        PCOR = (PG(2,N)+PATM)/PATM
        SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
        VMCP = SG(2,N)*PORD(2,N)
        DFP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
        FCGP = YG(N,NSL)/(VMCP+SMALL)
        TCOR = (T(2,NW)+TABS)/TSPRF
        PCOR = (PG(2,NW)+PATM)/PATM
        SDFGW = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
        VMCW = SG(2,NW)*PORD(2,NW)
        DFW = TORG(2,NW)*(SG(2,NW)-SGT(2,NW))*PORD(2,NW)*SDFGW
        FCGW = YG(NW,NSL)/(VMCW+SMALL)
        DCG = C(NW,NSL)*FCGW - C(N,NSL)*FCGP
        IF( IDISP .EQ. 1 ) THEN
          CALL ADVW( PORD,SG,UG,VG,WG,UGX,VGX,WGX,N,M )
          UGWSQ = UGX*UGX
          VGWSQ = VGX*VGX
          WGWSQ = WGX*WGX
          ZVW = SQRT(UGWSQ+VGWSQ+WGWSQ)
          INDX = 17
          DPLW = DIFMN(DISPL(IZ(NW)),DISPL(IZ(N)),DXGF(NW),DXGF(N),
     &      UGX,INDX)
          DPTW = DIFMN(DISPT(IZ(NW)),DISPT(IZ(N)),DXGF(NW),DXGF(N),
     &      UGX,INDX)
          DPW = (DPLW*UGWSQ + DPTW*(VGWSQ+WGWSQ))/(ZVW+SMALL)
        ELSE
          DPW = 0.D+0
        ENDIF
        INDX = 16
        DFW = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),UG(1,NPX),INDX)
        DDW = (DFW+DPW)/DXGP(NPX)
        IF( ISLC(1).GE.1 )  THEN
          UC(NPX,NSL) = UC(NPX,NSL) + DDW*(C(NW,NSL)*FCGW-C(N,NSL)*FCGP)
        ELSE
          AG = MAX( UG(1,NPX),ZERO ) +
     &      DDW*MAX( (ONE-(TENTH*ABS(UG(1,NPX))/(DDW+SMALL)))**5,ZERO )
          AGP = MAX( -UG(1,NPX),ZERO ) +
     &      DDW*MAX( (ONE-(TENTH*ABS(UG(1,NPX))/(DDW+SMALL)))**5,ZERO )
         UC(NPX,NSL) = UC(NPX,NSL)+(C(NW,NSL)*AG*FCGW-C(N,NSL)*AGP*FCGP)
        ENDIF
      ENDDO
!
!---  Y-direction solute flux gas-phase, excluding boundaries
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
        PCOR = (PG(2,N)+PATM)/PATM
        SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
        VMCP = SG(2,N)*PORD(2,N)
        DFP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
        FCGP = YG(N,NSL)/(VMCP+SMALL)
        TCOR = (T(2,NS)+TABS)/TSPRF
        PCOR = (PG(2,NS)+PATM)/PATM
        SDFGS = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
        VMCS = SG(2,NS)*PORD(2,NS)
        DFS = TORG(2,NS)*(SG(2,NS)-SGT(2,NS))*PORD(2,NS)*SDFGS
        FCGS = YG(NS,NSL)/(VMCS+SMALL)
        DCG = C(NS,NSL)*FCGS - C(N,NSL)*FCGP
        IF( IDISP .EQ. 1 ) THEN
          CALL ADVS( PORD,SG,UG,VG,WG,UGSX,VGSX,WGSX,N,M )
          UGSSQ = UGSX*UGSX
          VGSSQ = VGSX*VGSX
          WGSSQ = WGSX*WGSX
          ZVS = SQRT(UGSSQ+VGSSQ+WGSSQ)
          INDX = 17
          DPLS = DIFMN(DISPL(IZ(NS)),DISPL(IZ(N)),DYGF(NS),DYGF(N),
     &      VGSX,INDX)
          DPTS = DIFMN(DISPT(IZ(NS)),DISPT(IZ(N)),DYGF(NS),DYGF(N),
     &      VGSX,INDX)
          DPS = (DPLS*VGSSQ + DPTS*(UGSSQ+WGSSQ))/(ZVS+SMALL)
        ELSE
          DPS = 0.D+0
        ENDIF
        INDX = 16
        DFS = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),VG(1,NPY),INDX)
        DDS = (DFS+DPS)/(DYGP(NPY)*RP(I))
        IF( ISLC(1).GE.1 )  THEN
          VC(NPY,NSL) = VC(NPY,NSL) + DDS*(C(NS,NSL)*FCGS-C(N,NSL)*FCGP)
        ELSE
          AG = MAX( VG(1,NPY),ZERO ) +
     &      DDS*MAX( (ONE-(TENTH*ABS(VG(1,NPY))/(DDS+SMALL)))**5,ZERO )
          AGP = MAX( -VG(1,NPY),ZERO ) +
     &      DDS*MAX( (ONE-(TENTH*ABS(VG(1,NPY))/(DDS+SMALL)))**5,ZERO )
         VC(NPY,NSL) = VC(NPY,NSL)+(C(NS,NSL)*AG*FCGS-C(N,NSL)*AGP*FCGP)
        ENDIF
      ENDDO
!
!---  Z-direction solute flux gas-phase, excluding boundaries
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
        PCOR = (PG(2,N)+PATM)/PATM
        SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
        VMCP = SG(2,N)*PORD(2,N)
        DFP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
        FCGP = YG(N,NSL)/(VMCP+SMALL)
        TCOR = (T(2,NB)+TABS)/TSPRF
        PCOR = (PG(2,NB)+PATM)/PATM
        SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
        VMCB = SG(2,NB)*PORD(2,NB)
        DFB = TORG(2,NB)*(SG(2,NB)-SGT(2,NB))*PORD(2,NB)*SDFGB
        FCGB = YG(NB,NSL)/(VMCB+SMALL)
        DCG = C(NB,NSL)*FCGB - C(N,NSL)*FCGP
        IF( IDISP .EQ. 1 ) THEN
          CALL ADVB( PORD,SG,UG,VG,WG,UGBX,VGBX,WGBX,N,M )
          UGBSQ = UGBX*UGBX
          VGBSQ = VGBX*VGBX
          WGBSQ = WGBX*WGBX
          ZVB = SQRT(UGBSQ+VGBSQ+WGBSQ)
          INDX = 17
          DPLB = DIFMN(DISPL(IZ(NB)),DISPL(IZ(N)),DZGF(NB),DZGF(N),
     &      WGBX,INDX)
          DPTB = DIFMN(DISPT(IZ(NB)),DISPT(IZ(N)),DZGF(NB),DZGF(N),
     &      WGBX,INDX)
          DPB = (DPLB*WGBSQ + DPTB*(VGBSQ+UGBSQ))/(ZVB+SMALL)
        ELSE
          DPB = 0.D+0
        ENDIF
        INDX = 16
        DFB = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX)
        DDB = (DFB+DPB)/DZGP(NPZ)
        IF( ISLC(1).GE.1 ) THEN
          WC(NPZ,NSL) = WC(NPZ,NSL) + DDB*(C(NB,NSL)*FCGB-C(N,NSL)*FCGP)
        ELSE
          AG = MAX( WG(1,NPZ),ZERO ) +
     &      DDB*MAX( (ONE-(TENTH*ABS(WG(1,NPZ))/(DDB+SMALL)))**5,ZERO )
          AGP = MAX( -WG(1,NPZ),ZERO ) +
     &      DDB*MAX( (ONE-(TENTH*ABS(WG(1,NPZ))/(DDB+SMALL)))**5,ZERO )
         WC(NPZ,NSL) = WC(NPZ,NSL)+(C(NB,NSL)*AG*FCGB-C(N,NSL)*AGP*FCGP)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFXG group
!
      RETURN
      END



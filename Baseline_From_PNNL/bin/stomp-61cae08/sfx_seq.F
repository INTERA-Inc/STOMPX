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
      DO 100 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 100
        DO NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
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
     &        UGX,INDX)
            DPTW = DIFMN(DISPT(IZ(NW)),DISPT(IZ(N)),DXGF(NW),DXGF(N),
     &        UGX,INDX)
            DPW = (DPLW*UGWSQ + DPTW*(VGWSQ+WGWSQ))/(ZVW+SMALL)
          ELSE
            DPW = 0.D+0
          ENDIF
          INDX = 16
          DFW = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),UG(1,NPX),INDX)
          DDW = (DFW+DPW)/DXGP(NPX)
          AG = MAX( UG(1,NPX),ZERO ) +
     &      DDW*MAX( (ONE-(TENTH*ABS(UG(1,NPX))/(DDW+SMALL)))**5,ZERO )
          AGP = MAX( -UG(1,NPX),ZERO ) +
     &      DDW*MAX( (ONE-(TENTH*ABS(UG(1,NPX))/(DDW+SMALL)))**5,ZERO )
         UC(NPX,NSL) = UC(NPX,NSL)+(C(NW,NSL)*AG*FCGW-C(N,NSL)*AGP*FCGP)
        ENDDO
  100 CONTINUE
!
!---  Y-direction solute flux gas-phase, excluding boundaries
!
      DO 200 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
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
     &        VGSX,INDX)
            DPTS = DIFMN(DISPT(IZ(NS)),DISPT(IZ(N)),DYGF(NS),DYGF(N),
     &        VGSX,INDX)
            DPS = (DPLS*VGSSQ + DPTS*(UGSSQ+WGSSQ))/(ZVS+SMALL)
          ELSE
            DPS = 0.D+0
          ENDIF
          INDX = 16
          DFS = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),VG(1,NPY),INDX)
          DDS = (DFS+DPS)/(DYGP(NPY)*RP(I))
          AG = MAX( VG(1,NPY),ZERO ) +
     &      DDS*MAX( (ONE-(TENTH*ABS(VG(1,NPY))/(DDS+SMALL)))**5,ZERO )
          AGP = MAX( -VG(1,NPY),ZERO ) +
     &      DDS*MAX( (ONE-(TENTH*ABS(VG(1,NPY))/(DDS+SMALL)))**5,ZERO )
         VC(NPY,NSL) = VC(NPY,NSL)+(C(NS,NSL)*AG*FCGS-C(N,NSL)*AGP*FCGP)
        ENDDO
  200 CONTINUE
!
!---  Z-direction solute flux gas-phase, excluding boundaries
!
      DO 300 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 300
        DO NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
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
     &        WGBX,INDX)
            DPTB = DIFMN(DISPT(IZ(NB)),DISPT(IZ(N)),DZGF(NB),DZGF(N),
     &        WGBX,INDX)
            DPB = (DPLB*WGBSQ + DPTB*(VGBSQ+UGBSQ))/(ZVB+SMALL)
          ELSE
            DPB = 0.D+0
          ENDIF
          INDX = 16
          DFB = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX)
          DDB = (DFB+DPB)/DZGP(NPZ)
          AG = MAX( WG(1,NPZ),ZERO ) +
     &      DDB*MAX( (ONE-(TENTH*ABS(WG(1,NPZ))/(DDB+SMALL)))**5,ZERO )
          AGP = MAX( -WG(1,NPZ),ZERO ) +
     &      DDB*MAX( (ONE-(TENTH*ABS(WG(1,NPZ))/(DDB+SMALL)))**5,ZERO )
         WC(NPZ,NSL) = WC(NPZ,NSL)+(C(NB,NSL)*AG*FCGB-C(N,NSL)*AGP*FCGP)
       ENDDO
  300 CONTINUE
!
!---  End of SFXG group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXL( NSL )
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
!     Compute solute transport flux aqueous-phase, excluding boundaries,
!     using either a Patankar scheme or a TVD scheme with  third-order
!     Leonard limiting for the  advective transport component.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, October.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/SFXL'
      M = 1
!
!---  X-direction solute flux aqueous-phase, excluding boundaries
!
      DO 100 N = 1,NFLD
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 100
        DO NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
          TCOR = (T(2,N)+TABS)/TSPRF
          SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
          VMCP = SL(2,N)*PORD(2,N)
          FCLP = YL(N,NSL)/(VMCP+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFP = TORL(2,N)*VMCP*SDFLP
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFP = SDCL(1,IZ(N),NSL)*SDCL(2,IZ(N),NSL)*
     &        EXP(VMCP*SDCL(3,IZ(N),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFP = TORL(2,N)*VMCP*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFP = SDCL(1,IZ(N),NSL)*SDCL(2,IZ(N),NSL)*
     &        VMCP**SDCL(3,IZ(N),NSL)
          ENDIF
          TCOR = (T(2,NW)+TABS)/TSPRF
          SDFLW = SMDL(NSL)*TCOR*(VISRL/VISL(2,NW))
          VMCW = SL(2,NW)*PORD(2,NW)
          FCLW = YL(NW,NSL)/(VMCW+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFW = TORL(2,NW)*VMCW*SDFLW
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFW = SDCL(1,IZ(NW),NSL)*SDCL(2,IZ(NW),NSL)*
     &        EXP(VMCW*SDCL(3,IZ(NW),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFW = TORL(2,NW)*VMCW*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFW = SDCL(1,IZ(NW),NSL)*SDCL(2,IZ(NW),NSL)*
     &        VMCW**SDCL(3,IZ(NW),NSL)
          ENDIF
          IF( IDISP .EQ. 1 ) THEN
            CALL ADVW( PORD,SL,UL,VL,WL,ULX,VLX,WLX,N,M )
            ULWSQ = ULX*ULX
            VLWSQ = VLX*VLX
            WLWSQ = WLX*WLX
            ZVW = SQRT(ULWSQ+VLWSQ+WLWSQ)
            INDX = 17
            DISPLW = DISPL(IZ(NW))*SMDEF(IZ(NW),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTW = DISPT(IZ(NW))*SMDEF(IZ(NW),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLW = DIFMN(DISPLW,DISPLP,DXGF(NW),DXGF(N),ULX,INDX)
            DPTW = DIFMN(DISPTW,DISPTP,DXGF(NW),DXGF(N),ULX,INDX)
            DPW = (DPLW*ULWSQ + DPTW*(VLWSQ+WLWSQ))/(ZVW+SMALL)
          ELSE
            DPW = 0.D+0
          ENDIF
          INDX = 16
          DFW = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),UL(1,NPX),INDX)
          DDW = (DFW+DPW)/DXGP(NPX)
          AL = MAX( UL(1,NPX),ZERO ) +
     &      DDW*MAX( (ONE-(TENTH*ABS(UL(1,NPX))/(DDW+SMALL)))**5,ZERO )
          ALP = MAX( -UL(1,NPX),ZERO ) +
     &      DDW*MAX( (ONE-(TENTH*ABS(UL(1,NPX))/(DDW+SMALL)))**5,ZERO )
         UC(NPX,NSL) = UC(NPX,NSL)+(C(NW,NSL)*AL*FCLW-C(N,NSL)*ALP*FCLP)
        ENDDO
  100 CONTINUE
!
!---  Y-direction solute flux aqueous-phase, excluding boundaries
!
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
          TCOR = (T(2,N)+TABS)/TSPRF
          SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
          VMCP = SL(2,N)*PORD(2,N)
          FCLP = YL(N,NSL)/(VMCP+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFP = TORL(2,N)*VMCP*SDFLP
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFP = SDCL(1,IZ(N),NSL)*SDCL(2,IZ(N),NSL)*
     &        EXP(VMCP*SDCL(3,IZ(N),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFP = TORL(2,N)*VMCP*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFP = SDCL(1,IZ(N),NSL)*SDCL(2,IZ(N),NSL)*
     &        VMCP**SDCL(3,IZ(N),NSL)
          ENDIF
          TCOR = (T(2,NS)+TABS)/TSPRF
          SDFLS = SMDL(NSL)*TCOR*(VISRL/VISL(2,NS))
          VMCS = SL(2,NS)*PORD(2,NS)
          FCLS = YL(NS,NSL)/(VMCS+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFS = TORL(2,NS)*VMCS*SDFLS
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFS = SDCL(1,IZ(NS),NSL)*SDCL(2,IZ(NS),NSL)*
     &        EXP(VMCS*SDCL(3,IZ(NS),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFS = TORL(2,NS)*VMCS*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFS = SDCL(1,IZ(NS),NSL)*SDCL(2,IZ(NS),NSL)*
     &        VMCS**SDCL(3,IZ(NS),NSL)
          ENDIF
          IF( IDISP .EQ. 1 ) THEN
            CALL ADVS( PORD,SL,UL,VL,WL,ULSX,VLSX,WLSX,N,M )
            ULSSQ = ULSX*ULSX
            VLSSQ = VLSX*VLSX
            WLSSQ = WLSX*WLSX
            ZVS = SQRT(ULSSQ+VLSSQ+WLSSQ)
            INDX = 17
            DISPLS = DISPL(IZ(NS))*SMDEF(IZ(NS),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTS = DISPT(IZ(NS))*SMDEF(IZ(NS),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLS = DIFMN(DISPLS,DISPLP,DYGF(NS),DYGF(N),VLSX,INDX)
            DPTS = DIFMN(DISPTS,DISPTP,DYGF(NS),DYGF(N),VLSX,INDX)
            DPS = (DPLS*VLSSQ + DPTS*(ULSSQ+WLSSQ))/(ZVS+SMALL)
          ELSE
            DPS = 0.D+0
          ENDIF
          INDX = 16
          DFS = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),VL(1,NPY),INDX)
          DDS = (DFS+DPS)/(DYGP(NPY)*RP(I))
          AL = MAX( VL(1,NPY),ZERO ) +
     &      DDS*MAX( (ONE-(TENTH*ABS(VL(1,NPY))/(DDS+SMALL)))**5,ZERO )
          ALP = MAX( -VL(1,NPY),ZERO ) +
     &      DDS*MAX( (ONE-(TENTH*ABS(VL(1,NPY))/(DDS+SMALL)))**5,ZERO )
         VC(NPY,NSL) = VC(NPY,NSL)+(C(NS,NSL)*AL*FCLS-C(N,NSL)*ALP*FCLP)
        ENDDO
  200 CONTINUE
!
!---  Z-direction solute flux aqueous-phase, excluding boundaries
!
      DO 300 N = 1,NFLD
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 300
        DO NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
          TCOR = (T(2,N)+TABS)/TSPRF
          SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
          VMCP = SL(2,N)*PORD(2,N)
          FCLP = YL(N,NSL)/(VMCP+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFP = TORL(2,N)*VMCP*SDFLP
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFP = SDCL(1,IZ(N),NSL)*SDCL(2,IZ(N),NSL)*
     &        EXP(VMCP*SDCL(3,IZ(N),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFP = TORL(2,N)*VMCP*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFP = SDCL(1,IZ(N),NSL)*SDCL(2,IZ(N),NSL)*
     &        VMCP**SDCL(3,IZ(N),NSL)
          ENDIF
          TCOR = (T(2,NB)+TABS)/TSPRF
          SDFLB = SMDL(NSL)*TCOR*(VISRL/VISL(2,NB))
          VMCB = SL(2,NB)*PORD(2,NB)
          FCLB = YL(NB,NSL)/(VMCB+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFB = TORL(2,NB)*VMCB*SDFLB
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFB = SDCL(1,IZ(NB),NSL)*SDCL(2,IZ(NB),NSL)*
     &        EXP(VMCB*SDCL(3,IZ(NB),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFB = TORL(2,NB)*VMCB*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFB = SDCL(1,IZ(NB),NSL)*SDCL(2,IZ(NB),NSL)*
     &        VMCB**SDCL(3,IZ(NB),NSL)
          ENDIF
          IF( IDISP .EQ. 1 ) THEN
            CALL ADVB( PORD,SL,UL,VL,WL,ULBX,VLBX,WLBX,N,M )
            ULBSQ = ULBX*ULBX
            VLBSQ = VLBX*VLBX
            WLBSQ = WLBX*WLBX
            ZVB = SQRT(ULBSQ+VLBSQ+WLBSQ)
            INDX = 17
            DISPLB = DISPL(IZ(NB))*SMDEF(IZ(NB),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTB = DISPT(IZ(NB))*SMDEF(IZ(NB),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLB = DIFMN(DISPLB,DISPLP,DZGF(NB),DZGF(N),WLBX,INDX)
            DPTB = DIFMN(DISPTB,DISPTP,DZGF(NB),DZGF(N),WLBX,INDX)
            DPB = (DPLB*WLBSQ + DPTB*(VLBSQ+ULBSQ))/(ZVB+SMALL)
          ELSE
            DPB = 0.D+0
          ENDIF
          INDX = 16
          DFB = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
          DDB = (DFB+DPB)/DZGP(NPZ)
          AL = MAX( WL(1,NPZ),ZERO ) +
     &      DDB*MAX( (ONE-(TENTH*ABS(WL(1,NPZ))/(DDB+SMALL)))**5,ZERO )
          ALP = MAX( -WL(1,NPZ),ZERO ) +
     &      DDB*MAX( (ONE-(TENTH*ABS(WL(1,NPZ))/(DDB+SMALL)))**5,ZERO )
         WC(NPZ,NSL) = WC(NPZ,NSL)+(C(NB,NSL)*AL*FCLB-C(N,NSL)*ALP*FCLP)
        ENDDO 
  300 CONTINUE
!
!---  End of SFXL group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

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
      DO 100 N = 1,NFLD
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 100
        DO NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
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
     &        UNWX,INDX)
            DPTW = DIFMN(DISPT(IZ(NW)),DISPT(IZ(N)),DXGF(NW),DXGF(N),
     &        UNWX,INDX)
            DPW = (DPLW*UNWSQ + DPTW*(VNWSQ+WNWSQ))/(ZVW+SMALL)
          ELSE
            DPW = 0.D+0
          ENDIF
          INDX = 16
          DFW = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),UN(1,NPX),INDX)
          DDW = (DFW+DPW)/DXGP(NPX)
          AN = MAX( UN(1,NPX),ZERO ) +
     &      DDW*MAX( (ONE-(TENTH*ABS(UN(1,NPX))/(DDW+SMALL)))**5,ZERO )
          ANP = MAX( -UN(1,NPX),ZERO ) +
     &      DDW*MAX( (ONE-(TENTH*ABS(UN(1,NPX))/(DDW+SMALL)))**5,ZERO )
         UC(NPX,NSL) = UC(NPX,NSL)+(C(NW,NSL)*AN*FCNW-C(N,NSL)*ANP*FCNP)
        ENDDO
  100 CONTINUE
!
!---  Y-direction solute flux NAPL-phase, excluding boundaries
!
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
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
     &        VNSX,INDX)
            DPTS = DIFMN(DISPT(IZ(NS)),DISPT(IZ(N)),DYGF(NS),DYGF(N),
     &        VNSX,INDX)
            DPS = (DPLS*VNSSQ + DPTS*(UNSSQ+WNSSQ))/(ZVS+SMALL)
          ELSE
            DPS = 0.D+0
          ENDIF
          INDX = 16
          DFS = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),VN(1,NPY),INDX)
          DDS = (DFS+DPS)/DYGP(NPY)/RP(I)
          AN = MAX( VN(1,NPY),ZERO ) +
     &      DDS*MAX( (ONE-(TENTH*ABS(VN(1,NPY))/(DDS+SMALL)))**5,ZERO )
          ANP = MAX( -VN(1,NPY),ZERO ) +
     &      DDS*MAX( (ONE-(TENTH*ABS(VN(1,NPY))/(DDS+SMALL)))**5,ZERO )
         VC(NPY,NSL) = VC(NPY,NSL)+(C(NS,NSL)*AN*FCNS-C(N,NSL)*ANP*FCNP)
        ENDDO
  200 CONTINUE
!
!---  Z-direction solute flux NAPL-phase, excludinN boundaries
!
      DO 300 N = 1,NFLD
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 300
        DO NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
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
     &        WNBX,INDX)
            DPTB = DIFMN(DISPT(IZ(NB)),DISPT(IZ(N)),DZGF(NB),DZGF(N),
     &        WNBX,INDX)
            DPB = (DPLB*WNBSQ + DPTB*(VNBSQ+UNBSQ))/(ZVB+SMALL)
          ELSE
            DPB = 0.D+0
          ENDIF
          INDX = 16
          DFB = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX)
          DDB = (DFB+DPB)/DZGP(NPZ)
          AN = MAX( WN(1,NPZ),ZERO ) +
     &      DDB*MAX( (ONE-(TENTH*ABS(WN(1,NPZ))/(DDB+SMALL)))**5,ZERO )
          ANP = MAX( -WN(1,NPZ),ZERO ) +
     &      DDB*MAX( (ONE-(TENTH*ABS(WN(1,NPZ))/(DDB+SMALL)))**5,ZERO )
         WC(NPZ,NSL) = WC(NPZ,NSL)+(C(NB,NSL)*AN*FCNB-C(N,NSL)*ANP*FCNP)
        ENDDO
  300 CONTINUE
!
!---  End of SFXN group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXGB( NSL )
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
!     Compute solute gas-phase fluxes on boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1995.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PORMED
      USE HYST
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 BCX(LSPBC+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFXGB'
      SVGB = 0.D+0
      FCGB = 0.D+0
!
!---  Loop over number of specified boundary conditions  ---
!
      NBCT = MIN( NSL+LUK,NSOLU+LUK+1 )
      DO 200 NB = 1,NBC
        TMZ = TM
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 200
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) GOTO 200
        IF( IBCM(NB) .EQ. 1 ) THEN
!
!---      Solute transport  ---
!
          IF( NSL.LE.NSOLU ) THEN
            BCX(1) = BC(NSL+LBCU,1,MB)
!
!---      Reactive species transport  ---
!
          ELSE
            BCX(1) = 0.D+0
            DO 10 NSPX = 1,IBCSP(1,NB)
              MX = NSOLU+LBCU+NSPX
              BCX(NSPX+1) = BC(MX,1,MB)
   10       CONTINUE
          ENDIF
        ELSE
          DO 100 M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
              DTBC = MIN( BC(1,M,MB)-TM,DT )
              TFBC = (TM-5.D-1*DTBC-BC(1,M-1,MB))/
     &          (BC(1,M,MB)-BC(1,M-1,MB))
!
!---          Solute transport  ---
!
              IF( NSL.LE.NSOLU ) THEN
                BCX(1) = BC(NSL+LBCU,M-1,MB) +
     &            TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
                IF( IBCT(NBCT,NB).EQ.12 ) BCX(1) = CBO(NB,NSL)
!
!---          Reactive species transport  ---
!
              ELSE
                BCX(1) = 0.D+0
                DO 20 NSPX = 1,IBCSP(1,NB)
                  MX = NSOLU+LBCU+NSPX
                  BCX(NSPX+1) = BC(MX,M-1,MB) +
     &              TFBC*(BC(MX,M,MB)-BC(MX,M-1,MB))
                  IF( IBCT(NBCT,NB).EQ.12 ) BCX(NSPX) = CBO(NB,NSL)
   20           CONTINUE
              ENDIF
              GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 200
        ENDIF
  110   CONTINUE
        N = IBCN(NB)
        MF = 1
        IZN = IZ(N)
        MP = IXP(N)
        I = ID(N)
        J = JD(N)
        K = KD(N)
!
!---    Compute adjacent node phase fractions  ---
!
        SVGP = SG(2,N)*PORD(2,N)
        FCGP = 0.D+0
        IF( SVGP.GT.SMALL ) FCGP = YG(N,NSL)/SVGP
!
!---    Solute transport only, skip calculations for reactive
!       species transport  ---
!
        IF( NSL.LE.NSOLU ) THEN
!
!---      Compute boundary phase fractions  ---
!
          IF( IPCL(NSL).EQ.2 ) THEN
            SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))*
     &        SLB(2,NB)
          ELSE
            SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))
          ENDIF
          SVLB = SLB(2,NB)*PORDB(2,NB)
          SVGB = SGB(2,NB)*PORDB(2,NB)
          SVNB = SNB(2,NB)*PORDB(2,NB)
!
!---      Constant gas-aqueous partition coefficient  ---
!
          IF( IPCGL(NSL).EQ.0 ) THEN
            PCGLX = PCGL(1,NSL)
!
!---      Temperature dependent gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.1 ) THEN
            TK = TB(2,NB)+TABS
            PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &        + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &        + PCGL(5,NSL)*TK**2 )
!
!---      Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.2 ) THEN
            PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
          ENDIF
          PCGLX = MAX( PCGLX,1.D-20 )
          PCGLX = MIN( PCGLX,1.D+20 )
!
!---      Phase-volumetric concentration ratios  ---
!
          FCLB = 1.D+0/(SVSB + SVLB + SVNB/PCLN(1,NSL)
     &       + SVGB*PCGLX)
          FCGB = 1.D+0/((SVSB + SVLB + SVNB)/PCGLX + SVGB)
          FCNB = 1.D+0/((SVSB + SVLB + SVGB*PCGLX)*PCLN(1,NSL) + SVNB)
! 
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = SVLB*FCLB
          YGB(NB,NSL) = SVGB*FCGB
          YNB(NB,NSL) = SVNB*FCNB
!
!---      Convert boundary concentrations  ---
!
          IF( IBCT(NBCT,NB).EQ.8 .OR. IBCT(NBCT,NB).EQ.14
     &      .OR. IBCT(NBCT,NB).EQ.23 ) THEN
            BCX(1) = BCX(1)/(FCLB+SMALL)
          ELSEIF( IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.15 ) THEN
            BCX(1) = BCX(1)/(FCGB+SMALL)
          ELSEIF( IBCT(NBCT,NB).EQ.10 .OR. IBCT(NBCT,NB).EQ.16 ) THEN
            BCX(1) = BCX(1)/(FCNB+SMALL)
          ENDIF
        ELSE
!
!---      Skip for initial condition type boundary condition  ---
!
          IF( IBCT(NBCT,NB).NE.12 ) THEN
!
!---        Convert species concentrations to total-component
!           concentrations  ---
!
            IF( NSL.LE.NSOLU+NEQC ) THEN
              NEQ = NSL-NSOLU
              DO 130 NSP = 1,IEQ_C(1,NEQ)
                DO 120 NSPX = 1,IBCSP(1,NB)
                  IF( IBCSP(NSPX+1,NB).EQ.IEQ_C(NSP+1,NEQ) ) THEN
                    BCX(1) = BCX(1) + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                  ENDIF
  120           CONTINUE
  130         CONTINUE
!
!---        Convert species concentrations to total-kinetic
!           concentrations  ---
!
            ELSEIF( NSL.LE.NSOLU+NEQC+NEQK ) THEN
              NEQ = NSL-NSOLU-NEQC
              DO 150 NSP = 1,IEQ_K(1,NEQ)
                DO 140 NSPX = 1,IBCSP(1,NB)
                  IF( IBCSP(NSPX+1,NB).EQ.IEQ_K(NSP+1,NEQ) ) THEN
                    BCX(1) = BCX(1) + EQ_K(NSP,NEQ)*BCX(NSPX+1)
                  ENDIF
  140           CONTINUE
  150         CONTINUE
            ENDIF
          ENDIF
          SVLB = SLB(2,NB)*PORDB(2,NB)
          YLB(NB,NSL) = 1.D+0
          FCLB = 0.D+0
          IF( SVLB.GT.SMALL ) FCLB = YLB(NB,NSL)/SVLB
!
!---      Convert boundary phase concentrations to
!         volumetric concentrations  ---
!
          IF( IBCT(NBCT,NB).EQ.8 .OR. IBCT(NBCT,NB).EQ.14
     &      .OR. IBCT(NBCT,NB).EQ.23 ) THEN
            BCX(1) = BCX(1)*SVLB
          ENDIF
        ENDIF
!
!---    Diffusion coefficients at nodes adjacent to boundaries  ---
!
        IF( IEDL(NSL).EQ.1 ) THEN
          TCOR = (T(2,N)+TABS)/TSPRF
          PCOR = (PG(2,N)+PATM)/PATM
          SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          DGP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
        ENDIF
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB) .EQ. -3 ) THEN
          NPZ = NSZ(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGBX,VGBX,WGBX,N,MF )
            CALL SHDP( WGBX,UGBX,VGBX,DISPL(IZN),DISPT(IZN),DPGB )
          ELSE
            DPGB = 0.D+0
          ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
            INDX = 16
            DGZ = DIFMN(DGB,DGP,DZGF(N),DZGF(N),WG(1,NPZ),INDX)
            DGZ = (DGZ+DPGB)/(5.D-1*DZGF(N))
            AB = MAX( WG(1,NPZ),ZERO ) +
     &       DGZ*MAX((ONE-(TENTH*ABS(WG(1,NPZ))/(DGZ+SMALL)))**5,ZERO)
            AP = MAX( -WG(1,NPZ),ZERO ) +
     &       DGZ*MAX((ONE-(TENTH*ABS(WG(1,NPZ))/(DGZ+SMALL)))**5,ZERO)
            WC(NPZ,NSL) = WC(NPZ,NSL) + (BCX(1)*AB*FCGB
     &        - C(N,NSL)*AP*FCGP)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            AP = MAX( -WG(1,NPZ),ZERO )
            WC(NPZ,NSL) = WC(NPZ,NSL) - C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            AB = MAX( WL(1,NPZ),ZERO )
            WC(NPZ,NSL) = WC(NPZ,NSL) + BCX(1)*AB*FCGB
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. -2 ) THEN
          NPY = NSY(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGSX,VGSX,WGSX,N,MF )
            CALL SHDP( VGSX,WGSX,UGSX,DISPL(IZN),DISPT(IZN),DPGS )
          ELSE
            DPGS = 0.D+0
          ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
            INDX = 16
            DGY = DIFMN(DGB,DGP,DYGF(N),DYGF(N),VG(1,NPY),INDX)
            DGY = (DGY+DPGS)/RP(I)/(5.D-1*DYGF(N))
            AS = MAX( VG(1,NPY),ZERO ) +
     &       DGY*MAX((ONE-(TENTH*ABS(VG(1,NPY))/(DGY+SMALL)))**5,ZERO)
            AP = MAX( -VG(1,NPY),ZERO ) +
     &       DGY*MAX((ONE-(TENTH*ABS(VG(1,NPY))/(DGY+SMALL)))**5,ZERO)
            VC(NPY,NSL) = VC(NPY,NSL) + (BCX(1)*AS*FCGB
     &        - C(N,NSL)*AP*FCGP)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            AP = MAX( -VG(1,NPY),ZERO )
            VC(NPY,NSL) = VC(NPY,NSL) - C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            AS = MAX( VG(1,NPY),ZERO )
            VC(NPY,NSL) = VC(NPY,NSL) - BCX(1)*AS*FCGB
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. -1 ) THEN
          NPX = NSX(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGX,VGX,WGX,N,MF )
            CALL SHDP( UGX,VGX,WGX,DISPL(IZN),DISPT(IZN),DPGW )
            ELSE
              DPGW = 0.D+0
            ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
            INDX = 16
            DGX = DIFMN(DGB,DGP,DXGF(N),DXGF(N),UG(1,NPX),INDX)
            DGX = (DGX+DPGW)/(5.D-1*DXGF(N))
            AW = MAX( UG(1,NPX),ZERO ) +
     &       DGX*MAX((ONE-(TENTH*ABS(UG(1,NPX))/(DGX+SMALL)))**5,ZERO)
            AP = MAX( -UG(1,NPX),ZERO ) +
     &       DGX*MAX((ONE-(TENTH*ABS(UG(1,NPX))/(DGX+SMALL)))**5,ZERO)
            UC(NPX,NSL) = UC(NPX,NSL) + (BCX(1)*AW*FCGB
     &        - C(N,NSL)*AP*FCGP)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            AP = MAX( -UG(1,NPX),ZERO )
            UC(NPX,NSL) = UC(NPX,NSL) -  C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            AW = MAX( UG(1,NPX),ZERO )
            UC(NPX,NSL) = UC(NPX,NSL) + BCX(1)*AW*FCGB
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. 1 ) THEN
          NQX = NSX(N)+1
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGEX,VGEX,WGEX,N,MF )
            CALL SHDP( UGEX,VGEX,WGEX,DISPL(IZN),DISPT(IZN),DPGE )
          ELSE
            DPGE = 0.D+0
          ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
            INDX = 16
            DGX = DIFMN(DGP,DGB,DXGF(N),DXGF(N),UG(1,NQX),INDX)
            DGX = (DGX+DPGE)/(5.D-1*DXGF(N))
            AE = MAX( -UG(1,NQX),ZERO ) +
     &       DGX*MAX((ONE-(TENTH*ABS(UG(1,NQX))/(DGX+SMALL)))**5,ZERO)
            AP = MAX( UG(1,NQX),ZERO ) +
     &       DGX*MAX((ONE-(TENTH*ABS(UG(1,NQX))/(DGX+SMALL)))**5,ZERO)
            UC(NQX,NSL) = UC(NQX,NSL)+(C(N,NSL)*AP*FCGP-BCX(1)*AE*FCGB)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            AP = MAX( UG(1,NQX),ZERO )
            UC(NQX,NSL) = UC(NQX,NSL) + C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            AE = MAX( -UG(1,NQX),ZERO )
            UC(NQX,NSL) = UC(NQX,NSL) - BCX(1)*AE*FCGB
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. 2 ) THEN
          NQY = NSY(N)+IFLD
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGNX,VGNX,WGNX,N,MF )
            CALL SHDP( VGNX,WGNX,UGNX,DISPL(IZN),DISPT(IZN),DPGN )
          ELSE
            DPGN = 0.D+0
          ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
            INDX = 16
            DGY = DIFMN(DGP,DGB,DYGF(N),DYGF(N),VG(1,NQY),INDX)
            DGY = (DGY+DPGN)/RP(I)/(5.D-1*DYGF(N))
            AN = MAX( -VG(1,NQY),ZERO ) +
     &       DGY*MAX((ONE-(TENTH*ABS(VG(1,NQY))/(DGY+SMALL)))**5,ZERO)
            AP = MAX( VG(1,NQY),ZERO ) +
     &       DGY*MAX((ONE-(TENTH*ABS(VG(1,NQY))/(DGY+SMALL)))**5,ZERO)
            VC(NQY,NSL) = VC(NQY,NSL) + (C(N,NSL)*AP*FCGP
     &        - BCX(1)*AN*FCGB)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            AP = MAX( VG(1,NQY),ZERO )
            VC(NQY,NSL) = VC(NQY,NSL) + C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            AN = MAX( -VG(1,NQY),ZERO )
            VC(NQY,NSL) = VC(NQY,NSL) - BCX(1)*AN*FCGB
          ENDIF
!
!---    Top boundary
!
        ELSEIF( IBCD(NB) .EQ. 3 ) THEN
          NQZ = NSZ(N)+IJFLD
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGTX,VGTX,WGTX,N,MF )
            CALL SHDP( WGTX,UGTX,VGTX,DISPL(IZN),DISPT(IZN),DPGT )
          ELSE
            DPGT = 0.D+0
          ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
            INDX = 16
            DGZ = DIFMN(DGP,DGB,DZGF(N),DZGF(N),WG(1,NQZ),INDX)
            DGZ = (DGZ+DPGT)/(5.D-1*DZGF(N))
            AT = MAX( -WG(1,NQZ),ZERO ) +
     &        DGZ*MAX((ONE-(TENTH*ABS(WG(1,NQZ))/(DGZ+SMALL)))**5,ZERO)
            AP = MAX( WG(1,NQZ),ZERO ) +
     &       DGZ*MAX((ONE-(TENTH*ABS(WG(1,NQZ))/(DGZ+SMALL)))**5,ZERO)
            WC(NQZ,NSL) = WC(NQZ,NSL) + (C(N,NSL)*AP*FCGP
     &        - BCX(1)*AT*FCGB)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            AP = MAX( WG(1,NQZ),ZERO )
             WC(NQZ,NSL) = WC(NQZ,NSL) + C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            AT = MAX( -WG(1,NQZ),ZERO )
            WC(NQZ,NSL) = WC(NQZ,NSL) - BCX(1)*AT*FCGB
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  End of SFXGB group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXLB( NSL )
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
!     Compute solute aqueous-phase fluxes on boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1995.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PORMED
      USE GRID
      USE FLUXP
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
      REAL*8 BCX(LSPBC+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFXLB'
!
!---  Loop over number of specified boundary conditions  ---
!
      NBCT = MIN( NSL+LUK,NSOLU+LUK+1 )
      DO 200 NB = 1,NBC
        TMZ = TM
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 200
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) GOTO 200
        IF( IBCM(NB).EQ.1 ) THEN
!
!---      Solute transport  ---
!
          IF( NSL.LE.NSOLU ) THEN
            BCX(1) = BC(NSL+LBCU,1,MB)
!
!---      Reactive species transport  ---
!
          ELSE
            BCX(1) = 0.D+0
            DO 10 NSPX = 1,IBCSP(1,NB)
              MX = NSOLU+LBCU+NSPX
              BCX(NSPX+1) = BC(MX,1,MB)
   10       CONTINUE
          ENDIF
        ELSE
          DO 100 M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
              DTBC = MIN( BC(1,M,MB)-TM,DT )
              TFBC = (TM-5.D-1*DTBC-BC(1,M-1,MB))/
     &          (BC(1,M,MB)-BC(1,M-1,MB))
!
!---          Solute transport  ---
!
              IF( NSL.LE.NSOLU ) THEN
                BCX(1) = BC(NSL+LBCU,M-1,MB) +
     &            TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
                IF( IBCT(NBCT,NB).EQ.12 ) BCX(1) = CBO(NB,NSL)
!
!---          Reactive species transport  ---
!
              ELSE
                BCX(1) = 0.D+0
                DO 20 NSPX = 1,IBCSP(1,NB)
                  MX = NSOLU+LBCU+NSPX
                  BCX(NSPX+1) = BC(MX,M-1,MB) +
     &              TFBC*(BC(MX,M,MB)-BC(MX,M-1,MB))
                  IF( IBCT(NBCT,NB).EQ.12 ) BCX(NSPX) = CBO(NB,NSL)
   20           CONTINUE
              ENDIF
              GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 200
        ENDIF
  110   CONTINUE
        N = IBCN(NB)
        MF = 1
        IZN = IZ(N)
        MP = IXP(N)
        I = ID(N)
        J = JD(N)
        K = KD(N)
!
!---  Compute adjacent node phase fractions  ---
!
        SVLP = SL(2,N)*PORD(2,N)
        FCLP = 0.D+0
        IF( SVLP.GT.SMALL ) FCLP = YL(N,NSL)/SVLP
!
!---    Solute transport only, skip calculations for reactive
!       species transport  ---
!
        IF( NSL.LE.NSOLU ) THEN
!
!---      Compute boundary phase fractions  ---
!
          IF( IPCL(NSL).EQ.2 ) THEN
            SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))*
     &        SLB(2,NB)
          ELSE
            SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))
          ENDIF
          SVLB = SLB(2,NB)*PORDB(2,NB)
          SVGB = SGB(2,NB)*PORDB(2,NB)
          SVNB = SNB(2,NB)*PORDB(2,NB)
!
!---      Constant gas-aqueous partition coefficient  ---
!
          IF( IPCGL(NSL).EQ.0 ) THEN
            PCGLX = PCGL(1,NSL)
!
!---      Temperature dependent gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.1 ) THEN
            TK = TB(2,NB)+TABS
            PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &        + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &        + PCGL(5,NSL)*TK**2 )
!
!---      Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.2 ) THEN
            PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
          ENDIF
          PCGLX = MAX( PCGLX,1.D-20 )
          PCGLX = MIN( PCGLX,1.D+20 )
!
!---    Phase-volumetric concentration ratios  ---
!
          FCLB = 1.D+0/(SVSB + SVLB + SVNB/PCLN(1,NSL) 
     &      + SVGB*PCGLX)
          FCGB = 1.D+0/((SVSB + SVLB + SVNB)/PCGLX + SVGB)
          FCNB = 1.D+0/((SVSB + SVLB + SVGB*PCGLX)*PCLN(1,NSL) + SVNB)
!
!---    Phase mole fractions  ---
!
          YLB(NB,NSL) = SVLB*FCLB
          YGB(NB,NSL) = SVGB*FCGB
          YNB(NB,NSL) = SVNB*FCNB
!
!---    Convert boundary concentrations  ---
!
          IF( IBCT(NBCT,NB).EQ.8 .OR. IBCT(NBCT,NB).EQ.14
     &      .OR. IBCT(NBCT,NB).EQ.23 ) THEN
            BCX(1) = BCX(1)/(FCLB+SMALL)
          ELSEIF( IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.15 ) THEN
            BCX(1) = BCX(1)/(FCGB+SMALL)
          ELSEIF( IBCT(NBCT,NB).EQ.10 .OR. IBCT(NBCT,NB).EQ.16 ) THEN
            BCX(1) = BCX(1)/(FCNB+SMALL)
          ENDIF
        ELSE
!
!---      Skip for initial condition type boundary condition  ---
!
          IF( IBCT(NBCT,NB).NE.12 ) THEN
!
!---        Convert species concentrations to total-component
!           concentrations  ---
!
            IF( NSL.LE.NSOLU+NEQC ) THEN
              NEQ = NSL-NSOLU
              DO 130 NSP = 1,IEQ_C(1,NEQ)
                DO 120 NSPX = 1,IBCSP(1,NB)
                  IF( IBCSP(NSPX+1,NB).EQ.IEQ_C(NSP+1,NEQ) ) THEN
                    BCX(1) = BCX(1) + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                  ENDIF
  120           CONTINUE
  130         CONTINUE
!
!---        Convert species concentrations to total-kinetic
!           concentrations  ---
!
            ELSEIF( NSL.LE.NSOLU+NEQC+NEQK ) THEN
              NEQ = NSL-NSOLU-NEQC
              DO 150 NSP = 1,IEQ_K(1,NEQ)
                DO 140 NSPX = 1,IBCSP(1,NB)
                  IF( IBCSP(NSPX+1,NB).EQ.IEQ_K(NSP+1,NEQ) ) THEN
                    BCX(1) = BCX(1) + EQ_K(NSP,NEQ)*BCX(NSPX+1)
                  ENDIF
  140           CONTINUE
  150         CONTINUE
            ENDIF
          ENDIF
          SVLB = SLB(2,NB)*PORDB(2,NB)
          YLB(NB,NSL) = 1.D+0
          FCLB = 0.D+0
          IF( SVLB.GT.SMALL ) FCLB = YLB(NB,NSL)/SVLB
!
!---      Convert boundary phase concentrations to
!         volumetric concentrations  ---
!
          IF( IBCT(NBCT,NB).EQ.8 .OR. IBCT(NBCT,NB).EQ.14
     &      .OR. IBCT(NBCT,NB).EQ.23 ) THEN
            BCX(1) = BCX(1)*SVLB
          ENDIF
        ENDIF
!
!---  Diffusion coefficients at nodes adjacent to boundaries  ---
!
        IF( IEDL(NSL).EQ.1 ) THEN
          TCOR = (T(2,N)+TABS)/TSPRF
          SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
          DLP = TORL(2,N)*SL(2,N)*PORD(2,N)*SDFLP
        ELSEIF( IEDL(NSL).EQ.2 ) THEN
          DLP = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &      EXP(SL(2,N)*PORD(2,N)*SDCL(3,IZN,NSL))
        ELSEIF( IEDL(NSL).EQ.3 ) THEN
          DLP = TORL(2,N)*SL(2,N)*PORD(2,N)*SMDL(NSL)
        ELSEIF( IEDL(NSL).EQ.4 ) THEN
          DLP = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &      (SL(2,N)*PORD(2,N))**SDCL(3,IZN,NSL)
        ENDIF
!
!---  Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
          NPZ = NSZ(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULBX,VLBX,WLBX,N,MF )
            CALL SHDP( WLBX,ULBX,VLBX,DISPL(IZN),DISPT(IZN),DPLB )
!            ULBX = (0.5D+0*(UL(1,NSX(N))+UL(1,NSX(N)+1)))**2
!            VLBX = (0.5D+0*(VL(1,NSY(N))+VL(1,NSY(N)+IFLD)))**2
!            WLBX = (WL(1,NPZ))**2
!            ZLB = SQRT(ULBX + VLBX + WLBX)
!            DPLB = (DISPL(IZN)*WLBX + DISPT(IZN)*(ULBX+VLBX))/
!     &        (ZLB+SMALL)
            DPLB = DPLB*SMDEF(IZN,NSL)
          ELSE
            DPLB = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            IF( IEDL(NSL).EQ.1 ) THEN
              TCOR = (TB(2,NB)+TABS)/TSPRF
              SDFLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
              DLB = TORLB(2,NB)*SVLB*SDFLB
            ELSEIF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SVLB*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SVLB*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          SVLB**SDCL(3,IZN,NSL)
            ENDIF
            INDX = 16
            DLZ = DIFMN(DLB,DLP,DZGF(N),DZGF(N),WL(1,NPZ),INDX)
            DLZ = (DLZ+DPLB)/(5.D-1*DZGF(N))
            AB = MAX( WL(1,NPZ),ZERO ) +
     &       DLZ*MAX((ONE-(TENTH*ABS(WL(1,NPZ))/(DLZ+SMALL)))**5,ZERO)
            AP = MAX( -WL(1,NPZ),ZERO ) +
     &       DLZ*MAX((ONE-(TENTH*ABS(WL(1,NPZ))/(DLZ+SMALL)))**5,ZERO)
            WC(NPZ,NSL) = WC(NPZ,NSL)+(BCX(1)*AB*FCLB-C(N,NSL)*AP*FCLP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(NBCT,NB).EQ.7 .OR.
     &      ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (WL(1,NPZ)/EPSL.LT.-EPSL)) ) THEN
            AP = MAX( -WL(1,NPZ),ZERO )
            WC(NPZ,NSL) = WC(NPZ,NSL) -  C(N,NSL)*AP*FCLP
!
!---  Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16
     &       .OR. ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (WL(1,NPZ)/EPSL.GT.EPSL)) ) THEN
            AB = MAX( WL(1,NPZ),ZERO )
            WC(NPZ,NSL) = WC(NPZ,NSL) + BCX(1)*AB*FCLB
          ENDIF
!
!---  South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          NPY = NSY(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULSX,VLSX,WLSX,N,MF )
            CALL SHDP( VLSX,WLSX,ULSX,DISPL(IZN),DISPT(IZN),DPLS )
!            ULSX = (0.5D+0*(UL(1,NSX(N))+UL(1,NSX(N)+1)))**2
!            VLSX = VL(1,NPY)**2
!            WLSX = (0.5D+0*(WL(1,NSZ(N))+WL(1,NSZ(N)+IJFLD)))**2
!            ZLS = SQRT(ULSX + VLSX + WLSX)
!            DPLS = (DISPL(IZN)*VLSX + DISPT(IZN)*(ULSX+WLSX))/
!     &        (ZLS+SMALL)
            DPLS = DPLS*SMDEF(IZN,NSL)
          ELSE
            DPLS = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR. 
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            IF( IEDL(NSL).EQ.1 ) THEN
              TCOR = (TB(2,NB)+TABS)/TSPRF
              SDFLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
              DLB = TORLB(2,NB)*SVLB*SDFLB
            ELSEIF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SVLB*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SVLB*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          SVLB**SDCL(3,IZN,NSL)
            ENDIF
            INDX = 16
            DLY = DIFMN(DLB,DLP,DYGF(N),DYGF(N),VL(1,NPY),INDX)
            DLY = (DLY+DPLS)/RP(I)/(5.D-1*DYGF(N))
            AS = MAX( VL(1,NPY),ZERO ) +
     &       DLY*MAX((ONE-(TENTH*ABS(VL(1,NPY))/(DLY+SMALL)))**5,ZERO)
            AP = MAX( -VL(1,NPY),ZERO ) +
     &       DLY*MAX((ONE-(TENTH*ABS(VL(1,NPY))/(DLY+SMALL)))**5,ZERO)
            VC(NPY,NSL) = VC(NPY,NSL)+(BCX(1)*AS*FCLB-C(N,NSL)*AP*FCLP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(NBCT,NB).EQ.7 .OR.
     &      ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (VL(1,NPY)/EPSL.LT.-EPSL)) ) THEN
            AP = MAX( -VL(1,NPY),ZERO )
            VC(NPY,NSL) = VC(NPY,NSL) - C(N,NSL)*AP*FCLP
!
!---  Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16
     &       .OR. ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (VL(1,NPY)/EPSL.GT.EPSL)) ) THEN
            AS = MAX( VL(1,NPY),ZERO )
            VC(NPY,NSL) = VC(NPY,NSL) - BCX(1)*AS*FCLB
          ENDIF
!
!---  West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          NPX = NSX(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULX,VLX,WLX,N,MF )
            CALL SHDP( ULX,VLX,WLX,DISPL(IZN),DISPT(IZN),DPLW )
!            ULX = UL(1,NPX)**2
!            VLX = (0.5D+0*(VL(1,NSY(N))+VL(1,NSY(N)+IFLD)))**2
!            WLX = (0.5D+0*(WL(1,NSZ(N))+WL(1,NSZ(N)+IJFLD)))**2
!            ZLW = SQRT(ULX + VLX + WLX)
!            DPLW = (DISPL(IZN)*ULX + DISPT(IZN)*(WLX+VLX))/
!     &          (ZLW+SMALL)
            DPLW = DPLW*SMDEF(IZN,NSL)
            ELSE
              DPLW = 0.D+0
            ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            IF( IEDL(NSL).EQ.1 ) THEN
              TCOR = (TB(2,NB)+TABS)/TSPRF
              SDFLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
              DLB = TORLB(2,NB)*SVLB*SDFLB
            ELSEIF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SVLB*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SVLB*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          SVLB**SDCL(3,IZN,NSL)
            ENDIF
            INDX = 16
            DLX = DIFMN(DLB,DLP,DXGF(N),DXGF(N),UL(1,NPX),INDX)
            DLX = (DLX+DPLW)/(5.D-1*DXGF(N))
            AW = MAX( UL(1,NPX),ZERO ) +
     &       DLX*MAX((ONE-(TENTH*ABS(UL(1,NPX))/(DLX+SMALL)))**5,ZERO)
            AP = MAX( -UL(1,NPX),ZERO ) +
     &       DLX*MAX((ONE-(TENTH*ABS(UL(1,NPX))/(DLX+SMALL)))**5,ZERO)
            UC(NPX,NSL) = UC(NPX,NSL)+(BCX(1)*AW*FCLB-C(N,NSL)*AP*FCLP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(NBCT,NB).EQ.7 .OR.
     &      ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (UL(1,NPX)/EPSL.LT.-EPSL)) ) THEN
            AP = MAX( -UL(1,NPX),ZERO )
            UC(NPX,NSL) = UC(NPX,NSL) -  C(N,NSL)*AP*FCLP
!
!---  Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16
     &       .OR. ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (UL(1,NPX)/EPSL.GT.EPSL)) ) THEN
            AW = MAX( UL(1,NPX),ZERO )
            UC(NPX,NSL) = UC(NPX,NSL) + BCX(1)*AW*FCLB
          ENDIF
!
!---  East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          NQX = NSX(N)+1
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULEX,VLEX,WLEX,N,MF )
            CALL SHDP( ULEX,VLEX,WLEX,DISPL(IZN),DISPT(IZN),DPLE )
!            ULEX = UL(1,NQX)**2
!            VLEX = (0.5D+0*(VL(1,NSY(N))+VL(1,NSY(N)+IFLD)))**2
!            WLEX = (0.5D+0*(WL(1,NSZ(N))+WL(1,NSZ(N)+IJFLD)))**2
!            ZLE = SQRT(ULEX + VLEX + WLEX)
!            DPLE = (DISPL(IZN)*ULEX + DISPT(IZN)*(WLEX+VLEX))/
!     &        (ZLE+SMALL)
            DPLE = DPLE*SMDEF(IZN,NSL)
          ELSE
            DPLE = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            IF( IEDL(NSL).EQ.1 ) THEN
              TCOR = (TB(2,NB)+TABS)/TSPRF
              SDFLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
              DLB = TORLB(2,NB)*SVLB*SDFLB
            ELSEIF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SVLB*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SVLB*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          SVLB**SDCL(3,IZN,NSL)
            ENDIF
            INDX = 16
            DLX = DIFMN(DLP,DLB,DXGF(N),DXGF(N),UL(1,NQX),INDX)
            DLX = (DLX+DPLE)/(5.D-1*DXGF(N))
            AE = MAX( -UL(1,NQX),ZERO ) +
     &       DLX*MAX((ONE-(TENTH*ABS(UL(1,NQX))/(DLX+SMALL)))**5,ZERO)
            AP = MAX( UL(1,NQX),ZERO ) +
     &       DLX*MAX((ONE-(TENTH*ABS(UL(1,NQX))/(DLX+SMALL)))**5,ZERO)
            UC(NQX,NSL) = UC(NQX,NSL)+(C(N,NSL)*AP*FCLP-BCX(1)*AE*FCLB)
!
!---  Outflow ---
!
          ELSEIF( IBCT(NBCT,NB).EQ.7 .OR.
     &      ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (UL(1,NQX)/EPSL.GT.EPSL)) ) THEN
            AP = MAX( UL(1,NQX),ZERO )
            UC(NQX,NSL) = UC(NQX,NSL) + C(N,NSL)*AP*FCLP
!
!---  Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16
     &       .OR. ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (UL(1,NQX)/EPSL.LT.-EPSL)) ) THEN
            AE = MAX( -UL(1,NQX),ZERO )
            UC(NQX,NSL) = UC(NQX,NSL) - BCX(1)*AE*FCLB
          ENDIF
!
!---  North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          NQY = NSY(N)+IFLD
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULNX,VLNX,WLNX,N,MF )
            CALL SHDP( VLNX,WLNX,ULNX,DISPL(IZN),DISPT(IZN),DPLN )
!            ULNX = (0.5D+0*(UL(1,NSX(N))+UL(1,NSX(N)+1)))**2
!            VLNX = VL(1,NQY)**2
!            WLNX = (0.5D+0*(WL(1,NSZ(N))+WL(1,NSZ(N)+IJFLD)))**2
!            ZLN = SQRT(ULNX + VLNX + WLNX)
!            DPLN = (DISPL(IZN)*VLNX + DISPT(IZN)*(ULNX+WLNX))/
!     &        (ZLN+SMALL)
            DPLN = DPLN*SMDEF(IZN,NSL)
          ELSE
            DPLN = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            IF( IEDL(NSL).EQ.1 ) THEN
              TCOR = (TB(2,NB)+TABS)/TSPRF
              SDFLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
              DLB = TORLB(2,NB)*SVLB*SDFLB
            ELSEIF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SVLB*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SVLB*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          SVLB**SDCL(3,IZN,NSL)
            ENDIF
            INDX = 16
            DLY = DIFMN(DLP,DLB,DYGF(N),DYGF(N),VL(1,NQY),INDX)
            DLY = (DLY+DPLN)/RP(I)/(5.D-1*DYGF(N))
            AN = MAX( -VL(1,NQY),ZERO ) +
     &       DLY*MAX((ONE-(TENTH*ABS(VL(1,NQY))/(DLY+SMALL)))**5,ZERO)
            AP = MAX( VL(1,NQY),ZERO ) +
     &       DLY*MAX((ONE-(TENTH*ABS(VL(1,NQY))/(DLY+SMALL)))**5,ZERO)
            VC(NQY,NSL) = VC(NQY,NSL)+(C(N,NSL)*AP*FCLP-BCX(1)*AN*FCLB)
!
!---  Outflow ---
!
          ELSEIF( IBCT(NBCT,NB).EQ.7 .OR.
     &      ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (VL(1,NQY)/EPSL.GT.EPSL)) ) THEN
            AP = MAX( VL(1,NQY),ZERO )
            VC(NQY,NSL) = VC(NQY,NSL) + C(N,NSL)*AP*FCLP
!
!---  Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16
     &       .OR. ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (VL(1,NQY)/EPSL.LT.-EPSL)) ) THEN
            AN = MAX( -VL(1,NQY),ZERO )
            VC(NQY,NSL) = VC(NQY,NSL) - BCX(1)*AN*FCLB
          ENDIF
!
!---  Top boundary
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          NQZ = NSZ(N)+IJFLD
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULTX,VLTX,WLTX,N,MF )
            CALL SHDP( WLTX,ULTX,VLTX,DISPL(IZN),DISPT(IZN),DPLT )
!            ULTX = (0.5D+0*(UL(1,NSX(N))+UL(1,NSX(N)+1)))**2
!            VLTX = (0.5D+0*(VL(1,NSY(N))+VL(1,NSY(N)+IFLD)))**2
!            WLTX = (WL(1,NQZ))**2
!            ZLT = SQRT(ULTX + VLTX + WLTX)
!            DPLT = (DISPL(IZN)*WLTX + DISPT(IZN)*(ULTX+VLTX))/
!     &        (ZLT+SMALL)
            DPLT = DPLT*SMDEF(IZN,NSL)
          ELSE
            DPLT = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            IF( IEDL(NSL).EQ.1 ) THEN
              TCOR = (TB(2,NB)+TABS)/TSPRF
              SDFLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
              DLB = TORLB(2,NB)*SVLB*SDFLB
            ELSEIF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SVLB*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SVLB*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          SVLB**SDCL(3,IZN,NSL)
            ENDIF
            INDX = 16
            DLZ = DIFMN(DLP,DLB,DZGF(N),DZGF(N),WL(1,NQZ),INDX)
            DLZ = (DLZ+DPLT)/(5.D-1*DZGF(N))
            AT = MAX( -WL(1,NQZ),ZERO ) +
     &       DLZ*MAX((ONE-(TENTH*ABS(WL(1,NQZ))/(DLZ+SMALL)))**5,ZERO)
            AP = MAX( WL(1,NQZ),ZERO ) +
     &       DLZ*MAX((ONE-(TENTH*ABS(WL(1,NQZ))/(DLZ+SMALL)))**5,ZERO)
            WC(NQZ,NSL) = WC(NQZ,NSL)+(C(N,NSL)*AP*FCLP-BCX(1)*AT*FCLB)
!
!---  Outflow ---
!
          ELSEIF( IBCT(NBCT,NB).EQ.7 .OR.
     &      ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (WL(1,NQZ)/EPSL.GT.EPSL)) ) THEN
            AP = MAX( WL(1,NQZ),ZERO )
            WC(NQZ,NSL) = WC(NQZ,NSL) + C(N,NSL)*AP*FCLP
!
!---  Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16
     &       .OR. ((IBCT(NBCT,NB).EQ.19 .OR. IBCT(NBCT,NB).EQ.23)
     &      .AND. (WL(1,NQZ)/EPSL.LT.-EPSL)) ) THEN
            AT = MAX( -WL(1,NQZ),ZERO )
            WC(NQZ,NSL) = WC(NQZ,NSL) - BCX(1)*AT*FCLB
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  End of SFXLB group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXNB( NSL )
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
!     Compute solute NAPL fluxes on boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1995.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE NAPL
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FDVP
      USE FDVN
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/SFXNB'
!
!---  Loop over number of specified boundary conditions  ---
!
      DO 200 NB = 1,NBC
        TMZ = TM
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 200
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) GOTO 200
        IF( IBCM(NB) .EQ. 1 ) THEN
          BCX = BC(NSL+LBCU,1,MB)
        ELSE
          DO 100 M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
             DTBC = MIN( BC(1,M,MB)-TM,DT )
             TFBC = (TM-5.D-1*DTBC-BC(1,M-1,MB))/
     &         (BC(1,M,MB)-BC(1,M-1,MB))
             BCX = BC(NSL+LBCU,M-1,MB) +
     &         TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 200
        ENDIF
  110   CONTINUE
        N = IBCN(NB)
        MF = 1
        IZN = IZ(N)
        MP = IXP(N)
        I = ID(N)
        J = JD(N)
        K = KD(N)
!
!---  Compute adjacent node phase fractions  ---
!
        SVNP = SN(2,N)*PORD(2,N)
        FCNP = 0.D+0
        IF( SVNP.GT.SMALL ) FCNP = YN(N,NSL)/SVNP
!
!---  Compute boundary phase fractions  ---
!
        IF( IPCL(NSL).EQ.2 ) THEN
          SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))*SLB(2,NB)
        ELSE
          SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))
        ENDIF
        SVLB = SLB(2,NB)*PORDB(2,NB)
        SVGB = SGB(2,NB)*PORDB(2,NB)
        SVNB = SNB(2,NB)*PORDB(2,NB)
!
!---    Constant gas-aqueous partition coefficient  ---
!
        IF( IPCGL(NSL).EQ.0 ) THEN
          PCGLX = PCGL(1,NSL)
!
!---    Temperature dependent gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.1 ) THEN
          TK = TB(2,NB)+TABS
          PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &      + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &      + PCGL(5,NSL)*TK**2 )
!
!---    Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.2 ) THEN
          PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
        ENDIF
        PCGLX = MAX( PCGLX,1.D-20 )
        PCGLX = MIN( PCGLX,1.D+20 )
!
!---  Phase-volumetric concentration ratios  ---
!
        FCLB = 1.D+0/(SVSB + SVLB + SVNB/PCLN(1,NSL) + SVGB*PCGLX)
        FCGB = 1.D+0/((SVSB + SVLB + SVNB)/PCGLX + SVGB)
        FCNB = 1.D+0/((SVSB + SVLB + SVGB*PCGLX)*PCLN(1,NSL) + SVNB)
!
!---  Phase mole fractions  ---
!
        YLB(NB,NSL) = SVLB*FCLB
        YGB(NB,NSL) = SVGB*FCGB
        YNB(NB,NSL) = SVNB*FCNB
!
!---  Convert boundary concentrations  ---
!
        IF( IBCT(IEQC,NB).EQ.8 .OR. IBCT(IEQC,NB).EQ.14 ) THEN
          BCX = BCX/(FCLB+SMALL)
        ELSEIF( IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.15 ) THEN
          BCX = BCX/(FCGB+SMALL)
        ELSEIF( IBCT(IEQC,NB).EQ.10 .OR. IBCT(IEQC,NB).EQ.16 ) THEN
          BCX = BCX/(FCNB+SMALL)
        ENDIF
!
!---  Diffusion coefficients at nodes adjacent to boundaries  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        SDFNP = SMDN(NSL)*TCOR*(VISRN/VISN(2,N))
        DNP = TORN(2,N)*SN(2,N)*PORD(2,N)*SDFNP
!
!---  Bottom boundary  ---
!
        IF( IBCD(NB) .EQ. -3 ) THEN
          NPZ = NSZ(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNBX,VNBX,WNBX,N,MF )
            CALL SHDP( WNBX,UNBX,VNBX,DISPL(IZN),DISPT(IZN),DPNB )
!            UNBX = (0.5D+0*(UN(1,NSX(N))+UN(1,NSX(N)+1)))**2
!            VNBX = (0.5D+0*(VN(1,NSY(N))+VN(1,NSY(N)+IFLD)))**2
!            WNBX = (WN(1,NPZ))**2
!            ZNB = SQRT(UNBX + VNBX + WNBX)
!            DPNB = (DISPL(IZN)*WNBX + DISPT(IZN)*(UNBX+VNBX))/
!     &        (ZNB+SMALL)
          ELSE
            DPNB = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNZ = DIFMN(DNB,DNP,DZGF(N),DZGF(N),WN(1,NPZ),INDX)
            DNZ = (DNZ+DPNB)/(5.D-1*DZGF(N))
            AB = MAX( WN(1,NPZ),ZERO ) +
     &       DNZ*MAX((ONE-(TENTH*ABS(WN(1,NPZ))/(DNZ+SMALL)))**5,ZERO)
            AP = MAX( -WN(1,NPZ),ZERO ) +
     &       DNZ*MAX((ONE-(TENTH*ABS(WN(1,NPZ))/(DNZ+SMALL)))**5,ZERO)
            WC(NPZ,NSL) = WC(NPZ,NSL) + (BCX*AB*FCNB - C(N,NSL)*AP*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            AP = MAX( -WN(1,NPZ),ZERO )
            WC(NPZ,NSL) = WC(NPZ,NSL) - C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            AB = MAX( WL(1,NPZ),ZERO )
             WC(NPZ,NSL) = WC(NPZ,NSL) + BCX*AB*FCNP
          ENDIF
!
!---  South boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. -2 ) THEN
          NPY = NSY(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNSX,VNSX,WNSX,N,MF )
            CALL SHDP( VNSX,WNSX,UNSX,DISPL(IZN),DISPT(IZN),DPNS )
!            UNSX = (0.5D+0*(UN(1,NSX(N))+UN(1,NSX(N)+1)))**2
!            VNSX = VN(1,NPY)**2
!            WNSX = (0.5D+0*(WN(1,NSZ(N))+WN(1,NSZ(N)+IJFLD)))**2
!            ZNS = SQRT(UNSX + VNSX + WNSX)
!            DPNS = (DISPL(IZN)*VNSX + DISPT(IZN)*(UNSX+WNSX))/
!     &        (ZNS+SMALL)
          ELSE
            DPNS = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNY = DIFMN(DNB,DNP,DYGF(N),DYGF(N),VN(1,NPY),INDX)
            DNY = (DNY+DPNS)/RP(I)/(5.D-1*DYGF(N))
            AS = MAX( VN(1,NPY),ZERO ) +
     &       DNY*MAX((ONE-(TENTH*ABS(VN(1,NPY))/(DNY+SMALL)))**5,ZERO)
            AP = MAX( -VN(1,NPY),ZERO ) +
     &       DNY*MAX((ONE-(TENTH*ABS(VN(1,NPY))/(DNY+SMALL)))**5,ZERO)
            VC(NPY,NSL) = VC(NPY,NSL) + (BCX*AS*FCNP - C(N,NSL)*AP*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            AP = MAX( -VN(1,NPY),ZERO )
            VC(NPY,NSL) = VC(NPY,NSL) - C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            AS = MAX( VN(1,NPY),ZERO )
            VC(NPY,NSL) = VC(NPY,NSL) - BCX*AS*FCNP
          ENDIF
!
!---  West boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. -1 ) THEN
          NPX = NSX(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNWX,VNWX,WNWX,N,MF )
            CALL SHDP( UNWX,VNWX,WNWX,DISPL(IZN),DISPT(IZN),DPNW )
!            UNWX = UN(1,NPX)**2
!            VNWX = (0.5D+0*(VN(1,NSY(N))+VN(1,NSY(N)+IFLD)))**2
!            WNWX = (0.5D+0*(WN(1,NSZ(N))+WN(1,NSZ(N)+IJFLD)))**2
!            ZNW = SQRT(UNWX + VNWX + WNWX)
!            DPNW = (DISPL(IZN)*UNWX + DISPT(IZN)*(WNWX+VNWX))/
!     &          (ZNW+SMALL)
            ELSE
              DPNW = 0.D+0
            ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNX = DIFMN(DNB,DNP,DXGF(N),DXGF(N),UN(1,NPX),INDX)
            DNX = (DNX+DPNW)/(5.D-1*DXGF(N))
            AW = MAX( UN(1,NPX),ZERO ) +
     &       DNX*MAX((ONE-(TENTH*ABS(UN(1,NPX))/(DNX+SMALL)))**5,ZERO)
            AP = MAX( -UN(1,NPX),ZERO ) +
     &       DNX*MAX((ONE-(TENTH*ABS(UN(1,NPX))/(DNX+SMALL)))**5,ZERO)
            UC(NPX,NSL) = UC(NPX,NSL) + (BCX*AW*FCNP - C(N,NSL)*AP*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            AP = MAX( -UN(1,NPX),ZERO )
            UC(NPX,NSL) = UC(NPX,NSL) -  C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            AW = MAX( UN(1,NPX),ZERO )
            UC(NPX,NSL) = UC(NPX,NSL) + BCX*AW*FCNP
          ENDIF
!
!---  East boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. 1 ) THEN
          NQX = NSX(N)+1
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNEX,VNEX,WNEX,N,MF )
            CALL SHDP( UNEX,VNEX,WNEX,DISPL(IZN),DISPT(IZN),DPNE )
!            UNEX = UN(1,NQX)**2
!            VNEX = (0.5D+0*(VN(1,NSY(N))+VN(1,NSY(N)+IFLD)))**2
!            WNEX = (0.5D+0*(WN(1,NSZ(N))+WN(1,NSZ(N)+IJFLD)))**2
!            ZNE = SQRT(UNEX + VNEX + WNEX)
!            DPNE = (DISPL(IZN)*UNEX + DISPT(IZN)*(WNEX+VNEX))/
!     &        (ZNE+SMALL)
          ELSE
            DPNE = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNX = DIFMN(DNP,DNB,DXGF(N),DXGF(N),UN(1,NQX),INDX)
            DNX = (DNX+DPNE)/(5.D-1*DXGF(N))
            AE = MAX( -UN(1,NQX),ZERO ) +
     &       DNX*MAX((ONE-(TENTH*ABS(UN(1,NQX))/(DNX+SMALL)))**5,ZERO)
            AP = MAX( UN(1,NQX),ZERO ) +
     &       DNX*MAX((ONE-(TENTH*ABS(UN(1,NQX))/(DNX+SMALL)))**5,ZERO)
            UC(NQX,NSL) = UC(NQX,NSL)+(C(N,NSL)*AP*FCNP-BCX*AE*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            AP = MAX( UN(1,NQX),ZERO )
            UC(NQX,NSL) = UC(NQX,NSL) + C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            AE = MAX( -UN(1,NQX),ZERO )
            UC(NQX,NSL) = UC(NQX,NSL) - BCX*AE*FCNP
          ENDIF
!
!---  North boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. 2 ) THEN
          NQY = NSY(N)+IFLD
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNNX,VNNX,WNNX,N,MF )
            CALL SHDP( VNNX,WNNX,UNNX,DISPL(IZN),DISPT(IZN),DPNN )
!            UNNX = (0.5D+0*(UN(1,NSX(N))+UN(1,NSX(N)+1)))**2
!            VNNX = VN(1,NQY)**2
!            WNNX = (0.5D+0*(WN(1,NSZ(N))+WN(1,NSZ(N)+IJFLD)))**2
!            ZNN = SQRT(UNNX + VNNX + WNNX)
!            DPNN = (DISPL(IZN)*VNNX + DISPT(IZN)*(UNNX+WNNX))/
!     &        (ZNN+SMALL)
          ELSE
            DPNN = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNY = DIFMN(DNP,DNB,DYGF(N),DYGF(N),VN(1,NQY),INDX)
            DNY = (DNY+DPNN)/RP(I)/(5.D-1*DYGF(N))
            AN = MAX( -VN(1,NQY),ZERO ) +
     &       DNY*MAX((ONE-(TENTH*ABS(VN(1,NQY))/(DNY+SMALL)))**5,ZERO)
            AP = MAX( VN(1,NQY),ZERO ) +
     &       DNY*MAX((ONE-(TENTH*ABS(VN(1,NQY))/(DNY+SMALL)))**5,ZERO)
            VC(NQY,NSL) = VC(NQY,NSL) + (C(N,NSL)*AP*FCNP - BCX*AN*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            AP = MAX( VN(1,NQY),ZERO )
            VC(NQY,NSL) = VC(NQY,NSL) + C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            AN = MAX( -VN(1,NQY),ZERO )
            VC(NQY,NSL) = VC(NQY,NSL) - BCX*AN*FCNP
          ENDIF
!
!---  Top boundary
!
        ELSEIF( IBCD(NB) .EQ. 3 ) THEN
          NQZ = NSZ(N)+IJFLD
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNTX,VNTX,WNTX,N,MF )
            CALL SHDP( WNTX,UNTX,VNTX,DISPL(IZN),DISPT(IZN),DPNT )
!            UNTX = (0.5D+0*(UN(1,NSX(N))+UN(1,NSX(N)+1)))**2
!            VNTX = (0.5D+0*(VN(1,NSY(N))+VN(1,NSY(N)+IFLD)))**2
!            WNTX = (WN(1,NQZ))**2
!            ZNT = SQRT(UNTX + VNTX + WNTX)
!            DPNT = (DISPL(IZN)*WNTX + DISPT(IZN)*(UNTX+VNTX))/
!     &        (ZNT+SMALL)
          ELSE
            DPNT = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNZ = DIFMN(DNP,DNB,DZGF(N),DZGF(N),WN(1,NQZ),INDX)
            DNZ = (DNZ+DPNT)/(5.D-1*DZGF(N))
            AT = MAX( -WN(1,NQZ),ZERO ) +
     &       DNZ*MAX((ONE-(TENTH*ABS(WN(1,NQZ))/(DNZ+SMALL)))**5,ZERO)
            AP = MAX( WN(1,NQZ),ZERO ) +
     &       DNZ*MAX((ONE-(TENTH*ABS(WN(1,NQZ))/(DNZ+SMALL)))**5,ZERO)
            WC(NQZ,NSL) = WC(NQZ,NSL) + (C(N,NSL)*AP*FCNP - BCX*AT*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            AP = MAX( WN(1,NQZ),ZERO )
            WC(NQZ,NSL) = WC(NQZ,NSL) + C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            AT = MAX( -WN(1,NQZ),ZERO )
            WC(NQZ,NSL) = WC(NQZ,NSL) - BCX*AT*FCNP
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  End of SFXNB group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END


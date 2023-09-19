!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VFLGB
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
!     Water-Air-Energy-VOC Transport Mode
!     Compute VOC transport fluxes on boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, April 8, 1994.
!     vflgb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VFLGB'
      M = 1
!
!---  Loop over boundary conditions  ---
!
      DO 100 NB = 1,NBC
        MB = IBCIN(NB)
        IF( TM.LT.BC(1,1,MB) ) GOTO 100
        IF( IBCM(NB).GT.1 .AND. TM.GT.BC(1,IBCM(NB),MB) ) GOTO 100
        N = IBCN(NB)
        IZN = IZ(N)
        MP = IXP(N)
!
!---  Compute diffusion coefficients  ---
!
        DGP = TORG(2,N)*PORD(2,N)*SG(2,N)*DFGO(2,N)
        DLP = TORL(2,N)*PORD(2,N)*SL(2,N)*DFLO(2,N)
!
!---  Compute phase fractions  ---
!
        FCL = 0.D+0
        FCG = 0.D+0
        IF( SL(2,N).GT.SMALL ) FCL = YL(N,1)/(SL(2,N)*PORD(2,N))
        IF( SG(2,N).GT.SMALL ) FCG = YG(N,1)/(SG(2,N)*PORD(2,N))
        FCLB = 0.D+0
        FCGB = 0.D+0
        IF( SLB(2,NB).GT.SMALL )FCLB = YLB(NB,1)/(SLB(2,NB)*PORDB(2,NB))
        IF( SGB(2,NB).GT.SMALL )FCGB = YGB(NB,1)/(SGB(2,NB)*PORDB(2,NB))
!
!---  Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGBX,VGBX,WGBX,N,M )
            CALL SHDP( WGBX,VGBX,UGBX,DISPL(IZN),DISPT(IZN),DPGB )
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULBX,VLBX,WLBX,N,M )
            CALL SHDP( WLBX,VLBX,ULBX,DISPL(IZN),DISPT(IZN),DPLB )
          ELSE
            DPGB = 0.D+0
            DPLB = 0.D+0
          ENDIF
!
!---  Dirichlet Boundary  ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &        IBCT(IEQC,NB).EQ.9 ) THEN
            K = KD(N)
            NPZ = NSZ(N)
            DGB = TORGB(2,NB)*PORDB(2,NB)*SGB(2,NB)*DFGOB(2,NB)
            DCG = CB(NB,1)*FCGB - C(N,1)*FCG
            INDX = 16
            DGZ = DIFMN( DGB,DGP,DZGF(N),DZGF(N),DCG,INDX)
            DGZ = (DGZ+DPGB)/(5.D-1*DZGF(N))
            AG = MAX( WG(1,NPZ),ZERO ) +
     &       DGZ*MAX((ONE-(TENTH*ABS(WG(1,NPZ))/(DGZ+SMALL)))**5,ZERO)
            AGP = MAX( -WG(1,NPZ),ZERO ) +
     &       DGZ*MAX((ONE-(TENTH*ABS(WG(1,NPZ))/(DGZ+SMALL)))**5,ZERO)
            DLB = TORLB(2,NB)*PORDB(2,NB)*SLB(2,NB)*DFLOB(2,NB)
            DCL = CB(NB,1)*FCLB - C(N,1)*FCL
            INDX = 16
            DLZ = DIFMN( DLB,DLP,DZGF(N),DZGF(N),DCL,INDX)
            DLZ = (DLZ+DPLB)/(5.D-1*DZGF(N))
            AL = MAX( WL(1,NPZ),ZERO ) +
     &       DLZ*MAX((ONE-(TENTH*ABS(WL(1,NPZ))/(DLZ+SMALL)))**5,ZERO)
            ALP = MAX( -WL(1,NPZ),ZERO ) +
     &       DLZ*MAX((ONE-(TENTH*ABS(WL(1,NPZ))/(DLZ+SMALL)))**5,ZERO)
            WC(NPZ,1) = AFZ(NPZ)*
     &     ( CB(NB,1)*(AL*FCLB + AG*FCGB) - C(N,1)*(ALP*FCL + AGP*FCG) )
!
!---  Outflow Boundary  ---
!
          ELSEIF( IBCT(IEQC,NB).EQ.7 ) THEN
            NPZ = NSZ(N)
            AG = MAX( WG(1,NPZ),ZERO )
            AGP = MAX( -WG(1,NPZ),ZERO )
            AL = MAX( WL(1,NPZ),ZERO )
            ALP = MAX( -WL(1,NPZ),ZERO )
            WC(NPZ,1) = AFZ(NPZ)*
     &     ( CB(NB,1)*(AL*FCLB + AG*FCGB) - C(N,1)*(ALP*FCL + AGP*FCG) )
          ENDIF
!
!---  South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGSX,VGSX,WGSX,N,M )
            CALL SHDP( VGSX,UGSX,WGSX,DISPL(IZN),DISPT(IZN),DPGS )
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULSX,VLSX,WLSX,N,M )
            CALL SHDP( VLSX,ULSX,WLSX,DISPL(IZN),DISPT(IZN),DPLS )
          ELSE
            DPGS = 0.D+0
            DPLS = 0.D+0
          ENDIF
!
!---  Dirichlet Boundary  ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &        IBCT(IEQC,NB).EQ.9 ) THEN
            J = JD(N)
            I = ID(N)
            NPY = NSY(N)
            DGB = TORGB(2,NB)*PORDB(2,NB)*SGB(2,NB)*DFGOB(2,NB)
            DCG = CB(NB,1)*FCGB - C(N,1)*FCG
            INDX = 16
            DGY = DIFMN(DGB,DGP,DYGF(N),DYGF(N),DCG,INDX)
            DGY = (DGY+DPGS)/((5.D-1*DYGF(N))*RP(I))
            AG = MAX( VG(1,NPY),ZERO ) +
     &       DGY*MAX((ONE-(TENTH*ABS(VG(1,NPY))/(DGY+SMALL)))**5,ZERO)
            AGP = MAX( -VG(1,NPY),ZERO ) +
     &       DGY*MAX((ONE-(TENTH*ABS(VG(1,NPY))/(DGY+SMALL)))**5,ZERO)
            DLB = TORLB(2,NB)*PORDB(2,NB)*SLB(2,NB)*DFLOB(2,NB)
            DCL = CB(NB,1)*FCLB - C(N,1)*FCL
            INDX = 16
            DLY = DIFMN(DLB,DLP,DYGF(N),DYGF(N),DCL,INDX)
            DLY = (DLY+DPLS)/((5.D-1*DYGF(N))*RP(I))
            AL = MAX( VL(1,NPY),ZERO ) +
     &       DLY*MAX((ONE-(TENTH*ABS(VL(1,NPY))/(DLY+SMALL)))**5,ZERO)
            ALP = MAX( -VL(1,NPY),ZERO ) +
     &       DLY*MAX((ONE-(TENTH*ABS(VL(1,NPY))/(DLY+SMALL)))**5,ZERO)
            VC(NPY,1) = AFY(NPY)*
     &     ( CB(NB,1)*(AL*FCLB + AG*FCGB) - C(N,1)*(ALP*FCL + AGP*FCG) )
!
!---  Outflow Boundary  ---
!
          ELSEIF( IBCT(IEQC,NB).EQ.7 ) THEN
            NPY = NSY(N)
            AG = MAX( VG(1,NPY),ZERO )
            AGP = MAX( -VG(1,NPY),ZERO )
            AL = MAX( VL(1,NPY),ZERO )
            ALP = MAX( -VL(1,NPY),ZERO )
            VC(NPY,1) = AFY(NPY)*
     &     ( CB(NB,1)*(AL*FCLB + AG*FCGB) - C(N,1)*(ALP*FCL + AGP*FCG) )
          ENDIF
!
!---  West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGX,VGX,WGX,N,M )
            CALL SHDP( UGX,VGX,WGX,DISPL(IZN),DISPT(IZN),DPGW )
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULX,VLX,WLX,N,M )
            CALL SHDP( ULX,VLX,WLX,DISPL(IZN),DISPT(IZN),DPLW )
          ELSE
            DPGW = 0.D+0
            DPLW = 0.D+0
          ENDIF
!
!---  Dirichlet Boundary  ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &        IBCT(IEQC,NB).EQ.9 ) THEN
            I = ID(N)
            NPX = NSX(N)
            DGB = TORGB(2,NB)*PORDB(2,NB)*SGB(2,NB)*DFGOB(2,NB)
            DCG = CB(NB,1)*FCGB - C(N,1)*FCG
            INDX = 16
            DGX = DIFMN( DGB,DGP,DXGF(N),DXGF(N),DCG,INDX)
            DGX = (DGX+DPGW)/(5.D-1*DXGF(N))
            AG = MAX( UG(1,NPX),ZERO ) +
     &       DGX*MAX((ONE-(TENTH*ABS(UG(1,NPX))/(DGX+SMALL)))**5,ZERO)
            AGP = MAX( -UG(1,NPX),ZERO ) +
     &       DGX*MAX((ONE-(TENTH*ABS(UG(1,NPX))/(DGX+SMALL)))**5,ZERO)
            DLB = TORLB(2,NB)*PORDB(2,NB)*SLB(2,NB)*DFLOB(2,NB)
            DCL = CB(NB,1)*FCLB - C(N,1)*FCL
            INDX = 16
            DLX = DIFMN( DLB,DLP,DXGF(N),DXGF(N),DCL,INDX)
            DLX = (DLX+DPLW)/(5.D-1*DXGF(N))
            AL = MAX( UL(1,NPX),ZERO ) +
     &       DLX*MAX((ONE-(TENTH*ABS(UL(1,NPX))/(DLX+SMALL)))**5,ZERO)
            ALP = MAX( -UL(1,NPX),ZERO ) +
     &       DLX*MAX((ONE-(TENTH*ABS(UL(1,NPX))/(DLX+SMALL)))**5,ZERO)
            UC(NPX,1) = AFX(NPX)*
     &     ( CB(NB,1)*(AL*FCLB + AG*FCGB) - C(N,1)*(ALP*FCL + AGP*FCG) )
!
!---  Outflow Boundary  ---
!
          ELSEIF( IBCT(IEQC,NB).EQ.7 ) THEN
            NPX = NSX(N)
            AG = MAX( UG(1,NPX),ZERO )
            AGP = MAX( -UG(1,NPX),ZERO )
            AL = MAX( UL(1,NPX),ZERO )
            ALP = MAX( -UL(1,NPX),ZERO )
            UC(NPX,1) = AFX(NPX)*
     &     ( CB(NB,1)*(AL*FCLB + AG*FCGB) - C(N,1)*(ALP*FCL + AGP*FCG) )
          ENDIF
!
!---  East boundary
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGEX,VGEX,WGEX,N,M )
            CALL SHDP( UGEX,VGEX,WGEX,DISPL(IZN),DISPT(IZN),DPGE )
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULEX,VLEX,WLEX,N,M )
            CALL SHDP( ULEX,VLEX,WLEX,DISPL(IZN),DISPT(IZN),DPLE )
          ELSE
            DPGE = 0.D+0
            DPLE = 0.D+0
          ENDIF
!
!---  Dirichlet Boundary  ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &        IBCT(IEQC,NB).EQ.9 ) THEN
            I = ID(N)
            NQX = NSX(N) + 1
            DGB = TORGB(2,NB)*PORDB(2,NB)*SGB(2,NB)*DFGOB(2,NB)
            DCG = C(N,1)*FCG - CB(NB,1)*FCGB
            INDX = 16
            DGX = DIFMN(DGP,DGB,DXGF(N),DXGF(N),DCG,INDX)
            DGX = (DGX+DPGE)/(5.D-1*DXGF(N))
            AG = MAX( -UG(1,NQX),ZERO ) +
     &       DGX*MAX((ONE-(TENTH*ABS(UG(1,NQX))/(DGX+SMALL)))**5,ZERO)
            AGP = MAX( UG(1,NQX),ZERO ) +
     &       DGX*MAX((ONE-(TENTH*ABS(UG(1,NQX))/(DGX+SMALL)))**5,ZERO)
            DLB = TORLB(2,NB)*PORDB(2,NB)*SLB(2,NB)*DFLOB(2,NB)
            DCL = C(N,1)*FCL - CB(NB,1)*FCLB
            INDX = 16
            DLX = DIFMN(DLP,DLB,DXGF(N),DXGF(N),DCL,INDX)
            DLX = (DLX+DPLE)/(5.D-1*DXGF(N))
            AL = MAX( -UL(1,NQX),ZERO ) +
     &       DLX*MAX((ONE-(TENTH*ABS(UL(1,NQX))/(DLX+SMALL)))**5,ZERO)
            ALP = MAX( UL(1,NQX),ZERO ) +
     &       DLX*MAX((ONE-(TENTH*ABS(UL(1,NQX))/(DLX+SMALL)))**5,ZERO)
            UC(NQX,1) = AFX(NQX)*
     &     ( C(N,1)*(ALP*FCL + AGP*FCG) - CB(NB,1)*(AL*FCLB + AG*FCGB) )
!
!---  Outflow Boundary  ---
!
          ELSEIF( IBCT(IEQC,NB).EQ.7 ) THEN
            NQX = NSX(N) + 1
            AG = MAX( -UG(1,NQX),ZERO )
            AGP = MAX( UG(1,NQX),ZERO )
            AL = MAX( -UL(1,NQX),ZERO )
            ALP = MAX( UL(1,NQX),ZERO )
            UC(NQX,1) = AFX(NQX)*
     &     ( C(N,1)*(ALP*FCL + AGP*FCG) - CB(NB,1)*(AL*FCLB + AG*FCGB) )
          ENDIF
!
!---  North boundary
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGNX,VGNX,WGNX,N,M )
            CALL SHDP( VGNX,UGNX,WGNX,DISPL(IZN),DISPT(IZN),DPGN )
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULNX,VLNX,WLNX,N,M )
            CALL SHDP( VLNX,ULNX,WLNX,DISPL(IZN),DISPT(IZN),DPLN )
          ELSE
            DPGN = 0.D+0
            DPLN = 0.D+0
          ENDIF
!
!---  Dirichlet Boundary  ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &        IBCT(IEQC,NB).EQ.9 ) THEN
            J = JD(N)
            I = ID(N)
            NQY = NSY(N) + IFLD
            DGB = TORGB(2,NB)*PORDB(2,NB)*SGB(2,NB)*DFGOB(2,NB)
            DCG = C(N,1)*FCG - CB(NB,1)*FCGB
            INDX = 16
            DGY = DIFMN(DGB,DGP,DYGF(N),DYGF(N),DCG,INDX)
            DGY = (DGY+DPGN)/(RP(I)*(5.D-1*DYGF(N)))
            AG = MAX( -VG(1,NQY),ZERO ) +
     &       DGY*MAX((ONE-(TENTH*ABS(VG(1,NQY))/(DGY+SMALL)))**5,ZERO)
            AGP = MAX( VG(1,NQY),ZERO ) +
     &       DGY*MAX((ONE-(TENTH*ABS(VG(1,NQY))/(DGY+SMALL)))**5,ZERO)
            DLB = TORLB(2,NB)*PORDB(2,NB)*SLB(2,NB)*DFLOB(2,NB)
            DCL = C(N,1)*FCL - CB(NB,1)*FCLB
            INDX = 16
            DLY = DIFMN(DLB,DLP,DYGF(N),DYGF(N),DCL,INDX)
            DLY = (DLY+DPLN)/(RP(I)*(5.D-1*DYGF(N)))
            AL = MAX( -VL(1,NQY),ZERO ) +
     &       DLY*MAX((ONE-(TENTH*ABS(VL(1,NQY))/(DLY+SMALL)))**5,ZERO)
            ALP = MAX( VL(1,NQY),ZERO ) +
     &       DLY*MAX((ONE-(TENTH*ABS(VL(1,NQY))/(DLY+SMALL)))**5,ZERO)
            VC(NQY,1) = AFY(NQY)*
     &     ( C(N,1)*(ALP*FCL + AGP*FCG) - CB(NB,1)*(AL*FCLB + AG*FCGB) )
!
!---  Outflow Boundary  ---
!
          ELSEIF( IBCT(IEQC,NB).EQ.7 ) THEN
            NQY = NSY(N) + IFLD
            AG = MAX( -VG(1,NQY),ZERO )
            AGP = MAX( VG(1,NQY),ZERO )
            AL = MAX( -VL(1,NQY),ZERO )
            ALP = MAX( VL(1,NQY),ZERO )
            VC(NQY,1) = AFY(NQY)*
     &     ( C(N,1)*(ALP*FCL + AGP*FCG) - CB(NB,1)*(AL*FCLB + AG*FCGB) )
          ENDIF
!
!---  Top boundary
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGTX,VGTX,WGTX,N,M )
            CALL SHDP( WGTX,VGTX,UGTX,DISPL(IZN),DISPT(IZN),DPGT )
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULTX,VLTX,WLTX,N,M )
            CALL SHDP( WLTX,VLTX,ULTX,DISPL(IZN),DISPT(IZN),DPLT )
          ELSE
            DPGT = 0.D+0
            DPLT = 0.D+0
          ENDIF
!
!---  Dirichlet Boundary  ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &        IBCT(IEQC,NB).EQ.9 ) THEN
            K = KD(N)
            NQZ = NSZ(N) + IJFLD
            DGB = TORGB(2,NB)*PORDB(2,NB)*SGB(2,NB)*DFGOB(2,NB)
            DCG = C(N,1)*FCG - CB(NB,1)*FCGB
            INDX = 16
            DGZ = DIFMN(DGB,DGP,DZGF(N),DZGF(N),DCG,INDX)
            DGZ = (DGZ+DPGT)/(5.D-1*DZGF(N))
            AG = MAX( -WG(1,NQZ),ZERO ) +
     &       DGZ*MAX((ONE-(TENTH*ABS(WG(1,NQZ))/(DGZ+SMALL)))**5,ZERO)
            AGP = MAX( WG(1,NQZ),ZERO ) +
     &       DGZ*MAX((ONE-(TENTH*ABS(WG(1,NQZ))/(DGZ+SMALL)))**5,ZERO)
            DLB = TORLB(2,NB)*PORDB(2,NB)*SLB(2,NB)*DFLOB(2,NB)
            DCL = C(N,1)*FCL - CB(NB,1)*FCLB
            INDX = 16
            DLZ = DIFMN(DLB,DLP,DZGF(N),DZGF(N),DCL,INDX)
            DLZ = (DLZ+DPLT)/(5.D-1*DZGF(N))
            AL = MAX( -WL(1,NQZ),ZERO ) +
     &       DLZ*MAX((ONE-(TENTH*ABS(WL(1,NQZ))/(DLZ+SMALL)))**5,ZERO)
            ALP = MAX( WL(1,NQZ),ZERO ) +
     &       DLZ*MAX((ONE-(TENTH*ABS(WL(1,NQZ))/(DLZ+SMALL)))**5,ZERO)
            WC(NQZ,1) = AFZ(NQZ)*
     &     ( C(N,1)*(ALP*FCL + AGP*FCG) - CB(NB,1)*(AL*FCLB + AG*FCGB) )
!
!---  Outflow Boundary  ---
!
          ELSEIF( IBCT(IEQC,NB).EQ.7 ) THEN
            NQZ = NSZ(N) + IJFLD
            AG = MAX( -WG(1,NQZ),ZERO )
            AGP = MAX( WG(1,NQZ),ZERO )
            AL = MAX( -WL(1,NQZ),ZERO )
            ALP = MAX( WL(1,NQZ),ZERO )
            WC(NQZ,1) = AFZ(NQZ)*
     &     ( C(N,1)*(ALP*FCL + AGP*FCG) - CB(NB,1)*(AL*FCLB + AG*FCGB) )
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VFLGB group  ---
!
      RETURN
      END


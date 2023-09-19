!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_BOREHOLE( INDX )
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
!
!     Define borehole nodes, determine trajectory points, and 
!     check for borehole trajectories within node surface planes
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 5 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_BH
      USE OUTPU
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
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
      REAL*8 XPX(5),YPX(5),ZPX(5)
      REAL*8 XIX(2),YIX(2),ZIX(2)
      REAL*8 AJ(3,3),BJ(3),VDX(3)
      REAL*8 XEX(4),YEX(4),ZEX(4)
      REAL*8, DIMENSION(3) :: PX
      INTEGER MSX(4,6),IJ(3)
      INTEGER N1X(4),N2X(4),NODX(8)
      CHARACTER*28 FORM1
      CHARACTER*33 FORM2
!
!----------------------Data Statements---------------------------------!
!
      DATA MSX / 1,2,4,3,1,5,6,2,1,3,7,5,2,6,8,4,3,4,8,7,5,7,8,6 /
      DATA N1X / 2,3,4,1 /
      DATA N2X / 1,2,3,4 /
      DATA FORM1 /'(A,I3,A,I3,A,I3,A,I3,A,I8,A)'/
      DATA FORM2 /'(A,I1,A,I1,A,I1,A,I1,A,I1,A,I1,A)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHK_BOREHOLE'
      EPSLX = 1.D-6
!
!---  Skip checks when reading preprocessed data  ---
!
      IF( INDX.EQ.0 ) THEN
!
!---  Loop over boreholes, checking for intersections of the
!     borehole trajectory with nodes  ---
!
      DO 600 NBH = 1,N_BH
        ID_BH(3,NBH) = NBN_BH+1
!
!---    Loop over number of borehole intervals  ---
!
        DO 490 NICW = ID_BH(1,NBH),ID_BH(2,NBH)
!
!---    Loop over active nodes to find borehole nodes and borehole
!       projections ---
!
          DO 480 N = 1,NFLD
            IF( IXP(N).EQ.0 ) GOTO 480
            I = ID(N)
            J = JD(N)
            K = KD(N)
            NC = 0
!
!---        Determine whether the transition points are within
!           the hexahedron node volume or on the hexahedron
!           surfaces  ---
!
            DO 190 NPT = 1,2
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1
     &          .AND. I.EQ.1 ) THEN
                XPX(1) = 0.D+0
                XPX(2) = 0.D+0
                YPX(1) = 0.D+0
                YPX(2) = 0.D+0
                ZPX(1) = ZE(1,N)
                ZPX(2) = ZE(5,N)
!
!---            Node height greater than EPSLX  ---
!
                IF( ABS(ZPX(1)-ZPX(2)).GT.EPSLX ) THEN
                  DZPX1 = ZTP_BH(NPT,NICW)-ZPX(1)
                  DZPX2 = ZPX(2)-ZTP_BH(NPT,NICW)
                  IF( ABS(DZPX1).LT.EPSLX ) DZPX1 = 0.D+0
                  IF( ABS(DZPX2).LT.EPSLX ) DZPX2 = 0.D+0
!
!---              Transition point within vertical limits of node  ---
!
                  IF( DZPX1.GE.0.D+0 .AND. DZPX2.GE.0.D+0 ) THEN
                    NC = NC+1
                    XIX(NC) = 0.D+0
                    YIX(NC) = 0.D+0
                    ZIX(NC) = ZTP_BH(NPT,NICW)
                  ENDIF
                ENDIF
                GOTO 190
              ENDIF
!
!---          Check for point with hexahedron  ---
!
              XMNX = XE(1,N)
              XMAX = XE(1,N)
              YMNX = YE(1,N)
              YMAX = YE(1,N)
              ZMNX = ZE(1,N)
              ZMAX = ZE(1,N)
              DO M = 2,8
                XMNX = MIN( XMNX,XE(M,N) )
                XMAX = MAX( XMAX,XE(M,N) )
                YMNX = MIN( YMNX,YE(M,N) )
                YMAX = MAX( YMAX,YE(M,N) )
                ZMNX = MIN( ZMNX,ZE(M,N) )
                ZMAX = MAX( ZMAX,ZE(M,N) )
              ENDDO
              XMNX = XMNX - 1.D-6
              XMAX = XMAX + 1.D-6
              YMNX = YMNX - 1.D-6
              YMAX = YMAX + 1.D-6
              ZMNX = ZMNX - 1.D-6
              ZMAX = ZMAX + 1.D-6
              ICWX = 0
              IF( XTP_BH(NPT,NICW).GE.XMNX .AND. 
     &          XTP_BH(NPT,NICW).LE.XMAX .AND. 
     &          YTP_BH(NPT,NICW).GE.YMNX .AND. 
     &          YTP_BH(NPT,NICW).LE.YMAX .AND. 
     &          ZTP_BH(NPT,NICW).GE.ZMNX .AND. 
     &          ZTP_BH(NPT,NICW).LE.ZMAX ) 
     &          CALL WITHIN( XTP_BH(NPT,NICW),YTP_BH(NPT,NICW),
     &            ZTP_BH(NPT,NICW),ICWX,N )
!
!---          Opposing rotations found or borehole point not within
!             hexahedron limits, point outside hexahedron  ---
!
              IF( ICWX.EQ.0 ) GOTO 190
!
!---          No opposing rotations found, point inside hexahedron  ---
!
              NC = NC+1
              XIX(NC) = XTP_BH(NPT,NICW)
              YIX(NC) = YTP_BH(NPT,NICW)
              ZIX(NC) = ZTP_BH(NPT,NICW)
  190       CONTINUE
!
!---        Both transition points inside hexahedron, skip
!           search for borehole path crossing hexahedron surfaces  ---
!
            IF( NC.EQ.2 ) GOTO 232
!
!---        Cylindrical coordinates with azimuthal symmetry,
!           centrally located wells  ---
!
            IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1
     &        .AND. I.EQ.1 ) THEN
!
!---          Interval crosses lower node surface  ---
!
              DZPX1 = ZPX(1)-ZTP_BH(1,NICW)
              IF( ABS(DZPX1).LT.EPSLX ) DZPX1 = 0.D+0
              DZPX2 = ZPX(1)-ZTP_BH(2,NICW)
              IF( ABS(DZPX2).LT.EPSLX ) DZPX2 = 0.D+0
              IF( (DZPX1*DZPX2).LT.-EPSLX ) THEN
                NC = NC+1
                XIX(NC) = 0.D+0
                YIX(NC) = 0.D+0
                ZIX(NC) = ZPX(1)
              ENDIF
!
!---          Interval crosses upper node surface  ---
!
              DZPX1 = ZPX(2)-ZTP_BH(1,NICW)
              IF( ABS(DZPX1).LT.EPSLX ) DZPX1 = 0.D+0
              DZPX2 = ZPX(2)-ZTP_BH(2,NICW)
              IF( ABS(DZPX2).LT.EPSLX ) DZPX2 = 0.D+0
              IF( (DZPX1*DZPX2).LT.-EPSLX ) THEN
                NC = NC+1
                XIX(NC) = 0.D+0
                YIX(NC) = 0.D+0
                ZIX(NC) = ZPX(2)
              ENDIF
              GOTO 232
            ENDIF
!
!---        Loop over node surfaces,
!           (bottom,south,west,east,north,top)  ---
!
            DO 230 NS = 1,6
!
!---          Define the five surface points, four corners
!             and one centroid---
!
              DO 200 NP = 1,4
                MX = MSX(NP,NS)
!
!---            Cylindrical coordinates---
!
                IF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
                  XPX(NP) = XE(MX,N)*COS(YE(MX,N))
                  YPX(NP) = XE(MX,N)*SIN(YE(MX,N))
                  ZPX(NP) = ZE(MX,N)
                ELSE
                  XPX(NP) = XE(MX,N)
                  YPX(NP) = YE(MX,N)
                  ZPX(NP) = ZE(MX,N)
                ENDIF
  200         CONTINUE
              NP = 4
              CALL PGCNTRD( NP,XPX(1),YPX(1),ZPX(1),
     &          XPX(5),YPX(5),ZPX(5) )
!
!
!---          Loop over the four triangular planes on the 
!             surface face  ---
!
              DO 220 NT = 1,4
!
!---            Load plane-line intersection matrix and problem vector  ---
!
                AJ(1,1) = XTP_BH(1,NICW)-XTP_BH(2,NICW)
                AJ(2,1) = YTP_BH(1,NICW)-YTP_BH(2,NICW)
                AJ(3,1) = ZTP_BH(1,NICW)-ZTP_BH(2,NICW)
                AJ(1,2) = XPX(N1X(NT))-XPX(5)
                AJ(2,2) = YPX(N1X(NT))-YPX(5)
                AJ(3,2) = ZPX(N1X(NT))-ZPX(5)
                AJ(1,3) = XPX(N2X(NT))-XPX(5)
                AJ(2,3) = YPX(N2X(NT))-YPX(5)
                AJ(3,3) = ZPX(N2X(NT))-ZPX(5)
                BJ(1) = XTP_BH(1,NICW)-XPX(5)
                BJ(2) = YTP_BH(1,NICW)-YPX(5)
                BJ(3) = ZTP_BH(1,NICW)-ZPX(5)
!
!---            Check for no intersection  ---
!
                DO 210 IP = 1,3
                  AAMAX = 0.D+0
                  DO 202 JP = 1,3
                    IF( ABS(AJ(IP,JP)).GT.AAMAX ) AAMAX = ABS(AJ(IP,JP))
  202             CONTINUE
!
!---              No intersection go to next triangle on the 
!                 surface  ---
!
                  IF( ABS(AAMAX)/EPSL.LT.EPSL ) GOTO 220
  210           CONTINUE
!
!---            Find plane-line intersection matrix inverse  ---
!
                JP = 3
                KP = 3
                CALL LUDCMP( AJ,JP,KP,IJ,DJ )
                CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---            Find plane-line intersection point  ---
!
                TX = BJ(1)
                UX = BJ(2)
                VX = BJ(3)
                IF( ABS(TX).LT.EPSL ) TX = 0.D+0
                IF( ABS(UX).LT.EPSL ) UX = 0.D+0
                IF( ABS(VX).LT.EPSL ) VX = 0.D+0
!
!---            Line crosses surface, within the triangle  ---
!
                IF( TX.GE.0.D+0 .AND. TX.LE.1.D+0 .AND.
     &            UX.GE.0.D+0 .AND. UX.LE.1.D+0 .AND. 
     &            VX.GE.0.D+0 .AND. VX.LE.1.D+0 .AND.
     &            (UX+VX).LE.1.D+0 ) THEN
                  XTX = XTP_BH(1,NICW) 
     &              + (XTP_BH(2,NICW)-XTP_BH(1,NICW))*TX
                  YTX = YTP_BH(1,NICW)
     &              + (YTP_BH(2,NICW)-YTP_BH(1,NICW))*TX
                  ZTX = ZTP_BH(1,NICW)
     &              + (ZTP_BH(2,NICW)-ZTP_BH(1,NICW))*TX
!
!---              Check for non-distinct points  ---
!
                  IF( NC.GE.1 ) THEN
                    DO 212 NDP = 1,NC
                      DPX = SQRT( ((XTX-XIX(NDP))**2) 
     &                  + ((YTX-YIX(NDP))**2) + ((ZTX-ZIX(NDP))**2) )
!
!---                  Duplicate point found  ---
!
                      IF( DPX.LT.EPSLX ) GOTO 214
  212               CONTINUE
                  ENDIF
                  IF( NC.EQ.2 ) THEN
                    INDX = 7
                    CHMSG = 'Three Distinct Coupled Well Points'
     &                //  ' at Node'
                    IMSG = N
                    CALL WRMSGS( INDX )
                  ENDIF
                  NC = NC + 1
                  XIX(NC) = XTX
                  YIX(NC) = YTX
                  ZIX(NC) = ZTX
                ENDIF
  214           CONTINUE
  220         CONTINUE
  230       CONTINUE
  232       CONTINUE
!
!---        Check that two borehole points are distinct  ---
!
            IF( NC.EQ.2 ) THEN
              DPX = SQRT( ((XIX(2)-XIX(1))**2) + ((YIX(2)-YIX(1))**2)
     &          + ((ZIX(2)-ZIX(1))**2) )
            ENDIF
!
!---        Two distinct borehole points within node define or update
!           borehole node trajectory points  ---
!
            IF( NC.EQ.2 .AND. DPX.GT.EPSLX ) THEN
!
!---          Cylindrical coordinates with azimuthal symmetry  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) GOTO 362
!
!---          Check if line between borehole points is contained in
!             a node surface plane, loop over node surfaces,
!             (bottom,south,west,east,north,top)  ---
!
              DO 360 NS = 1,6
!
!---            Define the five surface points, four corners
!               and one centroid---
!
                DO 352 NP = 1,4
                  MX = MSX(NP,NS)
!
!---              Cylindrical coordinates---
!
                  IF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
                    XPX(NP) = XE(MX,N)*COS(YE(MX,N))
                    YPX(NP) = XE(MX,N)*SIN(YE(MX,N))
                    ZPX(NP) = ZE(MX,N)
                  ELSE
                    XPX(NP) = XE(MX,N)
                    YPX(NP) = YE(MX,N)
                    ZPX(NP) = ZE(MX,N)
                  ENDIF
  352           CONTINUE
                NP = 4
                CALL PGCNTRD( NP,XPX(1),YPX(1),ZPX(1),
     &            XPX(5),YPX(5),ZPX(5) )
!
!---            Loop over the four triangular planes on the 
!               surface face  ---
!
                DO 358 NT = 1,4
!
!---              Loop over trajectory points  ---
!
                  DO 356 NPT = 1,2
!
!---                Load point-plane matrix  ---
!
                    AJ(1,1) = XIX(NPT)-XPX(5)
                    AJ(2,1) = XIX(NPT)-XPX(N1X(NT))
                    AJ(3,1) = XIX(NPT)-XPX(N2X(NT))
                    AJ(1,2) = YIX(NPT)-YPX(5)
                    AJ(2,2) = YIX(NPT)-YPX(N1X(NT))
                    AJ(3,2) = YIX(NPT)-YPX(N2X(NT))
                    AJ(1,3) = ZIX(NPT)-ZPX(5)
                    AJ(2,3) = ZIX(NPT)-ZPX(N1X(NT))
                    AJ(3,3) = ZIX(NPT)-ZPX(N2X(NT))
!
!---                Check for singular matrix  ---
!
                    DO 310 IP = 1,3
                      AAMAX = 0.D+0
                      DO 302 JP = 1,3
                        IF( ABS(AJ(IP,JP)).GT.AAMAX ) 
     &                  AAMAX = ABS(AJ(IP,JP))
  302                 CONTINUE
!
!---                  Singular matrix, trajectory point is within
!                     the surface plane  ---
!
                      IF( ABS(AAMAX)/EPSL.LT.EPSL ) GOTO 356
  310               CONTINUE
!
!---                Find matrix determinant  ---
!
                    JP = 3
                    KP = 3
                    CALL LUDCMP( AJ,JP,KP,IJ,DJ )
                    DO 354 M = 1,JP
                      DJ = DJ*AJ(M,M)
  354               CONTINUE
!
!---                If determinant equals zero trajectory point
!                   is within the surface plane  ---
!
                    IF( ABS(DJ).GT.EPSL ) GOTO 358
  356             CONTINUE
                  INDX = 7
                  CHMSG = 'Both Well Trajectory Points within '
     &              //  'Surface Plane at Node'
                  IMSG = N
                  CALL WRMSGS( INDX )
  358           CONTINUE
  360         CONTINUE
  362         CONTINUE
!
!---          Loop over borehole nodes in current borehole  ---
!
              DO 370 NWN = ID_BH(3,NBH),NBN_BH
!
!---            Well node previously counted  ---
!
                IF( IBN_BH(NWN).EQ.N ) THEN
                  GOTO 372
                ENDIF
  370         CONTINUE
  372         CONTINUE
!
!---          Define a new borehole node and set the trajectory
!             points  ---
!
              NBN_BH = NBN_BH + 1
              IF( NBN_BH.GT.LBN_BH ) THEN
               INDX = 5
               CHMSG = 'Number of Borehole Nodes ' //
     &           '> Parameter LBN_BH'
               CALL WRMSGS( INDX )
              ENDIF
              XP_BH(1,NBN_BH) = XIX(1)
              YP_BH(1,NBN_BH) = YIX(1)
              ZP_BH(1,NBN_BH) = ZIX(1)
              XP_BH(2,NBN_BH) = XIX(2)
              YP_BH(2,NBN_BH) = YIX(2)
              ZP_BH(2,NBN_BH) = ZIX(2)
              IBN_BH(NBN_BH) = N
              INV_BH(NBN_BH) = NICW
!
!---          Well projections on local grid  ---
!
              NWN = NBN_BH
              CALL PROJ_BH( PLX_BH(NWN),PLY_BH(NWN),PLZ_BH(NWN),
     &          XP_BH(1,NWN),YP_BH(1,NWN),ZP_BH(1,NWN),N )
!
!---          Double the number borehole nodes
!             for coaxial boreholes  ---
!
              IF( IT_BH(1,NBH).GE.10000 ) THEN
                NBN_BH = NBN_BH + 1
                IF( NBN_BH.GT.LBN_BH ) THEN
                 INDX = 5
                 CHMSG = 'Number of Borehole Nodes ' //
     &             '> Parameter LBN_BH'
                 CALL WRMSGS( INDX )
                ENDIF
                XP_BH(1,NBN_BH) = XIX(1)
                YP_BH(1,NBN_BH) = YIX(1)
                ZP_BH(1,NBN_BH) = ZIX(1)
                XP_BH(2,NBN_BH) = XIX(2)
                YP_BH(2,NBN_BH) = YIX(2)
                ZP_BH(2,NBN_BH) = ZIX(2)
                IBN_BH(NBN_BH) = N
                INV_BH(NBN_BH) = NICW
!
!---            Well projections on local grid  ---
!
                NWN = NBN_BH
                CALL PROJ_BH( PLX_BH(NWN),PLY_BH(NWN),PLZ_BH(NWN),
     &            XP_BH(1,NWN),YP_BH(1,NWN),ZP_BH(1,NWN),N )
              ENDIF
            ENDIF
  480     CONTINUE
  490   CONTINUE
        ID_BH(4,NBH) = NBN_BH
!
!---    Sequence borehole nodes according to their distance from
!       the previous borehole, beginning with the starting point 
!       of the borehole  ---
!
        I1X = ID_BH(1,NBH)
        I3X = ID_BH(3,NBH)
        I4X = ID_BH(4,NBH)
        XTP_BHX = XTP_BH(1,I1X)
        YTP_BHX = YTP_BH(1,I1X)
        ZTP_BHX = ZTP_BH(1,I1X)
  500   CONTINUE
        DMNX = 1.D+20
        DO 510 NBN = I3X,I4X
          XWPX = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
          YWPX = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
          ZWPX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
          DISTX = SQRT( (XTP_BHX-XWPX)**2 +
     &      (YTP_BHX-YWPX)**2 + (ZTP_BHX-ZWPX)**2 )
          IF( DISTX.LT.DMNX ) THEN
            DMNX = DISTX
            IMNX = NBN
            IWNX = IBN_BH(NBN)
            INVX = INV_BH(NBN)
            XIX(1) = XP_BH(1,NBN)
            XIX(2) = XP_BH(2,NBN)
            YIX(1) = YP_BH(1,NBN)
            YIX(2) = YP_BH(2,NBN)
            ZIX(1) = ZP_BH(1,NBN)
            ZIX(2) = ZP_BH(2,NBN)
            PLX_BHX = PLX_BH(NBN)
            PLY_BHX = PLY_BH(NBN)
            PLZ_BHX = PLZ_BH(NBN)
          ENDIF
  510   CONTINUE
        DO 520 JCW = I4X,I3X,-1
          IF( JCW.LT.IMNX ) THEN
            IBN_BH(JCW+1) = IBN_BH(JCW)
            INV_BH(JCW+1) = INV_BH(JCW)
            XP_BH(1,JCW+1) = XP_BH(1,JCW)
            XP_BH(2,JCW+1) = XP_BH(2,JCW)
            YP_BH(1,JCW+1) = YP_BH(1,JCW)
            YP_BH(2,JCW+1) = YP_BH(2,JCW)
            ZP_BH(1,JCW+1) = ZP_BH(1,JCW)
            ZP_BH(2,JCW+1) = ZP_BH(2,JCW)
            PLX_BH(JCW+1) = PLX_BH(JCW)
            PLY_BH(JCW+1) = PLY_BH(JCW)
            PLZ_BH(JCW+1) = PLZ_BH(JCW)
          ENDIF
  520   CONTINUE
        IBN_BH(I3X) = IWNX
        INV_BH(I3X) = INVX
        XP_BH(1,I3X) = XIX(1)
        XP_BH(2,I3X) = XIX(2)
        YP_BH(1,I3X) = YIX(1)
        YP_BH(2,I3X) = YIX(2)
        ZP_BH(1,I3X) = ZIX(1)
        ZP_BH(2,I3X) = ZIX(2)
        PLX_BH(I3X) = PLX_BHX
        PLY_BH(I3X) = PLY_BHX
        PLZ_BH(I3X) = PLZ_BHX
        XTP_BHX = XP_BH(2,I3X)
        YTP_BHX = YP_BH(2,I3X)
        ZTP_BHX = ZP_BH(2,I3X)
        I3X = I3X+1
        IF( I3X.LT.I4X ) GOTO 500
!
!---    Sequence borehole nodes according to their distance from
!       the previous borehole, beginning with the starting point 
!       of the borehole  ---
!
        I1X = ID_BH(1,NBH)
        I3X = ID_BH(3,NBH)
        I4X = ID_BH(4,NBH)
        XTP_BHX = XTP_BH(1,I1X)
        YTP_BHX = YTP_BH(1,I1X)
        ZTP_BHX = ZTP_BH(1,I1X)
        DO 530 NBN = I3X,I4X
          DIST1X = SQRT( (XTP_BHX-XP_BH(1,NBN))**2 +
     &      (YTP_BHX-YP_BH(1,NBN))**2 + 
     &      (ZTP_BHX-ZP_BH(1,NBN))**2 )
          DIST2X = SQRT( (XTP_BHX-XP_BH(2,NBN))**2 +
     &      (YTP_BHX-YP_BH(2,NBN))**2 + 
     &      (ZTP_BHX-ZP_BH(2,NBN))**2 )
          IF( DIST1X.GT.DIST2X ) THEN
            XIX(2) = XP_BH(1,NBN)
            YIX(2) = YP_BH(1,NBN)
            ZIX(2) = ZP_BH(1,NBN)
            XP_BH(1,NBN) = XP_BH(2,NBN)
            YP_BH(1,NBN) = YP_BH(2,NBN)
            ZP_BH(1,NBN) = ZP_BH(2,NBN)
            XP_BH(2,NBN) = XIX(2)
            YP_BH(2,NBN) = YIX(2)
            ZP_BH(2,NBN) = ZIX(2)
          ENDIF
          XTP_BHX = XP_BH(2,NBN)
          YTP_BHX = YP_BH(2,NBN)
          ZTP_BHX = ZP_BH(2,NBN)       
  530   CONTINUE
  600 CONTINUE
!
!---  Loop over boreholes, re-numbering borehole nodes for coaxial
!     boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Coaxial borehole  ---
!
        IF( IT_BH(1,NBH).GE.10000 ) THEN
          I1X = ID_BH(3,NBH)
          I2X = (ID_BH(4,NBH)-ID_BH(3,NBH))/2 + ID_BH(3,NBH)
          DO NBN = I1X,I2X
            NBN1 = (NBN-ID_BH(3,NBH))*2 + 1
            XP_BH(1,NBN) = XP_BH(1,NBN1)
            XP_BH(2,NBN) = XP_BH(2,NBN1)
            YP_BH(1,NBN) = YP_BH(1,NBN1)
            YP_BH(2,NBN) = YP_BH(2,NBN1)
            ZP_BH(1,NBN) = ZP_BH(1,NBN1)
            ZP_BH(2,NBN) = ZP_BH(2,NBN1)
            INV_BH(NBN) = INV_BH(NBN1)
            IBN_BH(NBN) = IBN_BH(NBN1)
          ENDDO
          DO NBN = I1X,I2X
            NBN1 = ID_BH(4,NBH) - NBN + 1
            XP_BH(1,NBN1) = XP_BH(2,NBN)
            XP_BH(2,NBN1) = XP_BH(1,NBN)
            YP_BH(1,NBN1) = YP_BH(2,NBN)
            YP_BH(2,NBN1) = YP_BH(1,NBN)
            ZP_BH(1,NBN1) = ZP_BH(2,NBN)
            ZP_BH(2,NBN1) = ZP_BH(1,NBN)
            INV_BH(NBN1) = INV_BH(NBN)
            IBN_BH(NBN1) = IBN_BH(NBN)
            IBN_BH(NBN) = 0
          ENDDO
        ENDIF
      ENDDO
!
!---  Closing ENDIF for preprocessed borehole data  ---
!
      ENDIF
!
!---  Record boreholes nodes to output file ---
!
      WRITE(IWR,'(//,A,/)') ' --- Borehole Node Record  ---'
!
!---  Loop over boreholes  ---
!
      DO 720 NBH = 1,N_BH
        NCH = INDEX( NM_BH(NBH),'  ') - 1
        WRITE(IWR,'(2X,A,I4,2A)') 'Borehole Node and Field Node ' //
     &    'Connections: ',NBH,' Borehole Name: ',NM_BH(NBH)(1:NCH)
!
!---    Loop over borehole nodes  ---
!
        DO 710 NBN = ID_BH(3,NBH),ID_BH(4,NBH)
          NX = IBN_BH(NBN)
          IF( IT_BH(1,NBH).GE.10000 .AND. NX.EQ.0 ) THEN
            NX = IBN_BH(ID_BH(4,NBH)-NBN+1)
            IX = ID(NX)
            JX = JD(NX)
            KX = KD(NX)
            WRITE(IWR,FORM1) '    Inner Coaxial Borehole Node: ',NBN,
     &        ' <=> Field Node: (',IX,',',JX,',',KX,':',NX,')'
          ELSEIF( IT_BH(1,NBH).GE.10000 .AND. NX.NE.0 ) THEN
            IX = ID(NX)
            JX = JD(NX)
            KX = KD(NX)
            WRITE(IWR,FORM1) '    Outer Coaxial Borehole Node: ',NBN,
     &        ' <=> Field Node: (',IX,',',JX,',',KX,':',NX,')'
          ELSE
            IX = ID(NX)
            JX = JD(NX)
            KX = KD(NX)
            WRITE(IWR,FORM1) '    Borehole Node: ',NBN,
     &        ' <=> Field Node: (',IX,',',JX,',',KX,':',NX,')'
          ENDIF
!
!---      Borehole node centroid  ---
!
          XP_BHX = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
          YP_BHX = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
          ZP_BHX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
!
!---      Write borehole node length  ---
!
          DISTX = SQRT( (XP_BH(2,NBN)-XP_BH(1,NBN))**2 + 
     &      (YP_BH(2,NBN)-YP_BH(1,NBN))**2 + 
     &      (ZP_BH(2,NBN)-ZP_BH(1,NBN))**2 )
          WRITE(IWR,'(6X,A,1PE11.4,A)') 'Borehole Node Length = ',DISTX,
     &      ', m'
!
!---      Divide field node into 512 sub-blocks and compute
!         average distance to borehole node centroid  ---
!
          DBN_BH(NBN) = 5.D-1*(VOL(NX)**(1.D+0/3.D+0))
          ADMBX = 0.D+0
          DO KX = 1,8
          DO JX = 1,8
          DO IX = 1,8
            XDX = (REAL(IX)-5.D-1)/8.D+0
            YDX = (REAL(JX)-5.D-1)/8.D+0
            ZDX = (REAL(KX)-5.D-1)/8.D+0
            CALL TRILINEAR( XDX,YDX,ZDX,XE(1,NX),PX(1) ) 
            CALL TRILINEAR( XDX,YDX,ZDX,YE(1,NX),PX(2) ) 
            CALL TRILINEAR( XDX,YDX,ZDX,ZE(1,NX),PX(3) )
            DISTX = SQRT( (PX(1)-XP_BHX)**2 + (PX(2)-YP_BHX)**2 +
     &        (PX(3)-ZP_BHX)**2 )
            ADMBX = ADMBX + DISTX
          ENDDO
          ENDDO
          ENDDO
          DBN_BH(NBN) = ADMBX/5.12D+2
          IF( ISLC(35).EQ.1 ) THEN
            WRITE(IWR,'(6X,A,I9,3A,1PE11.4,A)') 
     &      'Node #: ',NX,'  Borehole to Node Distance, ',
     &      UNLN(1:NCH),' : ',DBN_BH(NBN),', m)'
          ENDIF
  710   CONTINUE
  720 CONTINUE
!
!---  Loop over boreholes  ---
!
      NC = 0
      DO NBH = 1,N_BH
!
!---    Loop over borehole nodes  ---
!
        DO NBN1 = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---      Borehole node volume ---
!
          INVX = INV_BH(NBN1)
          DIST1X = SQRT( (XP_BH(2,NBN1)-XP_BH(1,NBN1))**2 + 
     &      (YP_BH(2,NBN1)-YP_BH(1,NBN1))**2 + 
     &      (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))**2 )
!    
!---      Inner coaxial borehole  ---
!    
          IF( IS_BH(INVX).EQ.10000 .AND. IBN_BH(NBN1).EQ.0 ) THEN
            RBX = MAX( PAR_BH(2,INVX),1.D-20 )
            VOL_BH(NBN1) = DIST1X*GPI*(RBX**2)
!    
!---      Outer coaxial borehole  ---
!    
          ELSEIF( IS_BH(INVX).EQ.10000 .AND. IBN_BH(NBN1).NE.0 ) THEN
            VOL_BH(NBN1) = DIST1X*GPI*
     &        ((PAR_BH(4,INVX)**2)-(PAR_BH(2,INVX)**2))
          ELSE
            RBX = MAX( PAR_BH(2,INVX),1.D-20 )    
            VOL_BH(NBN1) = DIST1X*GPI*(RBX**2)
          ENDIF
!
!---      First borehole node, use the two endpoints
!         to find the starting plane  ---
!
          ID_BH1 = ID_BH(3,NBH)
          ID_BH2 = ID_BH(3,NBH)
!
!---      Coaxial boreholes  ---
!
          IF( IS_BH(INVX).EQ.10000 ) 
     &      ID_BH2 = (ID_BH(4,NBH)-ID_BH(3,NBH))/2 + ID_BH(3,NBH) + 1
          IF( NBN1.EQ.ID_BH1 .OR. NBN1.EQ.ID_BH2 ) THEN
            DISTX = SQRT( (XP_BH(2,NBN1)-XP_BH(1,NBN1))**2 + 
     &        (YP_BH(2,NBN1)-YP_BH(1,NBN1))**2 + 
     &        (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))**2 )
            VDX(1) = (XP_BH(2,NBN1)-XP_BH(1,NBN1))/DISTX
            VDX(2) = (YP_BH(2,NBN1)-YP_BH(1,NBN1))/DISTX
            VDX(3) = (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))/DISTX
            CALL VERTICES_BH( RBX,VDX,XEX,YEX,ZEX )
            DO M = 1,4
              XE_BH(M,NBN1) = XEX(M) + XP_BH(1,NBN1)
              YE_BH(M,NBN1) = YEX(M) + YP_BH(1,NBN1)
              ZE_BH(M,NBN1) = ZEX(M) + ZP_BH(1,NBN1)
            ENDDO
          ELSE
            NBN2 = NBN1-1
            XP1X = 5.D-1*(XP_BH(2,NBN1)+XP_BH(1,NBN1))
            YP1X = 5.D-1*(YP_BH(2,NBN1)+YP_BH(1,NBN1))
            ZP1X = 5.D-1*(ZP_BH(2,NBN1)+ZP_BH(1,NBN1))
            XP2X = 5.D-1*(XP_BH(2,NBN2)+XP_BH(1,NBN2))
            YP2X = 5.D-1*(YP_BH(2,NBN2)+YP_BH(1,NBN2))
            ZP2X = 5.D-1*(ZP_BH(2,NBN2)+ZP_BH(1,NBN2))
            DISTX = SQRT( (XP1X-XP2X)**2 + (YP1X-YP2X)**2 + 
     &        (ZP1X-ZP2X)**2 )
            VDX(1) = (XP1X-XP2X)/DISTX
            VDX(2) = (YP1X-YP2X)/DISTX
            VDX(3) = (ZP1X-ZP2X)/DISTX
            NC = NC + 1
            DO M = 1,3
              SRFN_BH(M,NC) = VDX(M)
            ENDDO
            CALL VERTICES_BH( RBX,VDX,XEX,YEX,ZEX )
            DO M = 1,4
              XE_BH(M,NBN1) = XEX(M) + XP_BH(1,NBN1)
              YE_BH(M,NBN1) = YEX(M) + YP_BH(1,NBN1)
              ZE_BH(M,NBN1) = ZEX(M) + ZP_BH(1,NBN1)
            ENDDO
          ENDIF
!
!---      Last borehole node, use the two endpoints
!         to find the ending plane  ---
!
          ID_BH1 = ID_BH(4,NBH)
          ID_BH2 = ID_BH(4,NBH)
!
!---      Coaxial boreholes  ---
!
          IF( IS_BH(INVX).EQ.10000 ) 
     &      ID_BH2 = (ID_BH(4,NBH)-ID_BH(3,NBH))/2 + ID_BH(3,NBH)
          IF( NBN1.EQ.ID_BH1 .OR. NBN1.EQ.ID_BH2 ) THEN
            DISTX = SQRT( (XP_BH(2,NBN1)-XP_BH(1,NBN1))**2 + 
     &        (YP_BH(2,NBN1)-YP_BH(1,NBN1))**2 + 
     &        (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))**2 )
            VDX(1) = (XP_BH(2,NBN1)-XP_BH(1,NBN1))/DISTX
            VDX(2) = (YP_BH(2,NBN1)-YP_BH(1,NBN1))/DISTX
            VDX(3) = (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))/DISTX
            CALL VERTICES_BH( RBX,VDX,XEX,YEX,ZEX )
            DO M = 1,4
              XE_BH(M+4,NBN1) = XEX(M) + XP_BH(2,NBN1)
              YE_BH(M+4,NBN1) = YEX(M) + YP_BH(2,NBN1)
              ZE_BH(M+4,NBN1) = ZEX(M) + ZP_BH(2,NBN1)
            ENDDO
          ELSE
            NBN2 = NBN1+1
            XP1X = 5.D-1*(XP_BH(2,NBN1)+XP_BH(1,NBN1))
            YP1X = 5.D-1*(YP_BH(2,NBN1)+YP_BH(1,NBN1))
            ZP1X = 5.D-1*(ZP_BH(2,NBN1)+ZP_BH(1,NBN1))
            XP2X = 5.D-1*(XP_BH(2,NBN2)+XP_BH(1,NBN2))
            YP2X = 5.D-1*(YP_BH(2,NBN2)+YP_BH(1,NBN2))
            ZP2X = 5.D-1*(ZP_BH(2,NBN2)+ZP_BH(1,NBN2))
            DISTX = SQRT( (XP1X-XP2X)**2 + (YP1X-YP2X)**2 + 
     &        (ZP1X-ZP2X)**2 )
            VDX(1) = (XP2X-XP1X)/DISTX
            VDX(2) = (YP2X-YP1X)/DISTX
            VDX(3) = (ZP2X-ZP1X)/DISTX
            NC = NC + 1
            DO M = 1,3
              SRFN_BH(M,NC) = VDX(M)
            ENDDO
            CALL VERTICES_BH( RBX,VDX,XEX,YEX,ZEX )
            DO M = 1,4
              XE_BH(M+4,NBN1) = XEX(M) + XP_BH(2,NBN1)
              YE_BH(M+4,NBN1) = YEX(M) + YP_BH(2,NBN1)
              ZE_BH(M+4,NBN1) = ZEX(M) + ZP_BH(2,NBN1)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
!
!---  Loop over boreholes  ---
!
      NC = 0
      DO NBH = 1,N_BH
!
!---    Include inner to outer node connections for coaxial
!       boreholes  ---
!
        IF( IT_BH(1,NBH).GE.10000 ) THEN
          I1X = ID_BH(3,NBH)
          I2X = (ID_BH(4,NBH)-ID_BH(3,NBH))/2 + ID_BH(3,NBH)
          DO NBN = I1X,I2X
            NBN1 = NBN
            NBN2 = ID_BH(4,NBH) - NBN1 + 1
            DIST1X = SQRT( (XP_BH(2,NBN1)-XP_BH(1,NBN1))**2 + 
     &        (YP_BH(2,NBN1)-YP_BH(1,NBN1))**2 + 
     &        (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))**2 )
            DIST2X = SQRT( (XP_BH(2,NBN2)-XP_BH(1,NBN2))**2 + 
     &        (YP_BH(2,NBN2)-YP_BH(1,NBN2))**2 + 
     &        (ZP_BH(2,NBN2)-ZP_BH(1,NBN2))**2 )
            IF( NBN.EQ.I1X .AND. NBN.EQ.I2X ) THEN
              NC = NC + 1
              IPB_BH(1,NBN1) = 0
              IPB_BH(2,NBN1) = NC
              IBCM_BH(NC) = NBN2
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST2X)
              DBBM_BH(NC) = 5.D-1*DIST1X
              IPB_BH(3,NBN1) = -2             
              NC = NC + 1
              IPB_BH(1,NBN2) = NC
              IPB_BH(2,NBN2) = 0
              IBCM_BH(NC) = NBN1
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST2X)
              DBBM_BH(NC) = 5.D-1*DIST2X
              IPB_BH(3,NBN2) = -1             
            ELSEIF( NBN.EQ.I1X ) THEN
              DIST12X = SQRT( (XP_BH(2,NBN1+1)-XP_BH(1,NBN1+1))**2 + 
     &          (YP_BH(2,NBN1+1)-YP_BH(1,NBN1+1))**2 + 
     &          (ZP_BH(2,NBN1+1)-ZP_BH(1,NBN1+1))**2 )
              DIST22X = SQRT( (XP_BH(2,NBN2-1)-XP_BH(1,NBN2-1))**2 + 
     &          (YP_BH(2,NBN2-1)-YP_BH(1,NBN2-1))**2 + 
     &          (ZP_BH(2,NBN2-1)-ZP_BH(1,NBN2-1))**2 )
              NC = NC + 1
              IPB_BH(1,NBN1) = 0
              IPB_BH(2,NBN1) = NC
              IBCM_BH(NC) = NBN1 + 1
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST12X)
              DBBM_BH(NC) = 5.D-1*DIST1X
              NC = NC + 1
              IPB_BH(3,NBN1) = NC
              IBCM_BH(NC) = NBN2
              NC = NC + 1
              IPB_BH(1,NBN2) = NC
              IPB_BH(2,NBN2) = 0
              IBCM_BH(NC) = NBN2 - 1
              DBB_BH(NC) = 5.D-1*(DIST2X+DIST22X)
              DBBM_BH(NC) = 5.D-1*DIST2X
              NC = NC + 1
              IPB_BH(3,NBN2) = NC
              IBCM_BH(NC) = NBN1
            ELSEIF( NBN.EQ.I2X ) THEN
              DIST12X = SQRT( (XP_BH(2,NBN1-1)-XP_BH(1,NBN1-1))**2 + 
     &          (YP_BH(2,NBN1-1)-YP_BH(1,NBN1-1))**2 + 
     &          (ZP_BH(2,NBN1-1)-ZP_BH(1,NBN1-1))**2 )
              DIST22X = SQRT( (XP_BH(2,NBN2-2)-XP_BH(1,NBN2-2))**2 + 
     &          (YP_BH(2,NBN2+1)-YP_BH(1,NBN2+1))**2 + 
     &          (ZP_BH(2,NBN2+1)-ZP_BH(1,NBN2+1))**2 )
              NC = NC + 1
              IPB_BH(1,NBN1) = NC
              IBCM_BH(NC) = NBN1 - 1
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST12X)
              DBBM_BH(NC) = 5.D-1*DIST1X
              NC = NC + 1
              IPB_BH(2,NBN1) = NC
              IBCM_BH(NC) = NBN2
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST2X)
              DBBM_BH(NC) = 5.D-1*DIST1X
              IPB_BH(3,NBN1) = -2
              NC = NC + 1
              IPB_BH(1,NBN2) = NC
              IBCM_BH(NC) = NBN1
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST2X)
              DBBM_BH(NC) = 5.D-1*DIST2X
              NC = NC + 1
              IPB_BH(2,NBN2) = NC
              IBCM_BH(NC) = NBN2 + 1
              DBB_BH(NC) = 5.D-1*(DIST2X+DIST22X)
              DBBM_BH(NC) = 5.D-1*DIST2X
              IPB_BH(3,NBN2) = -1
            ELSE
              DIST12X = SQRT( (XP_BH(2,NBN1-1)-XP_BH(1,NBN1-1))**2 + 
     &          (YP_BH(2,NBN1-1)-YP_BH(1,NBN1-1))**2 + 
     &          (ZP_BH(2,NBN1-1)-ZP_BH(1,NBN1-1))**2 )
              NC = NC + 1
              IPB_BH(1,NBN1) = NC
              IBCM_BH(NC) = NBN1 - 1
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST12X)
              DBBM_BH(NC) = 5.D-1*DIST1X
              DIST12X = SQRT( (XP_BH(2,NBN1+1)-XP_BH(1,NBN1+1))**2 + 
     &          (YP_BH(2,NBN1+1)-YP_BH(1,NBN1+1))**2 + 
     &          (ZP_BH(2,NBN1+1)-ZP_BH(1,NBN1+1))**2 )
              NC = NC + 1
              IPB_BH(2,NBN1) = NC
              IBCM_BH(NC) = NBN1 + 1
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST12X)
              DBBM_BH(NC) = 5.D-1*DIST1X
              NC = NC + 1
              IPB_BH(3,NBN1) = NC
              IBCM_BH(NC) = NBN2
              DIST22X = SQRT( (XP_BH(2,NBN2-1)-XP_BH(1,NBN2-1))**2 + 
     &          (YP_BH(2,NBN2-1)-YP_BH(1,NBN2-1))**2 + 
     &          (ZP_BH(2,NBN2-1)-ZP_BH(1,NBN2-1))**2 )
              NC = NC + 1
              IPB_BH(1,NBN2) = NC
              IBCM_BH(NC) = NBN2 - 1
              DBB_BH(NC) = 5.D-1*(DIST2X+DIST22X)
              DBBM_BH(NC) = 5.D-1*DIST2X
              DIST22X = SQRT( (XP_BH(2,NBN2+1)-XP_BH(1,NBN2+1))**2 + 
     &          (YP_BH(2,NBN2+1)-YP_BH(1,NBN2+1))**2 + 
     &          (ZP_BH(2,NBN2+1)-ZP_BH(1,NBN2+1))**2 )
              NC = NC + 1
              IPB_BH(2,NBN2) = NC
              IBCM_BH(NC) = NBN2 + 1
              DBB_BH(NC) = 5.D-1*(DIST2X+DIST22X)
              DBBM_BH(NC) = 5.D-1*DIST2X
              NC = NC + 1
              IPB_BH(3,NBN2) = NC
              IBCM_BH(NC) = NBN1
            ENDIF
          ENDDO
        ELSE
!
!---      Loop over borehole nodes  ---
!
          DO NBN1 = ID_BH(3,NBH),ID_BH(4,NBH)
            DIST1X = SQRT( (XP_BH(2,NBN1)-XP_BH(1,NBN1))**2 + 
     &        (YP_BH(2,NBN1)-YP_BH(1,NBN1))**2 + 
     &        (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))**2 )
!
!---        Borehole node to borehole node connection map  ---
!
            IF( NBN1.EQ.ID_BH(3,NBH) .AND. NBN1.EQ.ID_BH(4,NBH) ) THEN
              IPB_BH(1,NBN1) = 0
              IPB_BH(2,NBN1) = 0
              IPB_BH(3,NBN1) = 0
            ELSEIF( NBN1.EQ.ID_BH(3,NBH) ) THEN
              NBN2 = NBN1 + 1
              NC = NC + 1
              IPB_BH(1,NBN1) = 0
              IPB_BH(2,NBN1) = NC
              IPB_BH(3,NBN1) = 0
              IBCM_BH(NC) = NBN2
              DIST2X = SQRT( (XP_BH(2,NBN2)-XP_BH(1,NBN2))**2 + 
     &          (YP_BH(2,NBN2)-YP_BH(1,NBN2))**2 + 
     &          (ZP_BH(2,NBN2)-ZP_BH(1,NBN2))**2 )
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST2X)
              DBBM_BH(NC) = 5.D-1*DIST1X
            ELSEIF( NBN1.EQ.ID_BH(4,NBH) ) THEN
              NBN2 = NBN1 - 1
              NC = NC + 1
              IPB_BH(1,NBN1) = NC
              IPB_BH(2,NBN1) = 0
              IPB_BH(3,NBN1) = 0
              IBCM_BH(NC) = NBN2
              DIST2X = SQRT( (XP_BH(2,NBN2)-XP_BH(1,NBN2))**2 + 
     &          (YP_BH(2,NBN2)-YP_BH(1,NBN2))**2 + 
     &          (ZP_BH(2,NBN2)-ZP_BH(1,NBN2))**2 )
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST2X)
              DBBM_BH(NC) = 5.D-1*DIST1X
            ELSE
              NBN2 = NBN1 - 1
              NC = NC + 1
              IPB_BH(1,NBN1) = NC
              IBCM_BH(NC) = NBN2
              DIST2X = SQRT( (XP_BH(2,NBN2)-XP_BH(1,NBN2))**2 + 
     &          (YP_BH(2,NBN2)-YP_BH(1,NBN2))**2 + 
     &          (ZP_BH(2,NBN2)-ZP_BH(1,NBN2))**2 )
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST2X)
              DBBM_BH(NC) = 5.D-1*DIST1X
              NBN2 = NBN1 + 1
              NC = NC + 1
              IPB_BH(2,NBN1) = NC
              IBCM_BH(NC) = NBN2
              DIST2X = SQRT( (XP_BH(2,NBN2)-XP_BH(1,NBN2))**2 + 
     &          (YP_BH(2,NBN2)-YP_BH(1,NBN2))**2 + 
     &          (ZP_BH(2,NBN2)-ZP_BH(1,NBN2))**2 )
              DBB_BH(NC) = 5.D-1*(DIST1X+DIST2X)
              DBBM_BH(NC) = 5.D-1*DIST1X
              IPB_BH(3,NBN1) = 0
            ENDIF
          ENDDO
        ENDIF
      ENDDO
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Loop over borehole nodes  ---
!
        DISTX = 0.D+0
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
          DIST1X = 5.D-1*SQRT( (XP_BH(2,NBN)-XP_BH(1,NBN))**2 + 
     &      (YP_BH(2,NBN)-YP_BH(1,NBN))**2 + 
     &      (ZP_BH(2,NBN)-ZP_BH(1,NBN))**2 )
          DISTX = DISTX + DIST1X
          DIST_BH(NBN) = DISTX
          DISTX = DISTX + DIST1X
        ENDDO
      ENDDO
!
!---  Loop over boreholes, checking for intersections of the
!     borehole nodes with fracture triangle surfaces  ---
!
      NC = 0
      DO NBH = 1,N_BH
        ID_BH(5,NBH) = NC+1
!
!---    Loop over number of borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---      Loop over fractures  ---
!
          DO NFX = 1,NF_FRC
!
!---        Loop over fracture triangles  ---
!
            DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---          Skip inactive triangles  ---
!
              IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---          Load plane-line intersection matrix and problem vector  ---
!
              AJ(1,1) = XP_BH(1,NBN)-XP_BH(2,NBN)
              AJ(2,1) = YP_BH(1,NBN)-YP_BH(2,NBN)
              AJ(3,1) = ZP_BH(1,NBN)-ZP_BH(2,NBN)
              AJ(1,2) = XE_FRC(2,NTX)-XE_FRC(1,NTX)
              AJ(2,2) = YE_FRC(2,NTX)-YE_FRC(1,NTX)
              AJ(3,2) = ZE_FRC(2,NTX)-ZE_FRC(1,NTX)
              AJ(1,3) = XE_FRC(3,NTX)-XE_FRC(1,NTX)
              AJ(2,3) = YE_FRC(3,NTX)-YE_FRC(1,NTX)
              AJ(3,3) = ZE_FRC(3,NTX)-ZE_FRC(1,NTX)
              BJ(1) = XP_BH(1,NBN)-XE_FRC(1,NTX)
              BJ(2) = YP_BH(1,NBN)-YE_FRC(1,NTX)
              BJ(3) = ZP_BH(1,NBN)-ZE_FRC(1,NTX)
!
!---          Check for no intersection  ---
!
              DO IP = 1,3
                AAMAX = 0.D+0
                DO JP = 1,3
                  IF( ABS(AJ(IP,JP)).GT.AAMAX ) AAMAX = ABS(AJ(IP,JP))
                ENDDO
!
!---            No intersection go to next triangle on the 
!               surface  ---
!
                IF( ABS(AAMAX)/EPSL.LT.EPSL ) GOTO 800
              ENDDO
!
!---          Find plane-line intersection matrix inverse  ---
!
              JP = 3
              KP = 3
              CALL LUDCMP( AJ,JP,KP,IJ,DJ )
              CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---          Find plane-line intersection point  ---
!
              TX = BJ(1)
              UX = BJ(2)
              VX = BJ(3)
              IF( ABS(TX).LT.EPSL ) TX = 0.D+0
              IF( ABS(UX).LT.EPSL ) UX = 0.D+0
              IF( ABS(VX).LT.EPSL ) VX = 0.D+0
!
!---          Line crosses surface, within the triangle  ---
!
              IF( TX.GE.0.D+0 .AND. TX.LE.1.D+0 .AND.
     &          UX.GE.0.D+0 .AND. UX.LE.1.D+0 .AND. 
     &          VX.GE.0.D+0 .AND. VX.LE.1.D+0 .AND.
     &          (UX+VX).LE.1.D+0 ) THEN
                NC = NC + 1
                XP_BF(NC) = XP_BH(1,NBN) 
     &            + (XP_BH(2,NBN)-XP_BH(1,NBN))*TX
                YP_BF(NC) = YP_BH(1,NBN)
     &            + (YP_BH(2,NBN)-YP_BH(1,NBN))*TX
                ZP_BF(NC) = ZP_BH(1,NBN)
     &            + (ZP_BH(2,NBN)-ZP_BH(1,NBN))*TX
                IBHT_FRC(NC) = NTX
                IBHN_FRC(NC) = NBN
              ENDIF
  800         CONTINUE
            ENDDO
          ENDDO
        ENDDO
        ID_BH(6,NBH) = NC
      ENDDO
!
!---  Record borehole-fracture triangle connections to output file ---
!
      WRITE(IWR,'(//,A,/)') ' --- Borehole Fracture Triangle ' // 
     &  'Record  ---'
      DO NBH = 1,N_BH
        NCH = INDEX( NM_BH(NBH),'  ') - 1
        WRITE(IWR,'(2X,A,I4,2A)') 'Fracture Triangles Traversed by ' //
     &    'Borehole Number: ',NBH,' Borehole Name: ',NM_BH(NBH)(1:NCH)
        DO NBTC = ID_BH(5,NBH),ID_BH(6,NBH)
          NTX = IBHT_FRC(NBTC)
          NBN = IBHN_FRC(NBTC)
          DO NFX = NF_FRC,1,-1
            NLTX = NTX - IP_FRC(1,NFX) + 1
            IF( NLTX.GT.0 ) EXIT
          ENDDO
          DO NBX = N_BH,1,-1
            NLBX = NBN - ID_BH(3,NBX) + 1
            IF( NLBX.GT.0 ) EXIT
          ENDDO
          WRITE(FORM2(5:5),'(I1)') ICOUNT(NFX)
          WRITE(FORM2(10:10),'(I1)') ICOUNT(NLTX)
          WRITE(FORM2(15:15),'(I1)') ICOUNT(NTX)
          WRITE(FORM2(20:20),'(I1)') ICOUNT(NBX)
          WRITE(FORM2(25:25),'(I1)') ICOUNT(NLBX)
          WRITE(FORM2(30:30),'(I1)') ICOUNT(NBN)
          WRITE(IWR,FORM2) ' (Fracture: ',NFX,
     &      ' Local Triangle: ',NLTX,
     &      ' Global Triangle: ',NTX,
     &      ' Borehole: ',NBX,
     &      ' Local Borehole Node: ',NLBX,
     &      ' Global Borehole Node: ',NBN,')'
        ENDDO
      ENDDO
!
!---  Open 'connect_bh' file for writing  ---
!
      OPEN(UNIT=27, FILE='connect_bh', STATUS='UNKNOWN', 
     &  FORM='FORMATTED')
      CLOSE(UNIT=27,STATUS='DELETE')
      OPEN(UNIT=27, FILE='connect_bh', STATUS='NEW', 
     &  FORM='FORMATTED')
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---      Lower surface  ---
!
          IF( NBN.EQ.ID_BH(3,NBH) ) THEN
            NODX(1) = (NBN-1)*8 + 1
            NODX(2) = (NBN-1)*8 + 2
            NODX(3) = (NBN-1)*8 + 4
            NODX(4) = (NBN-1)*8 + 3
          ELSE          
            NODX(1) = (NBN-2)*8 + 5
            NODX(2) = (NBN-2)*8 + 6
            NODX(3) = (NBN-2)*8 + 8
            NODX(4) = (NBN-2)*8 + 7
          ENDIF
!
!---      Upper surface  ---
!
          NODX(5) = (NBN-1)*8 + 5
          NODX(6) = (NBN-1)*8 + 6
          NODX(7) = (NBN-1)*8 + 8
          NODX(8) = (NBN-1)*8 + 7
          WRITE(27,'(8(I9,1X))') (NODX(I),I=1,8)
        ENDDO
      ENDDO
!
!---  Close 'connect_bh' file for writing  ---
!  
      CLOSE(UNIT=27)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_BOREHOLE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CONST_INTER( P1X,Q1X,R1X,P2X,Q2X,R2X,VN1X,VN2X,
     &  PE1X,PE2X,ITX )
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
!     Construct intersection
!
!     This subroutine is called when the triangles surely intersect
!     It constructs the segment of intersection of the two triangles
!     if they are not coplanar.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 March 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 P1X(3),Q1X(3),R1X(3)
      REAL*8 P2X(3),Q2X(3),R2X(3)
      REAL*8 V1X(3),V2X(3),VX(3)
      REAL*8 VNX(3),VN1X(3),VN2X(3)
      REAL*8 PE1X(3),PE2X(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CONST_INTER'
      CALL VSUB( Q1X,P1X,V1X )
      CALL VSUB( R2X,P1X,V2X )
      CALL VCROSSP( V1X,V2X,VNX )
      CALL VSUB( P2X,P1X,VX )
      IF( VDOTP(VX,VNX).GT.0.D+0 ) THEN
        CALL VSUB( R1X,P1X,V1X )
        CALL VCROSSP( V1X,V2X,VNX )
        IF( VDOTP(VX,VNX).LE.0.D+0 ) THEN
          CALL VSUB( Q2X,P1X,V2X )
          CALL VCROSSP( V1X,V2X,VNX )
          IF( VDOTP(VX,VNX).GT.0.D+0 ) THEN
            CALL VSUB( P1X,P2X,V1X )
            CALL VSUB( P1X,R1X,V2X )
            ALPHAX = VDOTP( V1X,VN2X )/VDOTP( V2X,VN2X )
            CALL VSCALAR( V2X,ALPHAX,V1X )
            CALL VSUB( P1X,V1X,PE1X )
            CALL VSUB( P2X,P1X,V1X )
            CALL VSUB( P2X,R2X,V2X )
            ALPHAX = VDOTP( V1X,VN1X )/VDOTP( V2X,VN1X )
            CALL VSCALAR( V2X,ALPHAX,V1X )
            CALL VSUB( P2X,V1X,PE2X )
            ITX = 1
          ELSE
            CALL VSUB( P2X,P1X,V1X )
            CALL VSUB( P2X,Q2X,V2X )
            ALPHAX = VDOTP( V1X,VN1X )/VDOTP( V2X,VN1X )
            CALL VSCALAR( V2X,ALPHAX,V1X )
            CALL VSUB( P2X,V1X,PE1X )
            CALL VSUB( P2X,P1X,V1X )
            CALL VSUB( P2X,R2X,V2X )
            ALPHAX = VDOTP( V1X,VN1X )/VDOTP( V2X,VN1X )
            CALL VSCALAR( V2X,ALPHAX,V1X )
            CALL VSUB( P2X,V1X,PE2X )
            ITX = 1
          ENDIF
        ELSE
          ITX = 0
        ENDIF
      ELSE
        CALL VSUB( Q2X,P1X,V2X )
        CALL VCROSSP( V1X,V2X,VNX )
        IF( VDOTP( VX,VNX ).LT.0.D+0 ) THEN
          ITX = 0
        ELSE
          CALL VSUB( R1X,P1X,V1X )
          CALL VCROSSP( V1X,V2X,VNX )
          IF( VDOTP( VX,VNX ).GE.0.D+0 ) THEN
            CALL VSUB( P1X,P2X,V1X )
            CALL VSUB( P1X,R1X,V2X )
            ALPHAX = VDOTP( V1X,VN2X )/VDOTP( V2X,VN2X )
            CALL VSCALAR( V2X,ALPHAX,V1X )
            CALL VSUB( P1X,V1X,PE1X )
            CALL VSUB( P1X,P2X,V1X )
            CALL VSUB( P1X,Q1X,V2X )
            ALPHAX = VDOTP( V1X,VN2X )/VDOTP( V2X,VN2X )
            CALL VSCALAR( V2X,ALPHAX,V1X )
            CALL VSUB( P1X,V1X,PE2X )
            ITX = 1
          ELSE
            CALL VSUB( P2X,P1X,V1X )
            CALL VSUB( P2X,Q2X,V2X )
            ALPHAX = VDOTP( V1X,VN1X )/VDOTP( V2X,VN1X )
            CALL VSCALAR( V2X,ALPHAX,V1X )
            CALL VSUB( P2X,V1X,PE1X )
            CALL VSUB( P1X,P2X,V1X )
            CALL VSUB( P1X,Q1X,V2X )
            ALPHAX = VDOTP( V1X,VN2X )/VDOTP( V2X,VN2X )
            CALL VSCALAR( V2X,ALPHAX,V1X )
            CALL VSUB( P1X,V1X,PE2X )
            ITX = 1
          ENDIF
        ENDIF
      ENDIF            
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CONST_INTER group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE COPL_TRI( PN1X,V1X,V2X,PTX,NPTX,ITX )
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
!     PN1X - normal vector to triangle V
!     V1X - x,y,z coordinates of vertices of triangle V
!     V2X - x,y,z coordinates of vertices of triangle U
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 PN1X(3)
      REAL*8 V1X(3,3),V2X(3,3)
      REAL*8 AX(3)
      REAL*8, DIMENSION(3,24) :: PTX
      INTEGER MX(3,2)
!
!----------------------Data Statements---------------------------------!
!
      DATA MX /1,2,3,2,3,1/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/COPL_TRI'
!
!---  First project onto an axis aligned plane that maximizes the
!     areas of the triangles, computing indices I1X and I2X  ---
!
      DO M = 1,3
        AX(M) = ABS(PN1X(M))
      ENDDO
      IF( AX(1).GT.AX(2) ) THEN
!
!---    AX(1) is the greatest  ---
!
        IF( AX(1).GT.AX(3) ) THEN
          I0X = 2
          I1X = 3
!
!---    AX(3) is the greatest  ---
!
        ELSE
          I0X = 1
          I1X = 2
        ENDIF
      ELSE
!
!---    AX(3) is the greatest  ---
!
        IF( AX(3).GT.AX(2) ) THEN
          I0X = 1
          I1X = 2
!
!---    AX(2) is the greatest  ---
!
        ELSE
          I0X = 1
          I1X = 3
        ENDIF
      ENDIF
!
!---  Check if vertices of triangle 1 are within triangle 2  ---
!
      NC1X = 0
      DO M = 1,3
        CALL PNT_TRI_CHK( V1X(1,M),V2X,I0X,I1X,ITX )
        IF( ITX.EQ.1 ) THEN
          NC1X = NC1X + 1
!
!---      Skip duplicate points  ---
!
          ISKPX = 0
          DO N = 1,NPTX
            IF( ABS(V1X(1,M)-PTX(1,N)).LT.EPSL .AND. 
     &        ABS(V1X(2,M)-PTX(2,N)).LT.EPSL .AND.
     &        ABS(V1X(3,M)-PTX(3,N)).LT.EPSL ) THEN
              ISKPX = 1
              EXIT
            ENDIF
          ENDDO
          IF( ISKPX.EQ.0 ) THEN
            NPTX = NPTX + 1
            PTX(1,NPTX) = V1X(1,M)
            PTX(2,NPTX) = V1X(2,M)
            PTX(3,NPTX) = V1X(3,M)
          ENDIF
        ENDIF
      ENDDO
!
!---  Check if vertices of triangle 2 are within triangle 1  ---
!
      NC2X = 0
      DO M = 1,3
        CALL PNT_TRI_CHK( V2X(1,M),V1X,I0X,I1X,ITX )
        IF( ITX.EQ.1 ) THEN
          NC2X = NC2X + 1
!
!---      Skip duplicate points  ---
!
          ISKPX = 0
          DO N = 1,NPTX
            IF( ABS(V2X(1,M)-PTX(1,N)).LT.EPSL .AND. 
     &        ABS(V2X(2,M)-PTX(2,N)).LT.EPSL .AND.
     &        ABS(V2X(3,M)-PTX(3,N)).LT.EPSL ) THEN
              ISKPX = 1
              EXIT
            ENDIF
          ENDDO
          IF( ISKPX.EQ.0 ) THEN
            NPTX = NPTX + 1
            PTX(1,NPTX) = V2X(1,M)
            PTX(2,NPTX) = V2X(2,M)
            PTX(3,NPTX) = V2X(3,M)
          ENDIF
        ENDIF
      ENDDO
!
!---  Check for singlet or doublet, indicating triangle-triangle
!     intersection  ---
!
      ITX = 0
      IF( NC1X.EQ.1.OR.NC1X.EQ.2 .OR. NC2X.EQ.1.OR.NC2X.EQ.2 ) THEN
        ITX = 1
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Check for intersections of line segments of triangle 1 with
!     line segments of triangle 2  ---
!
      DO M1X = 1,3
        DO M2X = 1,3
          CALL LSEG_CHK( V1X(1,MX(M1X,1)),V1X(1,MX(M1X,2)),
     &      V2X(1,MX(M2X,1)),V2X(1,MX(M2X,2)),I0X,I1X,ITX )
          IF( ITX.EQ.1 ) THEN
            ISUB_LOG = ISUB_LOG-1
            RETURN
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of COPL_TRI group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE COPL_INT( PN1X,V1X,V2X,PTX,NPTX )
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
!     Find the x,y,z coordinates of the four intersection points of 
!     two triangles along the intersection line.
!
!     PN1X - normal vector to triangle 1 or plane Pi1
!     V1X - x,y,z coordinates of vertices of triangle 1
!     V2X - x,y,z coordinates of vertices of triangle 2
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 August 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(3,24) :: PTX
      REAL*8 PN1X(3)
      REAL*8 XPNCX(3),YPNCX(3),ZPNCX(3)
      REAL*8 XNCX(3),YNCX(3),ZNCX(3)
      REAL*8 PIX(2,6)
      REAL*8 V1X(3,3),V2X(3,3)
      REAL*8 V1PX(2,3),V2PX(2,3)
      REAL*8 VCX(3),VDX(3)
      REAL*8 TRNSMX(3,3)
      INTEGER MX(3,2)
!
!----------------------Data Statements---------------------------------!
!
      DATA MX /1,2,3,2,3,1/
      DATA XNCX /1,0,0/
      DATA YNCX /0,1,0/
      DATA ZNCX /0,0,1/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/COPL_INT'
!
!---  Create a x'y' coordinate system in the plane with vertex 1 of
!     triangle 1 as the center, and the vector between vertex 1 and 2
!     of triangle 1 as the x' coordinate direction. The y' coordinate 
!     direction will be normal to the plane surface normal, and the 
!     x' coordinate direction  ---
!
      CALL VSUB( V1X(1,2),V1X(1,1),XPNCX )
      DXPCX = SQRT( (V1X(1,2)-V1X(1,1))**2 + (V1X(2,2)-V1X(2,1))**2
     &  + (V1X(3,2)-V1X(3,1))**2 )
      DZPCX = SQRT( PN1X(1)**2 + PN1X(2)**2 + PN1X(3)**2 )
      DO M = 1,3
        XPNCX(M) = XPNCX(M)/DXPCX
        ZPNCX(M) = PN1X(M)/DZPCX
      ENDDO
      CALL VCROSSP( ZPNCX,XPNCX,YPNCX )
!
!---  Convert triangles to the x'y' coordinate system  ---
!
      DO M = 1,3
        CALL VSUB( V1X(1,M),V1X(1,1),VDX )
        V1PX(1,M) = VDOTP(VDX,XPNCX)
        V1PX(2,M) = VDOTP(VDX,YPNCX)
        CALL VSUB( V2X(1,M),V1X(1,1),VDX )
        V2PX(1,M) = VDOTP(VDX,XPNCX)
        V2PX(2,M) = VDOTP(VDX,YPNCX)
      ENDDO
!
!---  Check for intersections of the triangle lines  ---
!
      NPIX = 0
      DO M1X = 1,3
        DO M2X = 1,3
          CALL LSEG_INT( V1PX(1,MX(M1X,1)),V1PX(1,MX(M1X,2)),
     &      V2PX(1,MX(M2X,1)),V2PX(1,MX(M2X,2)),PIX,NPIX )
        ENDDO
      ENDDO
!
!---  Convert intersection points to x,y,z coordinate system  ---
!
      TRNSMX(1,1) = VDOTP( XNCX,XPNCX )
      TRNSMX(1,2) = VDOTP( XNCX,YPNCX )
      TRNSMX(1,3) = VDOTP( XNCX,ZPNCX )
      TRNSMX(2,1) = VDOTP( YNCX,XPNCX )
      TRNSMX(2,2) = VDOTP( YNCX,YPNCX )
      TRNSMX(2,3) = VDOTP( YNCX,ZPNCX )
      TRNSMX(3,1) = VDOTP( ZNCX,XPNCX )
      TRNSMX(3,2) = VDOTP( ZNCX,YPNCX )
      TRNSMX(3,3) = VDOTP( ZNCX,ZPNCX )
      DO N = 1,NPIX
        VDX(1) = PIX(1,N)
        VDX(2) = PIX(2,N)
        VDX(3) = 0.D+0
!
!---    Rotation between x',y',z' and x,y,z  ---
!
        CALL MATMUL( TRNSMX,VDX,VCX,3,3,1 )
!
!---    Translation from origin of x',y',z'  ---
!
        DO M = 1,3
          VCX(M) = VCX(M) + V1X(M,1)
        ENDDO
!
!---    Skip duplicate points  ---
!
        ISKPX = 0
        DO M = 1,NPTX
          IF( ABS(VCX(1)-PTX(1,M)).LT.EPSL .AND. 
     &      ABS(VCX(2)-PTX(2,M)).LT.EPSL .AND.
     &      ABS(VCX(3)-PTX(3,M)).LT.EPSL ) THEN
            ISKPX = 1
            EXIT
          ENDIF
        ENDDO
        IF( ISKPX.EQ.0 ) THEN
          NPTX = NPTX + 1
          PTX(1,NPTX) = VCX(1)
          PTX(2,NPTX) = VCX(2)
          PTX(3,NPTX) = VCX(3)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of COPL_INT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CPTT3D( P1X,Q1X,R1X,P2X,Q2X,R2X,VN1X,VN2X,CPX,ICPX )
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
!     Find the intersection points of two coplanar triangles in 3D.
!     Triangle 1 (i.e., P1X,Q1X,R1X ) is the fracture/fault triangle
!     Triangle 2 (i.e., P2X,Q2X,R2X ) is the grid surface triangle
!
!     Triangle 1 TV1X(3,2) is 2D array of the fracture/fault triangle
!     Triangle 2 TV2X(3,2) is 2D array of the grid surface triangle
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 August 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VX(3)
      REAL*8 P1X(3),Q1X(3),R1X(3)
      REAL*8 P2X(3),Q2X(3),R2X(3)
      REAL*8 CP1X(3),CQ1X(3),CR1X(3)
      REAL*8 CP2X(3),CQ2X(3),CR2X(3)
      REAL*8 VN1X(3),VN2X(3),VN3X(3),VN4X(3)
      REAL*8 CPX(3,6),CPTX(2,6)
      REAL*8 ROTMATX(3,3),VMATX(3,3),V2MATX(3,3)
      REAL*8 TV1X(3,2),TV2X(3,2)
      INTEGER MP(2,3)
      LOGICAL SX_NEG,TX_NEG
!
!----------------------Data Statements---------------------------------!
!
      DATA MP / 1,2,2,3,3,1 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CPTT3D'
!
!---  Rotate both triangle vertices onto a principal axis, such that the
!     the area is maximized  ---
!
      VN_X = ABS(VN1X(1))
      VN_Y = ABS(VN1X(2))
      VN_Z = ABS(VN1X(3))
      SX = SQRT( VN1X(1)**2 + VN1X(2)**2 + VN1X(3)**2 )
      DO M = 1,3
        VN3X(M) = VN1X(M)/SX
      ENDDO
!
!     Rotate onto yz plane
!
      IF( VN_X.GT.VN_Z .AND. VN_X.GE.VN_Y ) THEN
        VN4X(1) = 1.D+0
        VN4X(2) = 0.D+0
        VN4X(3) = 0.D+0
!
!     Rotate onto xz plane
!
      ELSEIF( VN_Y.GT.VN_Z .AND. VN_Y.GE.VN_X ) THEN
        VN4X(1) = 0.D+0
        VN4X(2) = 1.D+0
        VN4X(3) = 0.D+0
!
!     Rotate onto xy plane
!
      ELSE
        VN4X(1) = 0.D+0
        VN4X(2) = 0.D+0
        VN4X(3) = 1.D+0
      ENDIF
!
!     Rotation matrix of triangle unit vector onto unit vector of
!     principal axis
!
      CALL VCROSSP( VN3X,VN4X,VX )
      SX = SQRT( VX(1)**2 + VX(2)**2 + VX(3)**2 )
      CX = VDOTP( VN3X,VN4X )
      IF( ABS(1.D+0+CX).LT.EPSL ) CX = 1.D+0
      DO I = 1,3
        DO J = 1,3
          ROTMATX(I,J) = 0.D+0
          VMATX(I,J) = 0.D+0
          IF( I.EQ.J ) ROTMATX(I,J) = 1.D+0
        ENDDO
      ENDDO
      VMATX(1,2) = -VX(3)
      VMATX(1,3) = VX(2)
      VMATX(2,1) = VX(3)
      VMATX(2,3) = -VX(1)
      VMATX(3,1) = -VX(2)
      VMATX(3,2) = VX(1)
      CALL MATMUL( VMATX,VMATX,V2MATX,3,3,3 )
      DO I = 1,3
        DO J = 1,3
          ROTMATX(I,J) = ROTMATX(I,J) + VMATX(I,J) + 
     &      V2MATX(I,J)*(1.D+0/(1.D+0+CX))
        ENDDO
      ENDDO
      CALL MATMUL( ROTMATX,P1X,CP1X,3,3,1 )
      CALL MATMUL( ROTMATX,Q1X,CQ1X,3,3,1 )
      CALL MATMUL( ROTMATX,R1X,CR1X,3,3,1 )
      CALL MATMUL( ROTMATX,P2X,CP2X,3,3,1 )
      CALL MATMUL( ROTMATX,Q2X,CQ2X,3,3,1 )
      CALL MATMUL( ROTMATX,R2X,CR2X,3,3,1 )
!
!     Rotate onto yz plane
!
      IF( VN_X.GT.VN_Z .AND. VN_X.GE.VN_Y ) THEN
        TV1X(1,1) = CP1X(2)
        TV1X(1,2) = CP1X(3)
        TV1X(2,1) = CQ1X(2)
        TV1X(2,2) = CQ1X(3)
        TV1X(3,1) = CR1X(2)
        TV1X(3,2) = CR1X(3)
        TV2X(1,1) = CP2X(2)
        TV2X(1,2) = CP2X(3)
        TV2X(2,1) = CQ2X(2)
        TV2X(2,2) = CQ2X(3)
        TV2X(3,1) = CR2X(2)
        TV2X(3,2) = CR2X(3)
!
!     Rotate onto xz plane
!
      ELSEIF( VN_Y.GT.VN_Z .AND. VN_Y.GE.VN_X ) THEN
        TV1X(1,1) = CP1X(1)
        TV1X(1,2) = CP1X(3)
        TV1X(2,1) = CQ1X(1)
        TV1X(2,2) = CQ1X(3)
        TV1X(3,1) = CR1X(1)
        TV1X(3,2) = CR1X(3)
        TV2X(1,1) = CP2X(1)
        TV2X(1,2) = CP2X(3)
        TV2X(2,1) = CQ2X(1)
        TV2X(2,2) = CQ2X(3)
        TV2X(3,1) = CR2X(1)
        TV2X(3,2) = CR2X(3)
!
!     Rotate onto xy plane
!
      ELSE
        TV1X(1,1) = CP1X(1)
        TV1X(1,2) = CP1X(2)
        TV1X(2,1) = CQ1X(1)
        TV1X(2,2) = CQ1X(2)
        TV1X(3,1) = CR1X(1)
        TV1X(3,2) = CR1X(2)
        TV2X(1,1) = CP2X(1)
        TV2X(1,2) = CP2X(2)
        TV2X(2,1) = CQ2X(1)
        TV2X(2,2) = CQ2X(2)
        TV2X(3,1) = CR2X(1)
        TV2X(3,2) = CR2X(2)
      ENDIF
!
!---  Check to see if all triangle 1 points are within triangle 2
!     boundaries  ---
!
      ICPX = 0
      NC = 0
      DO I = 1,3
        SX = TV2X(1,2)*TV2X(3,1) - TV2X(1,1)*TV2X(3,2) + 
     &    (TV2X(3,2)-TV2X(1,2))*TV1X(I,1) + 
     &    (TV2X(1,1)-TV2X(3,1))*TV1X(I,2)
        TX = TV2X(1,1)*TV2X(2,2) - TV2X(1,2)*TV2X(2,1) + 
     &    (TV2X(1,2)-TV2X(2,2))*TV1X(I,1) + 
     &    (TV2X(2,1)-TV2X(1,1))*TV1X(I,2)
        SX_NEG = (SX.LT.0.D+0)
        TX_NEG = (TX.LT.0.D+0)
        IF( SX_NEG.NEQV.TX_NEG ) CYCLE
        AX = -TV2X(2,2)*TV2X(3,1) + TV2X(1,2)*(TV2X(3,1)-TV2X(2,1)) +
     &    TV2X(1,1)*(TV2X(2,2)-TV2X(3,2)) + TV2X(2,1)*TV2X(3,2)
        IF( AX.LT.0.D+0 ) THEN
          IF( SX.LE.0.D+0 .AND. (SX+TX).GE.AX ) THEN
            NC = NC + 1
            ICPX = ICPX + 1
            CPTX(1,ICPX) = TV1X(I,1)
            CPTX(2,ICPX) = TV1X(I,2)
          ELSE
            CYCLE
          ENDIF
        ELSE
          IF( SX.GE.0.D+0 .AND. (SX+TX).LE.AX ) THEN
            NC = NC + 1
            ICPX = ICPX + 1
            CPTX(1,ICPX) = TV1X(I,1)
            CPTX(2,ICPX) = TV1X(I,2)
          ELSE
            CYCLE
          ENDIF
        ENDIF
      ENDDO
      IF( NC.EQ.3 ) THEN
        DO M = 1,3
          CPX(M,1) = P1X(M)
          CPX(M,2) = Q1X(M)
          CPX(M,3) = R1X(M)
        ENDDO
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Check to see if all triangle 2 points are within triangle 1
!     boundaries  ---
!
      NC = 0
      DO I = 1,3
        SX = TV1X(1,2)*TV1X(3,1) - TV1X(1,1)*TV1X(3,2) + 
     &    (TV1X(3,2)-TV1X(1,2))*TV2X(I,1) + 
     &    (TV1X(1,1)-TV1X(3,1))*TV2X(I,2)
        TX = TV1X(1,1)*TV1X(2,2) - TV1X(1,2)*TV1X(2,1) + 
     &    (TV1X(1,2)-TV1X(2,2))*TV2X(I,1) + 
     &    (TV1X(2,1)-TV1X(1,1))*TV2X(I,2)
        SX_NEG = (SX.LT.0.D+0)
        TX_NEG = (TX.LT.0.D+0)
        IF( SX_NEG.NEQV.TX_NEG ) CYCLE
        AX = -TV1X(2,2)*TV1X(3,1) + TV1X(1,2)*(TV1X(3,1)-TV1X(2,1)) +
     &    TV1X(1,1)*(TV1X(2,2)-TV1X(3,2)) + TV1X(2,1)*TV1X(3,2)
        IF( AX.LT.0.D+0 ) THEN
          IF( SX.LE.0.D+0 .AND. (SX+TX).GE.AX ) THEN
            NC = NC + 1
            ICPX = ICPX + 1
            CPTX(1,ICPX) = TV2X(I,1)
            CPTX(2,ICPX) = TV2X(I,2)
          ELSE
            CYCLE
          ENDIF
        ELSE
          IF( SX.GE.0.D+0 .AND. (SX+TX).LE.AX ) THEN
            NC = NC + 1
            ICPX = ICPX + 1
            CPTX(1,ICPX) = TV2X(I,1)
            CPTX(2,ICPX) = TV2X(I,2)
          ELSE
            CYCLE
          ENDIF
        ENDIF
      ENDDO
      IF( NC.EQ.3 ) THEN
        ICPX = 3
        DO M = 1,3
          CPX(M,1) = P2X(M)
          CPX(M,2) = Q2X(M)
          CPX(M,3) = R2X(M)
        ENDDO
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Loop over the three line segments of triangle 1
!     checking for intersections with triangle 2  ---
!
      DO I = 1,3
        M1 = MP(1,I)
        M2 = MP(2,I)
        X1 = TV1X(M1,1)
        Y1 = TV1X(M1,2)
        X2 = TV1X(M2,1)
        Y2 = TV1X(M2,2)
        DO J = 1,3
          M3 = MP(1,J)
          M4 = MP(2,J)
          X3 = TV2X(M3,1)
          Y3 = TV2X(M3,2)
          X4 = TV2X(M4,1)
          Y4 = TV2X(M4,2)
          TX = ((X1-X3)*(Y3-Y4) - (Y1-Y3)*(X3-X4))/
     &      ((X1-X2)*(Y3-Y4)-(Y1-Y2)*(X3-X4)+SMALL)
          UX = ((X2-X1)*(Y1-Y3) - (Y2-Y1)*(X1-X3))/
     &      ((X1-X2)*(Y3-Y4)-(Y1-Y2)*(X3-X4)+SMALL)
          IF( TX.GE.0.D+0 .AND. TX.LE.1.D+0 .AND. 
     &      UX.GE.0.D+0 .AND. UX.LE.1.D+0 ) THEN
            ICPX = ICPX + 1
            CPTX(1,ICPX) = X1 + TX*(X2-X1)
            CPTX(2,ICPX) = Y1 + TX*(Y2-Y1)
          ENDIF
        ENDDO
      ENDDO
!
!---  Inverse rotation matrix (i.e., transpose of original
!     rotation matrix)  ---
!
      DO J = 1,3
        DO I = 1,3
          VMATX(J,I) = ROTMATX(I,J)
        ENDDO
      ENDDO
!
!---  Rotate 2D point back to 3D points  ---
!
      DO I = 1,ICPX
!
!       Rotate onto yz plane
!
        IF( VN_X.GT.VN_Z .AND. VN_X.GE.VN_Y ) THEN
          CP1X(2) = CPTX(1,I)
          CP1X(3) = CPTX(2,I)
!
!       Rotate onto xz plane
!
        ELSEIF( VN_Y.GT.VN_Z .AND. VN_Y.GE.VN_X ) THEN
          CP1X(1) = CPTX(1,I)
          CP1X(3) = CPTX(2,I)
!
!       Rotate onto xy plane
!
        ELSE
          CP1X(1) = CPTX(1,I)
          CP1X(2) = CPTX(2,I)
        ENDIF
        CALL MATMUL( VMATX,CP1X,P1X,3,3,1 )
        DO M = 1,3
          CPX(M,I) = P1X(M)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CPTT3D group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CRN_LIM_FRC( NSL )
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
!     Compute a time step that globally satisfies the Courant
!     number limit for fracture flow.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 3 July 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE GEOM_FRC
      USE FDVP_FRC
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
      SUB_LOG(ISUB_LOG) = '/CRN_LIM_FRC'
!
!---  Compute the maximum Courant number ---
!
      CRNLX = 0.D+0
      CRNNX = 0.D+0
      CRNGX = 0.D+0
!
!---  Standard Courant number control  ---
!
      IF( ISLC(17).EQ.1 ) THEN
!
!---    Loop over fractures  ---
!
        DO NFX = 1,NF_FRC
!
!---      Loop over fracture triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---        Check that fracture triangle has nonzero solute
!           concentrations  ---
!
            ICHK = 0
            DO NS = 1,NSOLU
              IF( C_FRC(NTX,NS).GT.1.D-20 ) ICHK = 1
            ENDDO
            IF( ICHK.EQ.0 ) CYCLE
!
!---        Active aqueous  ---
!
            IF( IAQU.EQ.1 ) CRNLX = MAX( CRNTL_FRC(NTX),CRNLX )
!
!---        Active gas  ---
!
            IF( IGAS.EQ.1 ) CRNGX = MAX( CRNTG_FRC(NTX),CRNGX )
!
!---        Active nonaqueous liquid  ---
!
            IF( INAPL.EQ.1 ) CRNNX = MAX( CRNTN_FRC(NTX),CRNNX )
          ENDDO
        ENDDO
!
!---  Special vadose zone Courant number control  ---
!
      ELSEIF( ISLC(17).EQ.2 .AND. IAQU.EQ.1 ) THEN
!
!---    Loop over fractures  ---
!
        DO NFX = 1,NF_FRC
!
!---      Loop over fracture triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---        Check for fully saturated conditions  ---
!
            IF( 1.D+0-SL_FRC(2,NTX).LT.EPSL ) CYCLE
            CLX = C_FRC(NTX,NSL)*YL_FRC(NTX,NSL)/(SL_FRC(2,NTX)+EPSL)
            IF( CLX.GT.CCL_CRN(NSL) ) THEN
              CRNLX = MAX( CRNTL_FRC(NTX),CRNLX )
            ENDIF
          ENDDO
        ENDDO
!
!---  Aqueous-only Courant number control  ---
!
      ELSEIF( ISLC(17).EQ.3 .AND. IAQU.EQ.1 ) THEN
!
!---    Loop over fractures  ---
!
        DO NFX = 1,NF_FRC
!
!---      Loop over fracture triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---        Check that fracture triangle has nonzero solute
!           concentrations  ---
!
            ICHK = 0
            DO NS = 1,NSOLU
              IF( C_FRC(NTX,NS).GT.1.D-20 ) ICHK = 1
            ENDDO
            IF( ICHK.EQ.0 ) CYCLE
!
!---        Active aqueous  ---
!
            IF( IAQU.EQ.1 ) CRNLX = MAX( CRNTL_FRC(NTX),CRNLX )
          ENDDO
        ENDDO
      ENDIF
      CRNMX = MAX( CRNLX,CRNGX,CRNNX )
      IF( CRNMX.GT.CRNTMXT ) THEN
        N_CRN(NSL) = MAX( INT(CRNMX/CRNTMXT) + 1,N_CRN(NSL) )
        REALX = REAL(N_CRN(NSL))
        DT = DT_CRN/REALX
        DTI = 1.D+0/(DT+EPSL)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CRN_LIM_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CRN_LIM_BH( NSL )
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
!     Compute a time step that globally satisfies the Courant
!     number limit for borehole flow.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE GEOM_BH
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/CRN_LIM_BH'
!
!---  Compute the maximum Courant number ---
!
      CRNLX = 0.D+0
      CRNNX = 0.D+0
      CRNGX = 0.D+0
!
!---  Standard Courant number control  ---
!
      IF( ISLC(17).EQ.1 ) THEN
!
!---    Loop over boreholes  ---
!
        DO NBH = 1,N_BH
!
!---      Skip for boundary condition type boreholes  ---
!
          IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---      Loop over borehole nodes  ---
!
          DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---        Check that borehole node has nonzero solute
!           concentrations  ---
!
            ICHK = 0
            DO NS = 1,NSOLU
              IF( C_BH(NBN,NS).GT.1.D-20 ) ICHK = 1
            ENDDO
            IF( ICHK.EQ.0 ) CYCLE
!
!---        Active aqueous  ---
!
            IF( IAQU.EQ.1 ) CRNLX = MAX( CRNTL_BH(NBN),CRNLX )
!
!---        Active gas  ---
!
            IF( IGAS.EQ.1 ) CRNGX = MAX( CRNTG_BH(NBN),CRNGX )
!
!---        Active nonaqueous liquid  ---
!
            IF( INAPL.EQ.1 ) CRNNX = MAX( CRNTN_BH(NBN),CRNNX )
          ENDDO
        ENDDO
!
!---  Special vadose zone Courant number control  ---
!
      ELSEIF( ISLC(17).EQ.2 .AND. IAQU.EQ.1 ) THEN
!
!---    Loop over boreholes  ---
!
        DO NBH = 1,N_BH
!
!---      Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---      Loop over borehole nodes  ---
!
          DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---        Check for fully saturated conditions  ---
!
            IF( 1.D+0-SL_BH(2,NBN).LT.EPSL ) CYCLE
            CLX = C_BH(NBN,NSL)*YL_BH(NBN,NSL)/(SL_BH(2,NBN)+EPSL)
            IF( CLX.GT.CCL_CRN(NSL) ) THEN
              CRNLX = MAX( CRNTL_BH(NBN),CRNLX )
            ENDIF
          ENDDO
        ENDDO
!
!---  Aqueous-only Courant number control  ---
!
      ELSEIF( ISLC(17).EQ.3 .AND. IAQU.EQ.1 ) THEN
!
!---    Loop over boreholes  ---
!
        DO NBH = 1,N_BH
!
!---      Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---      Loop over borehole nodes  ---
!
          DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---        Check that fracture triangle has nonzero solute
!           concentrations  ---
!
            ICHK = 0
            DO NS = 1,NSOLU
              IF( C_BH(NBN,NS).GT.1.D-20 ) ICHK = 1
            ENDDO
            IF( ICHK.EQ.0 ) CYCLE
!
!---        Active aqueous  ---
!
            IF( IAQU.EQ.1 ) CRNLX = MAX( CRNTL_BH(NBN),CRNLX )
          ENDDO
        ENDDO
      ENDIF
      CRNMX = MAX( CRNLX,CRNGX,CRNNX )
      IF( CRNMX.GT.CRNTMXT ) THEN
        N_CRN(NSL) = MAX( N_CRN(NSL),INT(CRNMX/CRNTMXT) + 1 )
        REALX = REAL(N_CRN(NSL))
        DT = DT_CRN/REALX
        DTI = 1.D+0/(DT+EPSL)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CRN_LIM_BH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CRNTNB_BF
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
!     Compute the local maximum Courant numbers for borehole and
!     fracture flow.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 3 July 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNS_BH
      USE SOLTN
      USE PARM_BH
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CRNTNB_BF'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fracture triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
          CRNTL_FRC(NTX) = 0.D+0
          CRNTG_FRC(NTX) = 0.D+0
          CRNTN_FRC(NTX) = 0.D+0
!
!---      Loop over fracture triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            IF( IAQU.GT.0 ) CRNTL_FRC(NTX) = MAX( CRNTL_FRC(NTX),
     &        UFFL(1,NCX)*DT/DFF_FRC(NCX) )
            IF( IGAS.GT.0 ) CRNTG_FRC(NTX) = MAX( CRNTG_FRC(NTX),
     &        UFFG(1,NCX)*DT/DFF_FRC(NCX) )
            IF( INAPL.GT.0 ) CRNTN_FRC(NTX) = MAX( CRNTN_FRC(NTX),
     &        UFFN(1,NCX)*DT/DFF_FRC(NCX) )
          ENDDO
        ENDDO
      ENDDO
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
          CRNTL_BH(NBN) = 0.D+0
          CRNTG_BH(NBN) = 0.D+0
          CRNTN_BH(NBN) = 0.D+0
!
!---      Loop over borehole node to node connections ---
!
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.LE.0 ) CYCLE
            IF( IAQU.GT.0 ) CRNTL_BH(NBN) = MAX( CRNTL_BH(NBN),
     &        UBBL(1,NCX)*DT/DBB_BH(NCX) )
            IF( IGAS.GT.0 ) CRNTG_BH(NBN) = MAX( CRNTG_BH(NBN),
     &        UBBG(1,NCX)*DT/DBB_BH(NCX) )
            IF( INAPL.GT.0 ) CRNTN_BH(NBN) = MAX( CRNTN_BH(NBN),
     &        UBBN(1,NCX)*DT/DBB_BH(NCX) )
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CRNTNB_BF group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FTTIT ( V1X,V2X,P1X,P2X,PN1X,PN2X,ICPX,ITX )
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
!     FTTIT A fast triangle-triangle intersection test.
!
!     Tomas Moller "A fast triangle-triangle intersection test."
!
!     V1X - vertices of the first triangle (fracture)
!       1,n x-coordinate of n vertex
!       2,n y-coordinate of n vertex
!       3,n z-coordinate of n vertex
!     V2X - vertices of the second triangle (grid surface)
!       1,n x-coordinate of n vertex
!       2,n y-coordinate of n vertex
!       3,n z-coordinate of n vertex
!     PN1X - normal vector to fracture triangle
!     PN2X - normal vector to grid surface triangle
!     ITX - itersection test
!       0 - no intersection
!       1 - intersection
!     ICPX - coplanar index
!       0 - triangles are not coplanar
!       1 - triangles are coplanar
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 August 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 V1X(3,3),V2X(3,3)
      REAL*8 VP1X(3),VP2X(3)
      REAL*8 AX(3),BX(3),DX(3)
      REAL*8 PN1X(3),PN2X(3)
      REAL*8 P1X(3),P2X(3)
      REAL*8 DV1X(3),DV2X(3)
      REAL*8 SECT1(2),SECT2(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FTTIT'
!
!---  Initialize points  ---
!
      DO N = 1,3
        P1X(N) = 0.D+0
        P2X(N) = 0.D+0
      ENDDO
!
!---  Plane equation Pi1 : N1 dot X + D1 = 0, where X is any point
!     on the plane is computed
!     N1 = (V1(2)-V1(1)) cross (V1(3)-V1(1)) 
!     D1 = -N1 dot V1(1) ---
!
      CALL VSUB( V1X(1,2),V1X(1,1),AX )
      CALL VSUB( V1X(1,3),V1X(1,1),BX )
      DO N = 1,3
        IF( ABS(AX(N)).LT.1.D-12 ) AX(N) = 0.D+0 
        IF( ABS(BX(N)).LT.1.D-12 ) BX(N) = 0.D+0 
      ENDDO
      CALL VCROSSP( AX,BX,PN1X )
      DO N = 1,3
        IF( ABS(PN1X(N)).LT.1.D-12 ) PN1X(N) = 0.D+0
      ENDDO
      PD1X = -VDOTP( PN1X,V1X(1,1) )
!
!---  Put V2 into plane equation Pi1 to compute signed distances
!     to plane  ---
!
      DO N = 1,3
        DV2X(N) = VDOTP( PN1X,V2X(1,N) )
        DV2X(N) = DV2X(N) + PD1X
      ENDDO
!
!---  Coplanarity robustness check  ---
!
      DO N = 1,3
        IF( ABS(DV2X(N)).LT.1.D-12 ) DV2X(N) = 0.D+0
      ENDDO
      DV2AX = DV2X(1)*DV2X(2)
      DV2BX = DV2X(1)*DV2X(3)
!
!---  Same sign on all of them, plus not equal to zero, therefore
!     no intersecion occurs  ---
!
      IF( DV2AX.GT.0.D+0 .AND. DV2BX.GT.0.D+0 ) THEN
        ITX = 0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Plane equation Pi2 : N2 dot X + D2 = 0, where X is any point
!     on the plane is computed
!     N2 = (V2(2)-V2(1)) cross (V2(3)-V2(1)) 
!     D2 = -N2 dot V2(1) ---
!
      CALL VSUB( V2X(1,2),V2X(1,1),AX )
      CALL VSUB( V2X(1,3),V2X(1,1),BX )
      DO N = 1,3
        IF( ABS(AX(N)).LT.1.D-12 ) AX(N) = 0.D+0 
        IF( ABS(BX(N)).LT.1.D-12 ) BX(N) = 0.D+0 
      ENDDO
      CALL VCROSSP( AX,BX,PN2X )
      DO N = 1,3
        IF( ABS(PN2X(N)).LT.1.D-12 ) PN2X(N) = 0.D+0
      ENDDO
      PD2X = -VDOTP( PN2X,V2X(1,1) )
!
!---  Put V1 into plane equation Pi2 to compute signed distances
!     to plane  ---
!
      DO N = 1,3
        DV1X(N) = VDOTP( PN2X,V1X(1,N) )
        DV1X(N) = DV1X(N) + PD2X
      ENDDO
!
!---  Coplanarity robustness check  ---
!
      DO N = 1,3
        IF( ABS(DV1X(N)).LT.1.D-12 ) DV1X(N) = 0.D+0
      ENDDO
      DV1AX = DV1X(1)*DV1X(2)
      DV1BX = DV1X(1)*DV1X(3)
!
!---  Same sign on all of them, plus not equal to zero, therefore
!     no intersecion occurs  ---
!
      IF( DV1AX.GT.0.D+0 .AND. DV1BX.GT.0.D+0 ) THEN
        ITX = 0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compute direction of intersection line, L  ---
!
      CALL VCROSSP( PN1X,PN2X,DX )
!
!---  Compute and index to the largest component of DX  ---
!
      DMAX = ABS(DX(1))
      INDX = 1
      DO N = 2,3
        IF( ABS(DX(N)).GT.DMAX ) THEN
          DMAX = DX(N)
          INDX = N
        ENDIF
      ENDDO
!
!---  This is the simplified projection onto L  ---
!
      DO N = 1,3
        VP1X(N) = V1X(INDX,N)
        VP2X(N) = V2X(INDX,N)
      ENDDO
!
!---  Compute interval for triangle 1  ---
!
      CALL INTERVALS( VP1X,DV1X,DV1AX,DV1BX,A1X,B1X,C1X,X10X,X11X,ICPX )
!
!---  Triangles are coplanar  ---
!
      IF( ICPX.NE.0 ) THEN
        ITX = 0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compute interval for triangle 2  ---
!
      CALL INTERVALS( VP2X,DV2X,DV2AX,DV2BX,A2X,B2X,C2X,X20X,X21X,ICPX )
!
!---  Triangles are coplanar  ---
!
      IF( ICPX.NE.0 ) THEN
        ITX = 0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Check for overlap of triangle intervals  ---
!
      X1X = X10X*X11X
      X2X = X20X*X21X
      X12X = X1X*X2X
      VARX = A1X*X12X
      SECT1(1) = VARX + B1X*X11X*X2X
      SECT1(2) = VARX + C1X*X10X*X2X
      VARX = A2X*X12X
      SECT2(1) = VARX + B2X*X21X*X1X
      SECT2(2) = VARX + C2X*X20X*X1X
      CALL SORT_VAR( SECT1(1),SECT1(2) )
      CALL SORT_VAR( SECT2(1),SECT2(2) )
!
!---  No intersection occurs  ---
!
!      IF( SECT1(2).LT.SECT2(1) .OR. SECT2(2).LT.SECT1(1) ) THEN
      IF( (SECT2(1)-SECT1(2)).GT.EPSL .OR. 
     &  (SECT1(1)-SECT2(2)).GT.EPSL ) THEN
        ITX = 0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Intersection occurs, find four points along intersection line  ---
!
      CALL LINE_INT( DX,DV1X,DV2X,P1X,P2X,PN1X,PN2X,V1X,V2X )
      ITX = 1
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FTTIT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      FUNCTION ICCWTTI2D( CP1X,CQ1X,CR1X,CP2X,CR2X,CQ2X )
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
!     2D triangle triangle intersection test
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 March 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 CP1X(2),CQ1X(2),CR1X(2)
      REAL*8 CP2X(2),CQ2X(2),CR2X(2)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ORIENT_2D( CP2X,CQ2X,CP1X ).GE.0.D+0 ) THEN
        IF( ORIENT_2D( CQ2X,CR2X,CP1X ).GE.0.D+0 ) THEN
          IF( ORIENT_2D( CR2X,CP2X,CP1X ).GE.0.D+0 ) THEN
            ICCWTTI2D = 1
          ELSE
            ICCWTTI2D = ITEDGE( CP1X,CQ1X,CR1X,CP2X,CQ2X,CR2X )
          ENDIF
        ELSE
          IF( ORIENT_2D( CR2X,CP2X,CP1X ).GE.0.D+0 ) THEN
            ICCWTTI2D = ITEDGE( CP1X,CQ1X,CR1X,CR2X,CP2X,CQ2X )
          ELSE
            ICCWTTI2D = ITVERTEX( CP1X,CQ1X,CR1X,CP2X,CQ2X,CR2X )
          ENDIF
        ENDIF
      ELSE
        IF( ORIENT_2D( CQ2X,CR2X,CP1X ).GE.0.D+0 ) THEN
          IF( ORIENT_2D( CR2X,CP2X,CP1X ).GE.0.D+0 ) THEN
            ICCWTTI2D = ITEDGE( CP1X,CQ1X,CR1X,CQ2X,CR2X,CP2X )
          ELSE
            ICCWTTI2D = ITVERTEX( CP1X,CQ1X,CR1X,CQ2X,CR2X,CP2X )
          ENDIF
        ELSE
          ICCWTTI2D = ITVERTEX( CP1X,CQ1X,CR1X,CR2X,CP2X,CQ2X )
        ENDIF
      ENDIF
!
!---  End of ICCWTTI2D group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      FUNCTION ICPTT3D( P1X,Q1X,R1X,P2X,Q2X,R2X,VN1X,VN2X )
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
!     Coplanar triangles
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 March 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 P1X(3),Q1X(3),R1X(3)
      REAL*8 P2X(3),Q2X(3),R2X(3)
      REAL*8 CP1X(2),CQ1X(2),CR1X(2)
      REAL*8 CP2X(2),CQ2X(2),CR2X(2)
      REAL*8 VN1X(3),VN2X(3)
!
!----------------------Executable Lines--------------------------------!
!
      VN_X = ABS(VN1X(1))
      VN_Y = ABS(VN1X(2))
      VN_Z = ABS(VN1X(3))
!
!     Projection of the triangles in 3D onto 2D such that the area of
!     the projection is maximized.
!
      IF( VN_X.GT.VN_Z .AND. VN_X.GE.VN_Y ) THEN
!
!       Project onto plane yz
!
        CP1X(1) = Q1X(3)
        CP1X(2) = Q1X(2)
        CQ1X(1) = P1X(3)
        CQ1X(2) = P1X(2)
        CR1X(1) = R1X(3)
        CR1X(2) = R1X(2)
        CP2X(1) = Q2X(3)
        CP2X(2) = Q2X(2)
        CQ2X(1) = P2X(3)
        CQ2X(2) = P2X(2)
        CR2X(1) = R2X(3)
        CR2X(2) = R2X(2)
      ELSEIF( VN_Y.GT.VN_Z .AND. VN_Y.GE.VN_X ) THEN
!
!       Project onto plane xz
!
        CP1X(1) = Q1X(1)
        CP1X(2) = Q1X(3)
        CQ1X(1) = P1X(1)
        CQ1X(2) = P1X(3)
        CR1X(1) = R1X(1)
        CR1X(2) = R1X(3)
        CP2X(1) = Q2X(1)
        CP2X(2) = Q2X(3)
        CQ2X(1) = P2X(1)
        CQ2X(2) = P2X(3)
        CR2X(1) = R2X(1)
        CR2X(2) = R2X(3)
      ELSE
!
!       Project onto plane xy
!
        CP1X(1) = P1X(1)
        CP1X(2) = P1X(2)
        CQ1X(1) = Q1X(1)
        CQ1X(2) = Q1X(2)
        CR1X(1) = R1X(1)
        CR1X(2) = R1X(2)
        CP2X(1) = P2X(1)
        CP2X(2) = P2X(2)
        CQ2X(1) = Q2X(1)
        CQ2X(2) = Q2X(2)
        CR2X(1) = R2X(1)
        CR2X(2) = R2X(2)
      ENDIF
      ICPTT3D = ITTOT2D( CP1X,CQ1X,CR1X,CP2X,CQ2X,CR2X )
!
!---  End of ICPTT3D group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      FUNCTION ITTOT2D( CP1X,CQ1X,CR1X,CP2X,CQ2X,CR2X )
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
!     2D triangle triangle overlap test
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 March 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 CP1X(2),CQ1X(2),CR1X(2)
      REAL*8 CP2X(2),CQ2X(2),CR2X(2)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ORIENT_2D( CP1X,CQ1X,CR1X ).LT.0.D+0 ) THEN
        IF( ORIENT_2D( CP2X,CQ2X,CR2X ).LT.0.D+0 ) THEN
          ITTOT2D = ICCWTTI2D( CP1X,CR1X,CQ1X,CP2X,CR2X,CQ2X )
        ELSE
          ITTOT2D = ICCWTTI2D( CP1X,CR1X,CQ1X,CP2X,CQ2X,CR2X )
        ENDIF
      ELSE
        IF( ORIENT_2D( CP2X,CQ2X,CR2X ).LT.0.D+0 ) THEN
          ITTOT2D = ICCWTTI2D( CP1X,CQ1X,CR1X,CP2X,CR2X,CQ2X )
        ELSE
          ITTOT2D = ICCWTTI2D( CP1X,CQ1X,CR1X,CP2X,CQ2X,CR2X )
        ENDIF
      ENDIF
!
!---  End of ITTOT2D group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      FUNCTION ITEDGE( CP1X,CQ1X,CR1X,CP2X,CQ2X,CR2X )
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
!     2D edge intersection test
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 March 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 CP1X(2),CQ1X(2),CR1X(2)
      REAL*8 CP2X(2),CQ2X(2),CR2X(2)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ORIENT_2D( CR2X,CP2X,CQ1X ).GE.0.D+0 ) THEN
        IF( ORIENT_2D( CP1X,CP2X,CQ1X ).GE.0.D+0 ) THEN
          IF( ORIENT_2D( CP1X,CQ1X,CR2X ).GE.0.D+0 ) THEN
            ITEDGE = 1
          ELSE
            ITEDGE = 0
          ENDIF
        ELSE
          IF( ORIENT_2D( CQ1X,CR1X,CP2X ).GE.0.D+0 ) THEN
            IF( ORIENT_2D( CR1X,CP1X,CP2X ).GE.0.D+0 ) THEN
              ITEDGE = 1
            ELSE
              ITEDGE = 0
            ENDIF
          ELSE
            ITEDGE = 0
          ENDIF
        ENDIF
      ELSE
        IF( ORIENT_2D( CR2X,CP2X,CR1X ).GE.0.D+0 ) THEN
          IF( ORIENT_2D( CP1X,CP2X,CR1X ).GE.0.D+0 ) THEN
            IF( ORIENT_2D( CP1X,CR1X,CR2X ).GE.0.D+0 ) THEN
              ITEDGE = 1
            ELSE
              IF( ORIENT_2D( CQ1X,CR1X,CR2X ).GE.0.D+0 ) THEN
                ITEDGE = 1
              ELSE
                ITEDGE = 0
              ENDIF
            ENDIF
          ELSE
            ITEDGE = 0
          ENDIF
        ELSE
          ITEDGE = 0
        ENDIF
      ENDIF
!
!---  End of ITEDGE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      FUNCTION ITVERTEX( CP1X,CQ1X,CR1X,CP2X,CQ2X,CR2X )
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
!     2D vertex intersection test
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 March 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 CP1X(2),CQ1X(2),CR1X(2)
      REAL*8 CP2X(2),CQ2X(2),CR2X(2)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ORIENT_2D( CR2X,CP2X,CQ1X ).GE.0.D+0 ) THEN
        IF( ORIENT_2D( CR2X,CQ2X,CQ1X ).LE.0.D+0 ) THEN
          IF( ORIENT_2D( CP1X,CP2X,CQ1X ).GT.0.D+0 ) THEN
            IF( ORIENT_2D( CP1X,CQ2X,CQ1X ).LE.0.D+0 ) THEN
              ITVERTEX = 1
            ELSE
              ITVERTEX = 0
            ENDIF
          ELSE
            IF( ORIENT_2D( CP1X,CP2X,CR1X ).GE.0.D+0 ) THEN
              IF( ORIENT_2D( CQ1X,CR1X,CP2X ).GE.0.D+0 ) THEN
                ITVERTEX = 1
              ELSE
                ITVERTEX = 0
              ENDIF
            ELSE
              ITVERTEX = 0
            ENDIF
          ENDIF
        ELSE
          IF( ORIENT_2D( CP1X,CQ2X,CQ1X ).LE.0.D+0 ) THEN
            IF( ORIENT_2D( CR2X,CQ2X,CR1X ).LE.0.D+0 ) THEN
              IF( ORIENT_2D( CQ1X,CR1X,CQ2X ).GE.0.D+0 ) THEN
                ITVERTEX = 1
              ELSE
                ITVERTEX = 0
              ENDIF
            ELSE
              ITVERTEX = 0
            ENDIF
          ELSE
            ITVERTEX = 0
          ENDIF
        ENDIF
      ELSE
        IF( ORIENT_2D( CR2X,CP2X,CR1X ).GE.0.D+0 ) THEN
          IF( ORIENT_2D( CQ1X,CR1X,CR2X ).GE.0.D+0 ) THEN
            IF( ORIENT_2D( CP1X,CP2X,CR1X ).GE.0.D+0 ) THEN
              ITVERTEX = 1
            ELSE
              ITVERTEX = 0
            ENDIF
          ELSE
            IF( ORIENT_2D( CQ1X,CR1X,CQ2X ).GE.0.D+0 ) THEN
              IF( ORIENT_2D( CR2X,CR1X,CQ2X ).GE.0.D+0 ) THEN
                ITVERTEX = 1
              ELSE
                ITVERTEX = 0
              ENDIF
            ELSE
              ITVERTEX = 0
            ENDIF
          ENDIF
        ENDIF
      ENDIF      
!
!---  End of ITVERTEX group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      FUNCTION ORIENT_2D( AX,BX,CX )
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
!     2D Orientation
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 March 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AX(2),BX(2),CX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ORIENT_2D'
      ORIENT_2D = ((AX(1)-CX(1))*(BX(2)-CX(2)) - 
     &  (AX(2)-CX(2))*(BX(1)-CX(1)))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ORIENT_2D group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TTIT3D( V1X,V2X,PE1X,PE2X,VN1X,VN2X,CPX,ICPX,ITX )
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
!     TTIT3D Triangle-triangle intersection test
!
!     This file contains C implementation of algorithms for                
!     performing two and three-dimensional triangle-triangle 
!     intersection test 
!     The algorithms and underlying theory are described in                    
!                                                                           
!     "Fast and Robust Triangle-Triangle Overlap Test 
!     Using Orientation Predicates"  P. Guigue - O. Devillers
!                                                 
!     Journal of Graphics Tools, 8(1), 2003                                    
!
!     V1X - vertices of the first triangle (fracture)
!       1,n x-coordinate of n vertex
!       2,n y-coordinate of n vertex
!       3,n z-coordinate of n vertex
!     V2X - vertices of the second triangle (grid surface)
!       1,n x-coordinate of n vertex
!       2,n y-coordinate of n vertex
!       3,n z-coordinate of n vertex
!     PE1X - first end point of the line segment of intersection
!     PE2X - second end point of the line segment of intersection
!     ITX - itersection test
!       0 - no intersection
!       1 - intersection
!     ICPX - coplanar index
!       0 - triangles are not coplanar
!       1 - triangles are coplanar
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 8 March 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 V1X(3,3),V2X(3,3)
      REAL*8 P1X(3),Q1X(3),R1X(3)
      REAL*8 P2X(3),Q2X(3),R2X(3)
      REAL*8 VV1X(3),VV2X(3)
      REAL*8 VN1X(3),VN2X(3)
      REAL*8 PE1X(3),PE2X(3)
      REAL*8 CPX(3,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TTIT3D'
!
!     Load local triangle vectors
!
      DO M = 1,3
        P1X(M) = V1X(M,1)
        Q1X(M) = V1X(M,2)
        R1X(M) = V1X(M,3)
        P2X(M) = V2X(M,1)
        Q2X(M) = V2X(M,2)
        R2X(M) = V2X(M,3)
      ENDDO
!
!     Compute distance signs  of p1, q1 and r1 
!     to the plane of triangle(p2,q2,r2)
!
      CALL VSUB( P2X,R2X,VV1X )
      CALL VSUB( Q2X,R2X,VV2X )
      CALL VCROSSP( VV1X,VV2X,VN2X )
      CALL VSUB( P1X,R2X,VV1X )
      DP1X = VDOTP( VV1X,VN2X )
      CALL VSUB( Q1X,R2X,VV1X )
      DQ1X = VDOTP( VV1X,VN2X )
      CALL VSUB( R1X,R2X,VV1X )
      DR1X = VDOTP( VV1X,VN2X )
      IF( ((DP1X*DQ1X).GT.0.D+0).AND.((DP1X*DR1X).GT.0.D+0) ) THEN
        ITX = 0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!     Compute distance signs  of p2, q2 and r2 
!     to the plane of triangle(p1,q1,r1)
!
      CALL VSUB( P1X,R1X,VV1X )
      CALL VSUB( Q1X,R1X,VV2X )
      CALL VCROSSP( VV1X,VV2X,VN1X )
      CALL VSUB( P2X,R1X,VV1X )
      DP2X = VDOTP( VV1X,VN1X )
      CALL VSUB( Q2X,R1X,VV1X )
      DQ2X = VDOTP( VV1X,VN1X )
      CALL VSUB( R2X,R1X,VV1X )
      DR2X = VDOTP( VV1X,VN1X )
      IF( ((DP2X*DQ2X).GT.0.D+0).AND.((DP2X*DR2X).GT.0.D+0) ) THEN
        ITX = 0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!     Permutation in a canonical form of T1's vertices
!
      IF( DP1X.GT.0.D+0 ) THEN
        IF( DQ1X.GT.0.D+0 ) THEN
          CALL TTI3D( R1X,P1X,Q1X,P2X,R2X,Q2X,DP2X,DR2X,DQ2X,
     &      VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
        ELSEIF( DR1X.GT.0.D+0 ) THEN
          CALL TTI3D( Q1X,R1X,P1X,P2X,R2X,Q2X,DP2X,DR2X,DQ2X,
     &      VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
        ELSE
          CALL TTI3D( P1X,Q1X,R1X,P2X,Q2X,R2X,DP2X,DQ2X,DR2X,
     &      VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
        ENDIF
      ELSEIF( DP1X.LT.0.D+0 ) THEN
        IF( DQ1X.LT.0.D+0 ) THEN
          CALL TTI3D( R1X,P1X,Q1X,P2X,Q2X,R2X,DP2X,DQ2X,DR2X,
     &      VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
        ELSEIF( DR1X.LT.0.D+0 ) THEN
          CALL TTI3D( Q1X,R1X,P1X,P2X,Q2X,R2X,DP2X,DQ2X,DR2X,
     &      VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
        ELSE
          CALL TTI3D( P1X,Q1X,R1X,P2X,R2X,Q2X,DP2X,DR2X,DQ2X,
     &      VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
        ENDIF
      ELSE
        IF( DQ1X.LT.0.D+0 ) THEN
          IF( DR1X.GE.0.D+0 ) THEN
            CALL TTI3D( Q1X,R1X,P1X,P2X,R2X,Q2X,DP2X,DR2X,DQ2X,
     &        VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
          ELSE
            CALL TTI3D( P1X,Q1X,R1X,P2X,Q2X,R2X,DP2X,DQ2X,DR2X,
     &        VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
          ENDIF
        ELSEIF( DQ1X.GT.0.D+0 ) THEN
          IF( DR1X.GT.0.D+0 ) THEN
            CALL TTI3D( P1X,Q1X,R1X,P2X,R2X,Q2X,DP2X,DR2X,DQ2X,
     &        VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
          ELSE
            CALL TTI3D( Q1X,R1X,P1X,P2X,Q2X,R2X,DP2X,DQ2X,DR2X,
     &        VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
          ENDIF
        ELSE
          IF( DR1X.GT.0.D+0 ) THEN
            CALL TTI3D( R1X,P1X,Q1X,P2X,Q2X,R2X,DP2X,DQ2X,DR2X,
     &        VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
          ELSEIF( DR1X.LT.0.D+0 ) THEN
            CALL TTI3D( R1X,P1X,Q1X,P2X,R2X,Q2X,DP2X,DR2X,DQ2X,
     &        VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
!
!         Triangles are co-planar
!
          ELSE
            ITX = 0
            CALL CPTT3D( P1X,Q1X,R1X,P2X,Q2X,R2X,VN1X,VN2X,CPX,ICPX )
          ENDIF
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TTIT3D group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TTI3D( P1X,Q1X,R1X,P2X,Q2X,R2X,DP2X,DQ2X,DR2X,
     &  VN1X,VN2X,PE1X,PE2X,ICPX,ITX )
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
!     TTI3D Triangle-triangle intersection
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 March 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 P1X(3),Q1X(3),R1X(3)
      REAL*8 P2X(3),Q2X(3),R2X(3)
      REAL*8 PE1X(3),PE2X(3)
      REAL*8 VN1X(3),VN2X(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TTI3D'
      IF( DP2X.GT.0.D+0 ) THEN
        IF( DQ2X.GT.0.D+0 ) THEN
          CALL CONST_INTER( P1X,R1X,Q1X,R2X,P2X,Q2X,VN1X,VN2X,
     &      PE1X,PE2X,ITX )
        ELSEIF( DR2X.GT.0.D+0 ) THEN
          CALL CONST_INTER( P1X,R1X,Q1X,Q2X,R2X,P2X,VN1X,VN2X,
     &      PE1X,PE2X,ITX )
        ELSE
          CALL CONST_INTER( P1X,Q1X,R1X,P2X,Q2X,R2X,VN1X,VN2X,
     &      PE1X,PE2X,ITX )
        ENDIF
      ELSEIF( DP2X.LT.0.D+0 ) THEN
        IF( DQ2X.LT.0.D+0 ) THEN
          CALL CONST_INTER( P1X,Q1X,R1X,R2X,P2X,Q2X,VN1X,VN2X,
     &      PE1X,PE2X,ITX )
        ELSEIF( DR2X.LT.0.D+0 ) THEN
          CALL CONST_INTER( P1X,Q1X,R1X,Q2X,R2X,P2X,VN1X,VN2X,
     &      PE1X,PE2X,ITX )
        ELSE
          CALL CONST_INTER( P1X,R1X,Q1X,P2X,Q2X,R2X,VN1X,VN2X,
     &      PE1X,PE2X,ITX )
        ENDIF
      ELSE
        IF( DQ2X.LT.0.D+0 ) THEN
          IF( DR2X.GE.0.D+0 ) THEN
            CALL CONST_INTER( P1X,R1X,Q1X,Q2X,R2X,P2X,VN1X,VN2X,
     &        PE1X,PE2X,ITX )
          ELSE
            CALL CONST_INTER( P1X,Q1X,R1X,P2X,Q2X,R2X,VN1X,VN2X,
     &        PE1X,PE2X,ITX )
          ENDIF
        ELSEIF( DQ2X.GT.0.D+0 ) THEN
          IF( DR2X.GT.0.D+0 ) THEN
            CALL CONST_INTER( P1X,R1X,Q1X,P2X,Q2X,R2X,VN1X,VN2X,
     &        PE1X,PE2X,ITX )
          ELSE
            CALL CONST_INTER( P1X,Q1X,R1X,Q2X,R2X,P2X,VN1X,VN2X,
     &        PE1X,PE2X,ITX )
          ENDIF
        ELSE
          IF( DR2X.GT.0.D+0 ) THEN
            CALL CONST_INTER( P1X,Q1X,R1X,R2X,P2X,Q2X,VN1X,VN2X,
     &        PE1X,PE2X,ITX )
          ELSEIF( DR2X.LT.0.D+0 ) THEN
            CALL CONST_INTER( P1X,R1X,Q1X,R2X,P2X,Q2X,VN1X,VN2X,
     &        PE1X,PE2X,ITX )
          ELSE
            ICPX = ICPTT3D( P1X,Q1X,R1X,P2X,Q2X,R2X,VN1X,VN2X )
          ENDIF
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TTI3D group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTERVALS( VPX,DVX,DVAX,DVBX,AX,BX,CX,X0X,X1X,ICPX )
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
!     Identify triangle intersection intervals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 August 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VPX(3),DVX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTERVALS'
!
!---  Here we know that DVBX <= 0.0, that is DVX(1) and DVX(2)
!     are on the same side, DVX(3) is on the other side or on
!     the plane  ---
!
      IF( DVAX.GT.0.D+0 ) THEN
        AX = VPX(3)
        BX = (VPX(1)-VPX(3))*DVX(3)
        CX = (VPX(2)-VPX(3))*DVX(3)
        X0X = DVX(3)-DVX(1)
        X1X = DVX(3)-DVX(2)
!
!---  Here we know that DVAX <= 0.0  ---
!
      ELSEIF( DVBX.GT.0.D+0 ) THEN
        AX = VPX(2)
        BX = (VPX(1)-VPX(2))*DVX(2)
        CX = (VPX(3)-VPX(2))*DVX(2)
        X0X = DVX(2)-DVX(1)
        X1X = DVX(2)-DVX(3)
!
!---  Here we know that DVX(2)*DVX(3) > 0.0 || DVX(1) /= 0.0  ---
!
      ELSEIF( DVX(2)*DVX(3).GT.0.D+0 .OR. 
     &  ABS(DVX(1))/EPSL.GT.EPSL ) THEN
        AX = VPX(1)
        BX = (VPX(2)-VPX(1))*DVX(1)
        CX = (VPX(3)-VPX(1))*DVX(1)
        X0X = DVX(1)-DVX(2)
        X1X = DVX(1)-DVX(3)
!
!---  Here we know that DVX(2) /= 0.0  ---
!
      ELSEIF( ABS(DVX(2))/EPSL.GT.EPSL ) THEN
        AX = VPX(2)
        BX = (VPX(1)-VPX(2))*DVX(2)
        CX = (VPX(3)-VPX(2))*DVX(2)
        X0X = DVX(2)-DVX(1)
        X1X = DVX(2)-DVX(3)
!
!---  Here we know that DVX(3) /= 0.0  ---
!
      ELSEIF( ABS(DVX(3))/EPSL.GT.EPSL ) THEN
        AX = VPX(3)
        BX = (VPX(1)-VPX(3))*DVX(3)
        CX = (VPX(2)-VPX(3))*DVX(3)
        X0X = DVX(3)-DVX(1)
        X1X = DVX(3)-DVX(2)
!
!---  Triangles are coplanar  ---
!
      ELSE
        ICPX = ICPX + 1
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTERVALS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LINE_INT( DX,DV1X,DV2X,P1X,P2X,PN1X,PN2X,V1X,V2X )
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
!     Find the x,y,z coordinates of the four intersection points of 
!     two triangles along the intersection line.
!
!     DX - direction vector of intersection line
!     DV1X - signed distances from triangle 1 to plane Pi2
!     DV2X - signed distances from triangle 2 to plane Pi1
!     PN1X - normal vector to triangle 1 or plane Pi1
!     PN2X - normal vector to triangle 2 or plane Pi2
!     P1X - x,y,z coordinates of point 1 (triangle 1)
!     P2X - x,y,z coordinates of point 2 (triangle 1)
!     P3X - x,y,z coordinates of point 1 (triangle 2)
!     P4X - x,y,z coordinates of point 2 (triangle 2)
!     V1X - x,y,z coordinates of vertices of triangle 1
!     V2X - x,y,z coordinates of vertices of triangle 2
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 August 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 DX(3),PN1X(3),PN2X(3)
      REAL*8 DV1X(3),DV2X(3)
      REAL*8 P1X(3),P2X(3),P3X(3),P4X(3)
      REAL*8 POX(3),PNX(3)
      REAL*8 AX(3),BX(3),CX(3),CABX(3),CCBX(3)
      REAL*8 V1X(3,3),V2X(3,3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LINE_INT'
!
!---  Pi1 and Pi2 plane equations  ---
!
      A1X = PN1X(1)
      B1X = PN1X(2)
      C1X = PN1X(3)
      D1X = A1X*V1X(1,3) + B1X*V1X(2,3) + C1X*V1X(3,3)
      A2X = PN2X(1)
      B2X = PN2X(2)
      C2X = PN2X(3)
      D2X = A2X*V2X(1,3) + B2X*V2X(2,3) + C2X*V2X(3,3)
!
!---  Intersection line of plane Pi1 and Pi2  ---
!
      IF( ABS(A1X).GT.EPSL .AND. ABS(B2X-B1X*A2X/A1X).GT.EPSL ) THEN
        POX(3) = (V1X(3,1)+V1X(3,2)+V1X(3,3)+V2X(3,1)+V2X(3,2)+V2X(3,3))
     &    /6.D+0
        POX(2) = (POX(3)*(C1X*A2X/A1X - C2X) + (D2X - D1X*A2X/A1X))/
     &    (B2X - B1X*A2X/A1X)
        POX(1) = (D1X - B1X*POX(2) - C1X*POX(3))/A1X
        PNX(1) = POX(1) + DX(1)
        PNX(2) = POX(2) + DX(2)
        PNX(3) = POX(3) + DX(3)
      ELSEIF( ABS(B1X).GT.EPSL .AND. ABS(C2X-C1X*B2X/B1X).GT.EPSL ) THEN
        POX(1) = (V1X(1,1)+V1X(1,2)+V1X(1,3)+V2X(1,1)+V2X(1,2)+V2X(1,3))
     &    /6.D+0
        POX(3) = (POX(1)*(A1X*B2X/B1X - A2X) + (D2X - D1X*B2X/B1X))/
     &    (C2X - C1X*B2X/B1X)
        POX(2) = (D1X - A1X*POX(1) - C1X*POX(3))/B1X
        PNX(1) = POX(1) + DX(1)
        PNX(2) = POX(2) + DX(2)
        PNX(3) = POX(3) + DX(3)
      ELSEIF( ABS(C1X).GT.EPSL .AND. ABS(A2X-A1X*C2X/C1X).GT.EPSL ) THEN
        POX(2) = (V1X(2,1)+V1X(2,2)+V1X(2,3)+V2X(2,1)+V2X(2,2)+V2X(2,3))
     &    /6.D+0
        POX(1) = (POX(2)*(B1X*C2X/C1X - B2X) + (D2X - D1X*C2X/C1X))/
     &    (A2X - A1X*C2X/C1X)
        POX(3) = (D1X - A1X*POX(1) - B1X*POX(2))/C1X
        PNX(1) = POX(1) + DX(1)
        PNX(2) = POX(2) + DX(2)
        PNX(3) = POX(3) + DX(3)
      ELSEIF( ABS(A2X).GT.EPSL .AND. ABS(B1X-B2X*A1X/A2X).GT.EPSL ) THEN
        POX(3) = (V2X(3,1)+V2X(3,2)+V2X(3,3)+V1X(3,1)+V1X(3,2)+V1X(3,3))
     &    /6.D+0
        POX(2) = (POX(3)*(C2X*A1X/A2X - C1X) + (D1X - D2X*A1X/A2X))/
     &    (B1X - B2X*A1X/A2X)
        POX(1) = (D2X - B2X*POX(2) - C2X*POX(3))/A2X
        PNX(1) = POX(1) + DX(1)
        PNX(2) = POX(2) + DX(2)
        PNX(3) = POX(3) + DX(3)
      ELSEIF( ABS(B2X).GT.EPSL .AND. ABS(C1X-C2X*B1X/B2X).GT.EPSL ) THEN
        POX(1) = (V2X(1,1)+V2X(1,2)+V2X(1,3)+V1X(1,1)+V1X(1,2)+V1X(1,3))
     &    /6.D+0
        POX(3) = (POX(1)*(A2X*B1X/B2X - A1X) + (D1X - D2X*B1X/B2X))/
     &    (C1X - C2X*B1X/B2X)
        POX(2) = (D2X - A2X*POX(1) - C2X*POX(3))/B2X
        PNX(1) = POX(1) + DX(1)
        PNX(2) = POX(2) + DX(2)
        PNX(3) = POX(3) + DX(3)
      ELSEIF( ABS(C2X).GT.EPSL .AND. ABS(A1X-A2X*C1X/C2X).GT.EPSL ) THEN
        POX(2) = (V2X(2,1)+V2X(2,2)+V2X(2,3)+V1X(2,1)+V1X(2,2)+V1X(2,3))
     &    /6.D+0
        POX(1) = (POX(2)*(B2X*C1X/C2X - B1X) + (D1X - D2X*C1X/C2X))/
     &    (A1X - A2X*C1X/C2X)
        POX(3) = (D2X - A2X*POX(1) - B2X*POX(2))/C2X
        PNX(1) = POX(1) + DX(1)
        PNX(2) = POX(2) + DX(2)
        PNX(3) = POX(3) + DX(3)
      ENDIF
!
!---  Vertices 1 and 2 of triangle 1 on the same side of plane Pi2, 
!     P1X and P2X are the fracture triangle points  ---
!
      IF( SIGN(1.D+0,DV1X(1)).EQ.SIGN(1.D+0,DV1X(2)) ) THEN
        DO M = 1,3
          AX(M) = V1X(M,1)-V1X(M,3)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V1X(M,3)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P1X(M) = V1X(M,3) + AX(M)*SX
        ENDDO
        T1X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T1X = T1X + (P1X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T1X = T1X/REAL(NC)
        DO M = 1,3
          AX(M) = V1X(M,2)-V1X(M,3)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V1X(M,3)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P2X(M) = V1X(M,3) + AX(M)*SX
        ENDDO
        T2X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T2X = T2X + (P2X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T2X = T2X/REAL(NC)
!
!---  Vertices 2 and 3 of triangle 1 on the same side of plane Pi2, 
!     P1X and P2X are the fracture triangle points  ---  ---
!
      ELSEIF( SIGN(1.D+0,DV1X(2)).EQ.SIGN(1.D+0,DV1X(3)) ) THEN
        DO M = 1,3
          AX(M) = V1X(M,3)-V1X(M,1)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V1X(M,1)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P1X(M) = V1X(M,1) + AX(M)*SX
        ENDDO
        T1X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T1X = T1X + (P1X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T1X = T1X/REAL(NC)
        DO M = 1,3
          AX(M) = V1X(M,2)-V1X(M,1)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V1X(M,1)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P2X(M) = V1X(M,1) + AX(M)*SX
        ENDDO
        T2X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T2X = T2X + (P2X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T2X = T2X/REAL(NC)
!
!---  Vertices 3 and 1 of triangle 1 on the same side of plane Pi2, 
!     P1X and P2X are the fracture triangle points  ---  ---
!
      ELSE
        DO M = 1,3
          AX(M) = V1X(M,1)-V1X(M,2)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V1X(M,2)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P1X(M) = V1X(M,2) + AX(M)*SX
        ENDDO
        T1X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T1X = T1X + (P1X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T1X = T1X/REAL(NC)
        DO M = 1,3
          AX(M) = V1X(M,3)-V1X(M,2)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V1X(M,2)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P2X(M) = V1X(M,2) + AX(M)*SX
        ENDDO
        T2X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T2X = T2X + (P2X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T2X = T2X/REAL(NC)
      ENDIF
!
!---  Vertices 1 and 2 of triangle 2 on the same side of plane Pi1, 
!     P3X and P4X are the grid-cell surface triangle points  ---
!
      IF( SIGN(1.D+0,DV2X(1)).EQ.SIGN(1.D+0,DV2X(2)) ) THEN
        DO M = 1,3
          AX(M) = V2X(M,1)-V2X(M,3)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V2X(M,3)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P3X(M) = V2X(M,3) + AX(M)*SX
        ENDDO
        T3X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T3X = T3X + (P3X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T3X = T3X/REAL(NC)
        DO M = 1,3
          AX(M) = V2X(M,2)-V2X(M,3)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V2X(M,3)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P4X(M) = V2X(M,3) + AX(M)*SX
        ENDDO
        T4X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T4X = T4X + (P4X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T4X = T4X/REAL(NC)
!
!---  Vertices 2 and 3 of triangle 2 on the same side of plane Pi1, 
!     P3X and P4X are the grid-cell surface triangle points  ---
!
      ELSEIF( SIGN(1.D+0,DV2X(2)).EQ.SIGN(1.D+0,DV2X(3)) ) THEN
        DO M = 1,3
          AX(M) = V2X(M,3)-V2X(M,1)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V2X(M,1)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P3X(M) = V2X(M,1) + AX(M)*SX
        ENDDO
        T3X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T3X = T3X + (P3X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T3X = T3X/REAL(NC)
        DO M = 1,3
          AX(M) = V2X(M,2)-V2X(M,1)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V2X(M,1)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P4X(M) = V2X(M,1) + AX(M)*SX
        ENDDO
        T4X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T4X = T4X + (P4X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T4X = T4X/REAL(NC)
!
!---  Vertices 3 and 1 of triangle 2 on the same side of plane Pi1, 
!     P3X and P4X are the grid-cell surface triangle points  ---
!
      ELSE
        DO M = 1,3
          AX(M) = V2X(M,1)-V2X(M,2)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V2X(M,2)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P3X(M) = V2X(M,2) + AX(M)*SX
        ENDDO
        T3X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T3X = T3X + (P3X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T3X = T3X/REAL(NC)
        DO M = 1,3
          AX(M) = V2X(M,3)-V2X(M,2)
          BX(M) = PNX(M)-POX(M)
          CX(M) = POX(M)-V2X(M,2)
        ENDDO
        CALL VCROSSP( AX,BX,CABX )
        CALL VCROSSP( CX,BX,CCBX )
        SX = VDOTP( CABX,CCBX )/(VDOTP( CABX,CABX ) + SMALL)
        DO M = 1,3
          P4X(M) = V2X(M,2) + AX(M)*SX
        ENDDO
        T4X = 0.D+0
        NC = 0
        DO M = 1,3
          IF( ABS(DX(M)).GT.EPSL ) THEN
            NC = NC + 1
            T4X = T4X + (P4X(M)-POX(M))/DX(M)
          ENDIF
        ENDDO
        T4X = T4X/REAL(NC)
      ENDIF
!
!---  Limit the intersection of the fracture triangle and the grid-cell
!     surface triangle to the bounds of the later  ---
!
      TLX = MAX( MIN( T3X,T4X ),MIN( T1X,T2X ) )
      TUX = MIN( MAX( T3X,T4X ),MAX( T1X,T2X ) )
!
!---  Intersection points on the intersection line, limited by
!     the fracture and grid-cell surface triangle extents  ---
!
      P1X(1) = POX(1) + DX(1)*TLX
      P1X(2) = POX(2) + DX(2)*TLX
      P1X(3) = POX(3) + DX(3)*TLX
      P2X(1) = POX(1) + DX(1)*TUX
      P2X(2) = POX(2) + DX(2)*TUX
      P2X(3) = POX(3) + DX(3)*TUX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LINE_INT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LSEG_CHK( P1X,P2X,P3X,P4X,I0X,I1X,ITX )
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
!     Check for the intersection of two line segments.
!     http://paulbourke.net/geometry/pointlineplane/pdb.c
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 August 2017.
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
      REAL*8 P1X(3),P2X(3),P3X(3),P4X(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LSEG_CHK'
      AX = P1X(I0X)-P2X(I0X)
      AY = P1X(I1X)-P2X(I1X)
      BX = P3X(I0X)-P4X(I0X)
      BY = P3X(I1X)-P4X(I1X)
      CX = P2X(I0X)-P3X(I0X)
      CY = P2X(I1X)-P3X(I1X)
      FX = AY*BX - AX*BY
      DX = BY*CX - BX*CY
      IF( (FX.GT.0.D+0 .AND. DX.GE.0.D+0 .AND. DX.LE.FX) .OR.
     &  (FX.LT.0.D+0 .AND. DX.LE.0.D+0 .AND. DX.GE.FX) ) THEN
        EX = AX*CY - AY*CX
        IF( FX.GT.0.D+0 ) THEN
          IF( EX.GE.0.D+0 .AND. EX.LE.FX ) ITX = 1
        ELSE
          IF( EX.LE.0.D+0 .AND. EX.GE.FX ) ITX = 1
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LSEG_CHK group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LSEG_INT( P1X,P2X,P3X,P4X,PIX,NPIX )
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
!     Determine the intersection point(s) of two line segments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 P1X(2),P2X(2),P3X(2),P4X(2)
      REAL*8 PIX(2,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LSEG_INT'
      DENOMX = (P4X(2)-P3X(2))*(P2X(1)-P1X(1)) - 
     &  (P4X(1)-P3X(1))*(P2X(2)-P1X(2))
      ANUMERX = (P4X(1)-P3X(1))*(P1X(2)-P3X(2)) - 
     &  (P4X(2)-P3X(2))*(P1X(1)-P3X(1))
      BNUMERX = (P2X(1)-P1X(1))*(P1X(2)-P3X(2)) - 
     &  (P2X(2)-P1X(2))*(P1X(1)-P3X(1))
!
!---  Line segments are coincidental  ---
!
      IF( ABS(ANUMERX).LT.EPSL .AND. ABS(BNUMERX).LT.EPSL. AND.
     &  ABS(DENOMX).LT.EPSL ) THEN
        NPIX = NPIX + 1
        DO M = 1,2
          PIX(M,NPIX) = P1X(M)
        ENDDO
        NPIX = NPIX + 1
        DO M = 1,2
          PIX(M,NPIX) = P2X(M)
        ENDDO
!
!---  Line segments are parallel, no intersection  ---
!
      ELSEIF( ABS(DENOMX).LT.EPSL ) THEN
        NPIX = NPIX + 0
!
!---  Check for intersection along line segments  ---
!
      ELSE
        AMUX = ANUMERX/DENOMX
        BMUX = BNUMERX/DENOMX
!
!---    No intersection along line segments  ---
!
        IF( AMUX.LT.0.D+0 .OR. AMUX.GT.1.D+0 .OR. BMUX.LT.0.D+0 
     &    .OR. BMUX.GT.1.D+0 ) THEN
          NPIX = NPIX + 0
        ELSE
          NPIX = NPIX + 1
          DO M = 1,2
            PIX(M,NPIX) = P1X(M) + AMUX*(P2X(M)-P1X(M))
          ENDDO
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LSEG_INT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROJ_BH( PLX_WX,PLY_WX,PLZ_WX,XP_WX,YP_WX,ZP_WX,N )
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
!
!     STOMP-CO2e
!
!     Well projections on local grid.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 01 October 2012.
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
      REAL*8 ROT_MAT(3,3),PLGCX(3),PLLCX(3)
      REAL*8 XVX(4),YVX(4),ZVX(4)
      REAL*8 XP_WX(2),YP_WX(2),ZP_WX(2)
      REAL*8 PLX_WX,PLY_WX,PLZ_WX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROJ_BH'
!
!---  Node surfaces ---
!
      NPX = NSX(N)
      NPY = NSY(N)
      NPZ = NSZ(N)
      NQX = NSX(N)+1
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
      NQY = NSY(N)+IFLD
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
      NQZ = NSZ(N)+IJFLD
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
!
!---  West surface centroid  ---
!
      XVX(1) = XE(1,N)
      XVX(2) = XE(3,N)
      XVX(3) = XE(7,N)
      XVX(4) = XE(5,N)
      YVX(1) = YE(1,N)
      YVX(2) = YE(3,N)
      YVX(3) = YE(7,N)
      YVX(4) = YE(5,N)
      ZVX(1) = ZE(1,N)
      ZVX(2) = ZE(3,N)
      ZVX(3) = ZE(7,N)
      ZVX(4) = ZE(5,N)
      NP = 4
      CALL PGCNTRD( NP,XVX,YVX,ZVX,XFW,YFW,ZFW )
!
!---  East surface centroid  ---
!
      XVX(1) = XE(2,N)
      XVX(2) = XE(4,N)
      XVX(3) = XE(8,N)
      XVX(4) = XE(6,N)
      YVX(1) = YE(2,N)
      YVX(2) = YE(4,N)
      YVX(3) = YE(8,N)
      YVX(4) = YE(6,N)
      ZVX(1) = ZE(2,N)
      ZVX(2) = ZE(4,N)
      ZVX(3) = ZE(8,N)
      ZVX(4) = ZE(6,N)
      NP = 4
      CALL PGCNTRD( NP,XVX,YVX,ZVX,XFE,YFE,ZFE )
!
!---  South surface centroid  ---
!
      XVX(1) = XE(1,N)
      XVX(2) = XE(2,N)
      XVX(3) = XE(6,N)
      XVX(4) = XE(5,N)
      YVX(1) = YE(1,N)
      YVX(2) = YE(2,N)
      YVX(3) = YE(6,N)
      YVX(4) = YE(5,N)
      ZVX(1) = ZE(1,N)
      ZVX(2) = ZE(2,N)
      ZVX(3) = ZE(6,N)
      ZVX(4) = ZE(5,N)
      NP = 4
      CALL PGCNTRD( NP,XVX,YVX,ZVX,XFS,YFS,ZFS )
!
!---  North surface centroid  ---
!
      XVX(1) = XE(3,N)
      XVX(2) = XE(4,N)
      XVX(3) = XE(8,N)
      XVX(4) = XE(7,N)
      YVX(1) = YE(3,N)
      YVX(2) = YE(4,N)
      YVX(3) = YE(8,N)
      YVX(4) = YE(7,N)
      ZVX(1) = ZE(3,N)
      ZVX(2) = ZE(4,N)
      ZVX(3) = ZE(8,N)
      ZVX(4) = ZE(7,N)
      NP = 4
      CALL PGCNTRD( NP,XVX,YVX,ZVX,XFN,YFN,ZFN )
!
!---  Bottom surface centroid  ---
!
      XVX(1) = XE(1,N)
      XVX(2) = XE(2,N)
      XVX(3) = XE(4,N)
      XVX(4) = XE(3,N)
      YVX(1) = YE(1,N)
      YVX(2) = YE(2,N)
      YVX(3) = YE(4,N)
      YVX(4) = YE(3,N)
      ZVX(1) = ZE(1,N)
      ZVX(2) = ZE(2,N)
      ZVX(3) = ZE(4,N)
      ZVX(4) = ZE(3,N)
      NP = 4
      CALL PGCNTRD( NP,XVX,YVX,ZVX,XFB,YFB,ZFB )
!
!---  Top surface centroid  ---
!
      XVX(1) = XE(5,N)
      XVX(2) = XE(6,N)
      XVX(3) = XE(8,N)
      XVX(4) = XE(7,N)
      YVX(1) = YE(5,N)
      YVX(2) = YE(6,N)
      YVX(3) = YE(8,N)
      YVX(4) = YE(7,N)
      ZVX(1) = ZE(5,N)
      ZVX(2) = ZE(6,N)
      ZVX(3) = ZE(8,N)
      ZVX(4) = ZE(7,N)
      NP = 4
      CALL PGCNTRD( NP,XVX,YVX,ZVX,XFT,YFT,ZFT )
!
!---  Load rotation matrix  ---
!
      DXFX = ABS( XFE-XFW )
      IF( DXFX.LT.EPSL ) DXFX = 0.D+0
      DYFX = ABS( YFE-YFW )
      IF( DYFX.LT.EPSL ) DYFX = 0.D+0
      DZFX = ABS( ZFE-ZFW )
      IF( DZFX.LT.EPSL ) DZFX = 0.D+0
      DLX = SQRT( DXFX**2 + DYFX**2 + DZFX**2 )
      ROT_MAT(1,1) = DXFX/(DLX+SMALL)
      ROT_MAT(2,1) = DYFX/(DLX+SMALL)
      ROT_MAT(3,1) = DZFX/(DLX+SMALL)
      DXFY = ABS( XFN-XFS )
      IF( DXFY.LT.EPSL ) DXFY = 0.D+0
      DYFY = ABS( YFN-YFS )
      IF( DYFY.LT.EPSL ) DYFY = 0.D+0
      DZFY = ABS( ZFN-ZFS )
      IF( DZFY.LT.EPSL ) DZFY = 0.D+0
      DLY = SQRT( DXFY**2 + DYFY**2 + DZFY**2 )
      ROT_MAT(1,2) = DXFY/(DLY+SMALL)
      ROT_MAT(2,2) = DYFY/(DLY+SMALL)
      ROT_MAT(3,2) = DZFY/(DLY+SMALL)
      DXFZ = ABS( XFT-XFB )
      IF( DXFZ.LT.EPSL ) DXFZ = 0.D+0
      DYFZ = ABS( YFT-YFB )
      IF( DYFZ.LT.EPSL ) DYFZ = 0.D+0
      DZFZ = ABS( ZFT-ZFB )
      IF( DZFZ.LT.EPSL ) DZFZ = 0.D+0
      DLZ = SQRT( DXFZ**2 + DYFZ**2 + DZFZ**2 )
      ROT_MAT(1,3) = DXFZ/(DLZ+SMALL)
      ROT_MAT(2,3) = DYFZ/(DLZ+SMALL)
      ROT_MAT(3,3) = DZFZ/(DLZ+SMALL)
!
!---  Well projections in global coordinate system  ---
!
      PLGCX(1) = ABS(XP_WX(2)-XP_WX(1))
      PLGCX(2) = ABS(YP_WX(2)-YP_WX(1))
      PLGCX(3) = ABS(ZP_WX(2)-ZP_WX(1))
!
!---  Well length in global coordinate system  ---
!
      WLGCX = SQRT( PLGCX(1)**2 + PLGCX(2)**2 + PLGCX(3)**2 )
!
!---  Well projections in local coordinate system  ---
!
      DO 20 I = 1,3
        PLLCX(I) = 0.D+0
        DO 10 J = 1,3
          PLLCX(I) = PLLCX(I) + ROT_MAT(I,J)*PLGCX(J)
   10   CONTINUE
   20 CONTINUE
!
!---  Well length in local coordinate system  ---
!
      WLLCX = SQRT( PLLCX(1)**2 + PLLCX(2)**2 + PLLCX(3)**2 )
!
!---  Well projections in local coordinate system, corrected for 
!     global well length  ---
!
      PLX_WX = PLLCX(1)*WLGCX/WLLCX
      PLY_WX = PLLCX(2)*WLGCX/WLLCX
      PLZ_WX = PLLCX(3)*WLGCX/WLLCX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROJ_BH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDGEOM_FRC( JNDX )
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
!     Read input file for fractures comprised of planar triangles
!
!     XE_FRC(NV,NT) x-vertex dimension (m) for vertex NV, 
!       and triangle NT
!     YE_FRC(NV,NT) y-vertex dimension (m) for vertex NV, 
!       and triangle NT
!     ZE_FRC(NV,NT) z-vertex dimension (m) for vertex NV, 
!       and triangle NT
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 31 July 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE OUTPU
      USE GRID
      USE GEOM_FRC
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
      CHARACTER*64 ADUM,BDUM,FDUM,UNTS
      CHARACTER*132 CTFILE
      CHARACTER*512 CHDUM,CHDUM0
      CHARACTER*9 FORM1
      CHARACTER*21 FORM2
      CHARACTER*13 FORM4
      CHARACTER*18 FORM5
      CHARACTER*6 SURFX(6)
      REAL*8 XPX(6),YPX(6),ZPX(6)
      REAL*8, DIMENSION(2) :: GX
!      REAL*8, DIMENSION(1) :: SIGNX,VAMX,VDMX
      REAL*8, DIMENSION(1) :: VAMX,VDMX
      REAL*8, DIMENSION(3) :: VNX,VPX,VQX
      REAL*8, DIMENSION(3) :: VAX,VBX,VCX,VDX
      REAL*8, DIMENSION(3) :: P1X,P2X,P3X,PN1X,PN2X,PX
      REAL*8, DIMENSION(3) :: VX0X,VX1X,VX2X
      REAL*8, DIMENSION(3,3) :: SIGX
      REAL*8, DIMENSION(3,3) :: V1X,V2X
      REAL*8, DIMENSION(3) :: VT1X,VT2X
      REAL*8, DIMENSION(3,24) :: PTX
      REAL*8, DIMENSION(3,6) :: CPX
      REAL*8, DIMENSION(24) :: XPTX,YPTX,ZPTX,ROTX
      REAL*8 :: JRCX
      REAL*8, DIMENSION(:), ALLOCATABLE :: XCTX,YCTX,ZCTX
!      REAL*8, DIMENSION(3,3) :: RM1,RM2,RM3
      LOGICAL FLG_EX,FCHK
      INTEGER IVPX(3),MSX(4,6),N1X(4),N2X(4)
      INTEGER, DIMENSION(3) :: NODX
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(8X,A,I9)'/
      DATA FORM2 /'(8X,A,I9,A,I9,A,I9,A)'/
      DATA FORM4 /'(A,I3,A,I3,A)'/
      DATA FORM5 /'(A,I3,A,I3,A,I3,A)'/
      DATA MSX / 1,2,4,3,1,5,6,2,1,3,7,5,2,6,8,4,3,4,8,7,5,7,8,6 /
      DATA N1X / 2,3,4,1 /
      DATA N2X / 1,2,3,4 /
      DATA SURFX / 'Bottom','South ','West  ','East  ',
     &  'North ','Top   ' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDGEOM_FRC'
!
!---  Create a stress tensor  ---
!
      SIGX(1,1) = 25.5D+06
      SIGX(1,2) = 0.D+0
      SIGX(1,3) = 0.D+0
      SIGX(2,1) = 0.D+0
      SIGX(2,2) = 17.7D+06
      SIGX(2,3) = 0.D+0
      SIGX(3,1) = 0.D+0
      SIGX(3,2) = 0.D+0
      SIGX(3,3) = 31.8D+06
!
!---  Check for pre-processed fracture data  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
      IF( INDEX(ADUM(1:),'preprocess').NE.0 .AND.
     &  INDEX(ADUM(1:),'fracture').NE.0 ) THEN
        INQUIRE( FILE='ppfrac.dat',EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          INDX = 4     
          CHMSG = 'Missing Preprocessed Fracture/Fault File: ppfrac.dat'
          CALL WRMSGS( INDX )
        ELSE
          OPEN(UNIT=35, FILE='ppfrac.dat', 
     &      STATUS='OLD', FORM='FORMATTED')
        ENDIF
!
!---    Read fracture parameters  ---
!
        READ(35,*) LF_FRCX,LFC_FRCX,LNC_FRCX,LT_FRCX,LTC_FRCX,LXP_FRCX
!
!---    Read number of fractures and number of inactive fracture
!       triangles  ---
!
        READ(35,*) NF_FRC,NXP_FRC
!
!---    Loop over fractures/faults  ---
!
        DO NFX = 1,NF_FRC
!
!---      Read fracture/fault preprocessed file  ---
!
          CALL RDPREP_FRC( NFX )
        ENDDO     
!
!---    Close fracture/fault preprocessed file  ---
!
        CLOSE(UNIT=35)
!
!---    Reset subroutine string sequence  ---
!
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ELSE
        BACKSPACE(UNIT=IRD)
      ENDIF
!!
!!---  Create a stress tensor in local x,y,z coordinates  ---
!!
!      SIGX(1,1) = 35.5D+06
!      SIGX(1,2) = 0.D+0
!      SIGX(1,3) = 0.D+0
!      SIGX(2,1) = 0.D+0
!      SIGX(2,2) = 21.7D+06
!      SIGX(2,3) = 0.D+0
!      SIGX(3,1) = 0.D+0
!      SIGX(3,2) = 0.D+0
!      SIGX(3,3) = 41.8D+06
!!
!!---  Rotation matrix, starting with 4˚ yaw rotation about 
!!     Homestake z axis  ---
!!
!      THETA = 4.D+0*GPI/180.D+0
!      RM1(1,1) = COS(THETA)
!      RM1(1,2) = -SIN(THETA)
!      RM1(1,3) = 0.D+0
!      RM1(2,1) = SIN(THETA)
!      RM1(2,2) = COS(THETA)
!      RM1(2,3) = 0.D+0
!      RM1(3,1) = 0.D+0
!      RM1(3,2) = 0.D+0
!      RM1(3,3) = 1.D+0
!!
!!---  Then 12˚ roll rotation about the x' axis  ---
!!
!      THETA = -12.D+0*GPI/180.D+0
!      RM2(1,1) = 1.D+0
!      RM2(1,2) = 0.D+0
!      RM2(1,3) = 0.D+0
!      RM2(2,1) = 0.D+0
!      RM2(2,2) = COS(THETA)
!      RM2(2,3) = -SIN(THETA)
!      RM2(3,1) = 0.D+0
!      RM2(3,2) = SIN(THETA)
!      RM2(3,3) = COS(THETA)
!!
!!---  General rotation  ---
!!
!      CALL MATMUL( RM1,RM2,RM3,3,3,3 )
!!
!!---  Transpose rotation matrix to create a rotation matrix from
!!     local x,y,z to global (i.e., Homestake) x,y,z coordinates
!!
!      DO N = 1,3
!        DO M = 1,3
!          RM1(N,M) = RM3(M,N)
!        ENDDO
!      ENDDO
!
!---  Write card information to ouput file  ---
!
      IF( JNDX.EQ.1 ) THEN
        CARD = 'Fracture Geometry Card'
      ELSEIF( JNDX.EQ.2 ) THEN
        CARD = 'Fault Geometry Card'
      ELSE
        INDX = 4
        CARD = 'Fracture/Fault Geometry Card'
        CHMSG = 'Unrecognized Card Option'
        CALL WRMSGS( INDX )
      ENDIF
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!      OPEN(UNIT=16, FILE='stomp_tnc.dat', STATUS='UNKNOWN', 
!     &  FORM='FORMATTED')
!      CLOSE(UNIT=16,STATUS='DELETE')
!      OPEN(UNIT=16, FILE='stomp_tnc.dat', STATUS='NEW', 
!     &  FORM='FORMATTED')
!
!---  Write fracture data to 'fracture'  ---
!
      IF( JNDX.EQ.1 ) THEN
        FDUM = 'fractures'
      ELSEIF( JNDX.EQ.2 ) THEN
        FDUM = 'faults'
      ENDIF
      OPEN(UNIT=36, FILE=TRIM(ADJUSTL(FDUM)), STATUS='UNKNOWN', 
     &  FORM='FORMATTED')
      CLOSE(UNIT=36,STATUS='DELETE')
      OPEN(UNIT=36, FILE=TRIM(ADJUSTL(FDUM)), STATUS='NEW', 
     &  FORM='FORMATTED')
!
!---  Read the number of fractures  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      IF( JNDX.EQ.1 ) THEN
        VARB = 'Number of Fractures'
      ELSEIF( JNDX.EQ.2 ) THEN
        VARB = 'Number of Faults'
      ENDIF
      CALL RDINT(ISTART,ICOMMA,CHDUM,NF_FRC)
      IF( JNDX.EQ.1 ) THEN
        WRITE(IWR,'(2X,A,I9)') 'Number of Fractures:   ',NF_FRC
        WRITE(36,'(2X,A,I9)') 'Number of Fractures:   ',NF_FRC
      ELSEIF( JNDX.EQ.2 ) THEN
        WRITE(IWR,'(2X,A,I9)') 'Number of Faults:   ',NF_FRC
        WRITE(36,'(2X,A,I9)') 'Number of Faults:   ',NF_FRC
      ENDIF
!
!---  Loop over the number of fractures  ---
!
      NTX = 0
      L1 : DO NFX = 1,NF_FRC
        IP_FRC(1,NFX) = NTX + 1
        IF( JNDX.EQ.1 ) THEN
          WRITE(IWR,'(/,2X,A,I9)') 'Fracture #:   ',NFX
          WRITE(36,'(/,2X,A,I9)') 'Fracture #:   ',NFX
        ELSEIF( JNDX.EQ.2 ) THEN
          WRITE(IWR,'(/,2X,A,I9)') 'Faults #:   ',NFX
          WRITE(36,'(/,2X,A,I9)') 'Faults #:   ',NFX
        ENDIF     
!
!---    Read new input line  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CHDUM0 = CHDUM
        CALL LCASE( CHDUM )
!
!---    Check for an external file  ---
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM0,FDUM)
          INQUIRE( FILE=TRIM(ADJUSTL(FDUM)),EXIST=FCHK )
          IF( .NOT.FCHK ) THEN
            INDX = 4     
            IF( JNDX.EQ.1 ) THEN
              CHMSG = 'Missing Fracture File: ' // TRIM(ADJUSTL(FDUM))
            ELSEIF( JNDX.EQ.2 ) THEN
              CHMSG = 'Missing Fault File: ' //  TRIM(ADJUSTL(FDUM))
            ENDIF
            CALL WRMSGS( INDX )
          ELSE
            OPEN(UNIT=35, FILE=TRIM(ADJUSTL(FDUM)), 
     &        STATUS='OLD', FORM='FORMATTED')
            IRDX = IRD
            IRD = 35
            IFILE = 1
            IF( JNDX.EQ.1 ) THEN
              ISTART = 1
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              WRITE(IWR,'(4X,2A)') 'Fracture Data Filename: ',
     &          TRIM(ADJUSTL(FDUM))
              WRITE(36,'(4X,2A)') 'Fracture Data Filename: ',
     &          TRIM(ADJUSTL(FDUM))
            ELSEIF( JNDX.EQ.2 ) THEN
              ISTART = 1
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              WRITE(IWR,'(4X,2A)') 'Fault Data Filename: ',
     &          TRIM(ADJUSTL(FDUM))
              WRITE(36,'(4X,2A)') 'Fault Data Filename: ',
     &          TRIM(ADJUSTL(FDUM))
            ENDIF
          ENDIF
        ELSE
          ISTART = 1
          IFILE = 0
        ENDIF
!
!---    Read the number of triangles in the fracture  ---
!
        IF( JNDX.EQ.1 ) THEN
          VARB = 'Number of Triangles in Fracture'
        ELSEIF( JNDX.EQ.2 ) THEN
          VARB = 'Number of Triangles in Faults'
        ENDIF     
        CALL RDINT(ISTART,ICOMMA,CHDUM,NTP_FRC(NFX))    
        IF( JNDX.EQ.1 ) THEN
          WRITE(IWR,'(4X,A,I9)') 'Number of Triangles in Fracture:   ',
     &      NTP_FRC(NFX)
          WRITE(36,'(4X,A,I9)') 'Number of Triangles in Fracture:   ',
     &      NTP_FRC(NFX)
        ELSEIF( JNDX.EQ.2 ) THEN
          WRITE(IWR,'(4X,A,I9)') 'Number of Triangles in Fault:   ',
     &      NTP_FRC(NFX)
          WRITE(36,'(4X,A,I9)') 'Number of Triangles in Fault:   ',
     &      NTP_FRC(NFX)
        ENDIF     
        NT_FRC = NT_FRC + NTP_FRC(NFX)
        IF( NT_FRC.GT.LT_FRC ) THEN
          INDX = 5
          IF( JNDX.EQ.1 ) THEN
            CHMSG = 'Number of Fracture Triangles > Parameter LT_FRC'
          ELSEIF( JNDX.EQ.2 ) THEN
            CHMSG = 'Number of Fault Triangles > Parameter LT_FRC'
          ENDIF     
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read fracture joint model  ---
!
        IF( JNDX.EQ.1 ) THEN
          VARB = 'Fracture Joint Model'
        ELSEIF( JNDX.EQ.2 ) THEN
          VARB = 'Fault Joint Model'
        ENDIF     
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX(ADUM(1:),'dynamic').NE.0 .AND.
     &    INDEX(ADUM(1:),'sneddon').NE.0 .AND.
     &    INDEX(ADUM(1:),'barton').NE.0 .AND.
     &    INDEX(ADUM(1:),'bandis').NE.0 ) THEN
          WRITE(IWR,'(4X,A)') 'Dynamic Aperture Model'
          WRITE(IWR,'(4X,A)') 'w/ Sneddon Open Model'
          WRITE(IWR,'(4X,A)') 'w/ Barton-Bandis Joint Model'
          WRITE(36,'(4X,A)') 'Dynamic Aperture Model'
          WRITE(36,'(4X,A)') 'w/ Sneddon Open Model'
          WRITE(36,'(4X,A)') 'w/ Barton-Bandis Joint Model'
          IJM_FRC(NFX) = 11
          VARB = 'Aperture at Zero Net Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DSBB_FRC(1,NFX))
          VARB = 'Aperture at Zero Net Pressure Units'
          UNTS = 'm'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IUNM = 1
          IDFLT = 1
          WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',DSBB_FRC(1,NFX)
          WRITE(36,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',DSBB_FRC(1,NFX)
          INDX = 0
          CALL RDUNIT(UNTS,DSBB_FRC(1,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',DSBB_FRC(1,NFX),', m)'
          WRITE(36,'(A,1PE11.4,A)') ' (',DSBB_FRC(1,NFX),', m)'
!          VARB = 'Young''s Modulus'
!          CALL RDDPR(ISTART,ICOMMA,CHDUM,DSBB_FRC(2,NFX))
!          VARB = 'Young''s Modulus Units'
!          UNTS = 'pa'
!          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
!          IUNM = -1
!          IUNKG = 1
!          IUNS = -2
!          IDFLT = 1
!          WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
!     &      UNTS(1:NCH),': ',DSBB_FRC(2,NFX)
!          INDX = 0
!          CALL RDUNIT(UNTS,DSBB_FRC(2,NFX),INDX)
!          WRITE(IWR,'(A,1PE11.4,A)') ' (',DSBB_FRC(2,NFX),', Pa)'
!          VARB = 'Poisson''s Ratio'
!          CALL RDDPR(ISTART,ICOMMA,CHDUM,DSBB_FRC(3,NFX))
!          WRITE(IWR,'(6X,2A,1PE11.4)') VARB(1:IVR),': ',
!     &      DSBB_FRC(3,NFX)
          VARB = 'Fracture Toughness'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DSBB_FRC(2,NFX))
          VARB = 'Fracture Toughness Units'
          UNTS = 'pa/m'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IUNM = -2
          IUNKG = 1
          IUNS = -2
          IDFLT = 1
          WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',DSBB_FRC(2,NFX)
          WRITE(36,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',DSBB_FRC(2,NFX)
          INDX = 0
          CALL RDUNIT(UNTS,DSBB_FRC(2,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',DSBB_FRC(2,NFX),', Pa/m)'
          WRITE(36,'(A,1PE11.4,A)') ' (',DSBB_FRC(2,NFX),', Pa/m)'
        ELSEIF( INDEX(ADUM(1:),'barton').NE.0 .AND.
     &    INDEX(ADUM(1:),'bandis').NE.0 ) THEN
          WRITE(IWR,'(4X,A)') 'Barton-Bandis Joint Model'
          WRITE(36,'(4X,A)') 'Barton-Bandis Joint Model'
          VARB = 'Barton-Bandis Joint Model Option'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'roughness').NE.0 .AND.
     &      INDEX(BDUM(1:),'coeff').NE.0 ) THEN
            WRITE(IWR,'(4X,A)') 'Roughness Coefficient Input Option'
            WRITE(36,'(4X,A)') 'Roughness Coefficient Input Option'
            IJM_FRC(NFX) = 1
          ELSEIF( INDEX(BDUM(1:),'unstressed').NE.0 .AND.
     &      INDEX(BDUM(1:),'aperture').NE.0 ) THEN
            WRITE(IWR,'(4X,A)') 'Unstressed Aperture Input Option'
            WRITE(36,'(4X,A)') 'Unstressed Aperture Input Option'
            IJM_FRC(NFX) = 2
          ELSEIF( INDEX(BDUM(1:),'residual').NE.0 .AND.
     &      INDEX(BDUM(1:),'aperture').NE.0 ) THEN
            WRITE(IWR,'(4X,A)') 'Residual Aperture Input Option'
            WRITE(36,'(4X,A)') 'Residual Aperture Input Option'
            IJM_FRC(NFX) = 3
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Barton-Bandis Joint Model ' // 
     &        'Option: ' // BDUM
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'constant').NE.0 .AND.
     &    INDEX(ADUM(1:),'aperture').NE.0 ) THEN
          WRITE(IWR,'(4X,A)') 'Constant Aperture Model'
          WRITE(36,'(4X,A)') 'Constant Aperture Model'
          IJM_FRC(NFX) = 0
        ELSEIF( INDEX(ADUM(1:),'sneddon').NE.0 .AND.
     &    INDEX(ADUM(1:),'aperture').NE.0 ) THEN
          WRITE(IWR,'(4X,A)') 'Sneddon Aperture Distribution'
          WRITE(36,'(4X,A)') 'Sneddon Aperture Distribution'
          IJM_FRC(NFX) = 10
          VARB = 'Flow Channel Radius'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,FCRX)
          VARB = 'Flow Channel Radius Units'
          UNTS = 'm'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IUNM = 1
          IDFLT = 1
          WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',FCRX
          WRITE(36,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',FCRX
          INDX = 0
          CALL RDUNIT(UNTS,FCRX,INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',FCRX,', m)'
          WRITE(36,'(A,1PE11.4,A)') ' (',FCRX,', m)'
          VARB = 'Maximum Aperture'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,BMAX)
          VARB = 'Maximum Aperture Units'
          UNTS = 'm'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IUNM = 1
          IDFLT = 1
          WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',BMAX
          WRITE(36,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',BMAX
          INDX = 0
          CALL RDUNIT(UNTS,BMAX,INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',BMAX,', m)'
          WRITE(36,'(A,1PE11.4,A)') ' (',BMAX,', m)'
          VARB = 'Channel Trajectory File'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,CTFILE)
          CTFILE = TRIM(ADJUSTL(CTFILE))
          INQUIRE( FILE=CTFILE, EXIST=FLG_EX )
          IF( .NOT.FLG_EX ) THEN
            NCH = INDEX( CTFILE(1:),'  ' )-1
            CHMSG = 'Missing Channel Trajectory File: ' // CTFILE(1:NCH)
            INDX = 4
            CALL WRMSGS( INDX )
          ELSE
            OPEN(UNIT=25, FILE=CTFILE, STATUS='OLD', FORM='FORMATTED')
          ENDIF
!  
!---      Read number channel trajectory points  ---
!
          READ(25,*) NCTPX
!  
!         Allocate memory for channel trajectory points
!
          ALLOCATE( XCTX(1:NCTPX),STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            CHMSG = 'Allocation Error: XCTX'
            INDX = 3
            CALL WRMSGS( INDX )
          ENDIF
          ALLOCATE( YCTX(1:NCTPX),STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            CHMSG = 'Allocation Error: YCTX'
            INDX = 3
            CALL WRMSGS( INDX )
          ENDIF 
          ALLOCATE( ZCTX(1:NCTPX),STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            CHMSG = 'Allocation Error: ZCTX'
            INDX = 3
            CALL WRMSGS( INDX )
          ENDIF  
          VARB = 'Channel Trajectory File Units'
          UNTS = 'm'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IUNM = 1
          IDFLT = 1
          WRITE(IWR,'(6X,3A)') VARB(1:IVR),', ',UNTS(1:NCH)
          WRITE(36,'(6X,3A)') VARB(1:IVR),', ',UNTS(1:NCH)
          INDX = 0
          VCTX = 1.D+0
          CALL RDUNIT(UNTS,VCTX,INDX)
!  
!---      Read x,y,z channel trajectory points and convert to 'm'
!
          DO N = 1,NCTPX
            READ(25,*) XCTX(N),YCTX(N),ZCTX(N)
            XCTX(N) = VCTX*XCTX(N)
            YCTX(N) = VCTX*YCTX(N)
            ZCTX(N) = VCTX*ZCTX(N)
          ENDDO
        ELSE
          INDX = 4
          IF( JNDX.EQ.1 ) THEN
            CHMSG = 'Unrecognized Fracture Joint Model: '//ADUM
          ELSEIF( JNDX.EQ.2 ) THEN
            CHMSG = 'Unrecognized Fault Joint Model: '//ADUM
          ENDIF     
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read the fracture name  ---
!
        IF( JNDX.EQ.1 ) THEN
          VARB = 'Fracture Name'
        ELSEIF( JNDX.EQ.2 ) THEN
          VARB = 'Fault Name'
        ENDIF     
        IDFLT = 1
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,NM_FRC(NFX))
        IF( JNDX.EQ.1 ) THEN
          WRITE (IWR,'(4X,2A)') 'Fracture Name: ',NM_FRC(NFX)
          WRITE (36,'(4X,2A)') 'Fracture Name: ',NM_FRC(NFX)
        ELSEIF( JNDX.EQ.2 ) THEN
          WRITE (IWR,'(4X,2A)') 'Fault Name: ',NM_FRC(NFX)
          WRITE (36,'(4X,2A)') 'Fault Name: ',NM_FRC(NFX)
        ENDIF     
!
!---    Read the fracture-fracture intersection skin factor  ---
!
        IF( JNDX.EQ.1 ) THEN
          VARB = 'Fracture-Fracture Intersection Skin Factor'
        ELSEIF( JNDX.EQ.2 ) THEN
          VARB = 'Fault-Fault Intersection Skin Factor'
        ENDIF     
        CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
        IF( INDX.EQ.1 ) THEN
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SKF_FRC(NFX))
          IF( JNDX.EQ.1 ) THEN
            WRITE (IWR,'(4X,A,1PE11.4)') 'Fracture-Fracture '
     &        // 'Intersection Skin Factor: ',SKF_FRC(NFX)
            WRITE (36,'(4X,A,1PE11.4)') 'Fracture-Fracture '
     &        // 'Intersection Skin Factor: ',SKF_FRC(NFX)
          ELSEIF( JNDX.EQ.2 ) THEN
            WRITE (IWR,'(4X,A,1PE11.4)') 'Fault-Fault Intersection '
     &        // 'Skin Factor: ',SKF_FRC(NFX)
            WRITE (36,'(4X,A,1PE11.4)') 'Fault-Fault Intersection '
     &        // 'Skin Factor: ',SKF_FRC(NFX)
          ENDIF     
        ENDIF
!
!---    Loop over the number of triangles in the fracture  ---
!
        DO MTX = 1,NTP_FRC(NFX)
          NTX = MTX + IP_FRC(1,NFX) - 1
          WRITE (36,'(/,4X,A,I9)') 'Triangle #:   ',NTX
!
!---      Read parameters and vertices for each triangle  ---
!
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
!
!---      Read triangle aperture  ---
!
          IF( IJM_FRC(NFX).EQ.0 ) THEN
            VARB = 'Triangle Aperture'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,JRC_FRC(NTX))
            VARB = 'Triangle Aperture Units'
            UNTS = 'm'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IUNM = 1
            IDFLT = 1
            WRITE(36,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',JRC_FRC(NTX)
            INDX = 0
            CALL RDUNIT(UNTS,JRC_FRC(NTX),INDX)
            WRITE(36,'(A,1PE11.4,A)') ' (',JRC_FRC(NTX),', m)'
!
!---      Read triangle joint roughness coefficient  ---
!
          ELSEIF( IJM_FRC(NFX).EQ.1 ) THEN
            VARB = 'Triangle Joint Roughness Coefficient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,JRC_FRC(NTX))
            IDFLT = 1
            WRITE(36,'(6X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        JRC_FRC(NTX)
!
!---      Read unstressed aperture, m  ---
!
          ELSEIF( IJM_FRC(NFX).EQ.2 ) THEN
            VARB = 'Triangle Joint Unstressed Aperture'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VARX)
            VARB = 'Triangle Joint Unstressed Aperture Units'
            UNTS = 'm'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IUNM = 1
            IDFLT = 1
            WRITE(36,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VARX
              INDX = 0
              CALL RDUNIT(UNTS,VARX,INDX)
            WRITE(36,'(A,1PE11.4,A)') ' (',VARX,', m)'
!
!---      Read residual aperture, m  ---
!
          ELSEIF( IJM_FRC(NFX).EQ.3 ) THEN
            VARB = 'Triangle Joint Residual Aperture'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VARX)
            VARB = 'Triangle Joint Residual Aperture Units'
            UNTS = 'm'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IUNM = 1
            IDFLT = 1
            WRITE(36,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VARX
              INDX = 0
              CALL RDUNIT(UNTS,VARX,INDX)
            WRITE(36,'(A,1PE11.4,A)') ' (',VARX,', m)'
          ENDIF
!
!---      Read triangle joint wall compressive strength  ---
!
          IF( IJM_FRC(NFX).GE.1 .AND. IJM_FRC(NFX).LE.3 ) THEN
            VARB = 'Triangle Joint Wall Compressive Strength'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,JCS_FRC(NTX))
            VARB = 'Triangle Joint Wall Compressive Strength Units'
              UNTS = 'pa'
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              IDFLT = 1
            WRITE(36,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',JCS_FRC(NTX)
            INDX = 0
            CALL RDUNIT(UNTS,JCS_FRC(NTX),INDX)
            WRITE(36,'(A,1PE11.4,A)') ' (',JCS_FRC(NTX),', Pa)'
!
!---        Read triangle unconfined compressive strength of rock
!             adjacent to joint wall  ---
!
            VARB = 'Triangle Unconfined Compressive Strength ' // 
     &          '(Adjac. Rock)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,UCS_FRC(NTX))
            VARB = 'Triangle Unconfined Compressive Strength Units'
              UNTS = 'pa'
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              IDFLT = 1
            WRITE(36,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',UCS_FRC(NTX)
            INDX = 0
            CALL RDUNIT(UNTS,UCS_FRC(NTX),INDX)
            WRITE(36,'(A,1PE11.4,A)') ' (',UCS_FRC(NTX),', Pa)'
          ENDIF
!
!---      Barton-Bandis Joint Model, with specified joint 
!         roughness coefficient  ---
!
          IF( IJM_FRC(NFX).EQ.1 ) THEN
!
!---        Initial joint aperture (mm), eqn. (5), 
!           Rock Mech Rock Eng. (2016) 49:837-853  ---
!
            AJX = 2.D-1*JRC_FRC(NTX)*(2.D-1*UCS_FRC(NTX)/
     &        JCS_FRC(NTX) - 1.D-1)
            IF( AJX.LT.0.D+0 ) THEN
              INDX = 7
              IMSG = NTX
              IF( JNDX.EQ.1 ) THEN
                CHMSG = 'Barton-Bandis Joint Model: Negative ' //
     &            'Initial Joint Aperture @ Fracture Triangle'
              ELSEIF( JNDX.EQ.2 ) THEN
                CHMSG = 'Barton-Bandis Joint Model: Negative ' //
     &            'Initial Joint Aperture @ Fault Triangle'
              ENDIF     
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Maximum joint closure (mm), eqn. (4), 
!           Rock Mech Rock Eng.(2016) 49:837-853, using joint 
!           compressive strength in untis of MPa  ---
!
            CSJX = 1.D-6*JCS_FRC(NTX)
            VMX = -2.96D-1 - 5.6D-3*JRC_FRC(NTX) + 
     &          2.241D+0*((CSJX/AJX)**(-2.45D-1))
            WRITE(36,'(6X,A,1PE11.4)') 'Unstressed Joint ' // 
     &        'Aperture, mm :',AJX
            WRITE(36,'(6X,A,1PE11.4)') 'Residual Joint ' // 
     &        'Aperture, mm :',MAX( AJX-VMX,0.D+0 )
!
!---      Barton-Bandis Joint Model, with specified unstressed
!         aperture  ---
!
          ELSEIF( IJM_FRC(NFX).EQ.2 ) THEN
            VARX = 1.D+3*VARX
!
!---        Initial joint aperture (mm), eqn. (5), 
!           Rock Mech Rock Eng. (2016) 49:837-853  ---
!
            JRC_FRC(NTX) = 5.D+0*VARX/(2.D-1*UCS_FRC(NTX)/
     &        JCS_FRC(NTX) - 1.D-1)
!
!---        Initial joint aperture (mm), eqn. (5), 
!           Rock Mech Rock Eng. (2016) 49:837-853  ---
!
            AJX = 2.D-1*JRC_FRC(NTX)*(2.D-1*UCS_FRC(NTX)/
     &        JCS_FRC(NTX) - 1.D-1)
            IF( AJX.LT.0.D+0 ) THEN
              INDX = 7
              IMSG = NTX
              IF( JNDX.EQ.1 ) THEN
                CHMSG = 'Barton-Bandis Joint Model: Negative ' //
     &            'Initial Joint Aperture @ Fracture Triangle'
              ELSEIF( JNDX.EQ.2 ) THEN
                CHMSG = 'Barton-Bandis Joint Model: Negative ' //
     &            'Initial Joint Aperture @ Fault Triangle'
              ENDIF     
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Maximum joint closure (mm), eqn. (4), 
!           Rock Mech Rock Eng.(2016) 49:837-853, using joint 
!           compressive strength in untis of MPa  ---
!
            CSJX = 1.D-6*JCS_FRC(NTX)
            VMX = -2.96D-1 - 5.6D-3*JRC_FRC(NTX) + 
     &          2.241D+0*((CSJX/AJX)**(-2.45D-1))
            WRITE(36,'(6X,A,1PE11.4)') 'Joint Roughness ' // 
     &       'Coefficient :',JRC_FRC(NTX)
            WRITE(36,'(6X,A,1PE11.4)') 'Residual Joint ' // 
     &         'Aperture, mm :',MAX( AJX-VMX,0.D+0 )
!
!---      Barton-Bandis Joint Model, with specified residual
!         aperture  ---
!
          ELSEIF( IJM_FRC(NFX).EQ.3 ) THEN
            VARX = 1.D+3*VARX
            NC = 0
            JRCX = 1.D+0
            DO
              NC = NC + 1
              DO M = 1,2
                IF( M.EQ.2 ) JRCX = JRCX + 1.D-6
                AJX = 2.D-1*JRCX*(2.D-1*UCS_FRC(NTX)/
     &          JCS_FRC(NTX) - 1.D-1)
                CSJX = 1.D-6*JCS_FRC(NTX)
                VMX = -2.96D-1 - 5.6D-3*JRCX + 
     &            2.241D+0*((CSJX/AJX)**(-2.45D-1))
                GX(M) = VARX - MAX( AJX-VMX,0.D+0 )
              ENDDO
              FX = GX(1)
              IF( ABS(FX).LT.EPSL ) EXIT
              DFX = (GX(2)-GX(1))/1.D-6
              DJRCX = -FX/DFX
              JRCX = JRCX + DJRCX
              IF( ABS(DJRCX).LT.1.D-6 ) EXIT
            ENDDO
            JRC_FRC(NTX) = JRCX
            WRITE(36,'(6X,A,1PE11.4)') 'Joint Roughness ' // 
     &       'Coefficient :',JRC_FRC(NTX)
!
!---        Initial joint aperture (mm), eqn. (5), 
!           Rock Mech Rock Eng. (2016) 49:837-853  ---
!
            AJX = 2.D-1*JRC_FRC(NTX)*(2.D-1*UCS_FRC(NTX)/
     &        JCS_FRC(NTX) - 1.D-1)
            IF( AJX.LT.0.D+0 ) THEN
              INDX = 7
              IMSG = NTX
              IF( JNDX.EQ.1 ) THEN
                CHMSG = 'Barton-Bandis Joint Model: Negative ' //
     &            'Initial Joint Aperture @ Fracture Triangle'
              ELSEIF( JNDX.EQ.2 ) THEN
                CHMSG = 'Barton-Bandis Joint Model: Negative ' //
     &            'Initial Joint Aperture @ Fault Triangle'
              ENDIF     
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Maximum joint closure (mm), eqn. (4), 
!           Rock Mech Rock Eng.(2016) 49:837-853, using joint 
!           compressive strength in untis of MPa  ---
!
            CSJX = 1.D-6*JCS_FRC(NTX)
            VMX = -2.96D-1 - 5.6D-3*JRC_FRC(NTX) + 
     &        2.241D+0*((CSJX/AJX)**(-2.45D-1))
            WRITE(36,'(6X,A,1PE11.4)') 'Unstressed Joint ' // 
     &         'Aperture, mm :',AJX
          ENDIF
!
!---      Loop over the three vertices in the triangle  ---
!
          DO NVX = 1,3
            VARB = 'X Vertex Dimension'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,XE_FRC(NVX,NTX))
            VARB = 'X Vertex Dimension Units'
            UNTS = 'm'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IUNM = 1
            IDFLT = 1
            INDX = 0
            CALL RDUNIT(UNTS,XE_FRC(NVX,NTX),INDX)
            VARB = 'Y Vertex Dimension'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,YE_FRC(NVX,NTX))
            VARB = 'Y Vertex Dimension Units'
            UNTS = 'm'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IUNM = 1
            IDFLT = 1
            INDX = 0
            CALL RDUNIT(UNTS,YE_FRC(NVX,NTX),INDX)
            VARB = 'Z Vertex Dimension'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,ZE_FRC(NVX,NTX))
            VARB = 'Z Vertex Dimension Units'
            UNTS = 'm'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IUNM = 1
            IDFLT = 1
            INDX = 0
            CALL RDUNIT(UNTS,ZE_FRC(NVX,NTX),INDX)
            IF( IJM_FRC(NFX).NE.10 ) 
     &        WRITE(36,'(6X,A,I1,3(A,1PE11.4))') 'Vertex # ',NVX,
     &        '; X, m: ',XE_FRC(NVX,NTX),
     &        '; Y, m: ',YE_FRC(NVX,NTX),
     &        '; Z, m: ',ZE_FRC(NVX,NTX)
          ENDDO
!
!---      Triangle area  ---
!
          CALL TRGAREA( XE_FRC(1,NTX),YE_FRC(1,NTX),
     &      ZE_FRC(1,NTX),XE_FRC(2,NTX),YE_FRC(2,NTX),
     &      ZE_FRC(2,NTX),XE_FRC(3,NTX),YE_FRC(3,NTX),
     &      ZE_FRC(3,NTX),APX )
          WRITE(36,'(6X,A,1PE11.4)') 'Triangle Area, m: ',APX
!
!---      Triangle aperture determined from Sneddon distribution  ---
!
          IF( IJM_FRC(NFX).EQ.10 ) THEN
            XIX = 1.D+20
            XOX = 0.D+0
!  
!---        Centroid of the triangle
!
            VX0X(1) = 0.D+0
            VX0X(2) = 0.D+0
            VX0X(3) = 0.D+0
            DO NVX = 1,3
              VX0X(1) = VX0X(1) + XE_FRC(NVX,NTX)
              VX0X(2) = VX0X(2) + YE_FRC(NVX,NTX)
              VX0X(3) = VX0X(3) + ZE_FRC(NVX,NTX)
            ENDDO
            VX0X(1) = VX0X(1)/3.D+0
            VX0X(2) = VX0X(2)/3.D+0
            VX0X(3) = VX0X(3)/3.D+0
            IF( VX0X(1).LT.XIX ) THEN
              XIX = VX0X(1)
              NTI = NTX
            ENDIF
            IF( VX0X(1).GT.XOX ) THEN
              XOX = VX0X(1)
              NTO = NTX
            ENDIF    
!  
!---        Find the minimum distance between the triangle centroid and 
!           channel trajectory  ---
!
            DISTX = 1.D+3
            DO M = 1,NCTPX-1
              VX1X(1) = XCTX(M)
              VX1X(2) = YCTX(M)
              VX1X(3) = ZCTX(M)
              VX2X(1) = XCTX(M+1)
              VX2X(2) = YCTX(M+1)
              VX2X(3) = ZCTX(M+1)
              DO N = 1,3
                VT1X(N) = VX1X(N)-VX0X(N)
                VT2X(N) = VX2X(N)-VX1X(N)
              ENDDO
              VRMX = VDOTP( VT1X,VT2X )
              VRDX = VT2X(1)**2 + VT2X(2)**2 + VT2X(3)**2
              TX = -VRMX/VRDX
              TX = MAX( MIN( TX,1.D+0 ),0.D+0 )
              DX = 0.D+0
              DO N = 1,3
                DX = DX + (VX1X(N)-VX0X(N))**2 
     &            + 2.D+0*TX*((VX2X(N)-VX1X(N))*(VX1X(N)-VX0X(N)))
     &            + (TX**2)*((VX2X(N)-VX1X(N))**2)
              ENDDO
              DISTX = MIN( DISTX,SQRT(DX) )
            ENDDO
!
!---        Sneddon function for aperture, m  ---
!
            RR = MIN( DISTX/FCRX,1.D+0 )
            JRC_FRC(NTX) = MAX( 1.D+3*BMAX*SQRT(1.D+0-(RR**2)),1.D-6 )
            VARB = 'Triangle Aperture'
            WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        'mm',': ',JRC_FRC(NTX)
            WRITE(36,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        'mm',': ',JRC_FRC(NTX)
            JRC_FRC(NTX) = 1.D-3*JRC_FRC(NTX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',JRC_FRC(NTX),', m)'
            WRITE(36,'(A,1PE11.4,A)') ' (',JRC_FRC(NTX),', m)'
            WRITE(IWR,'(6X,A,I1,3(A,1PE11.4))') 'Vertex # ',NVX,
     &        '; X, m: ',XE_FRC(NVX,NTX),
     &        '; Y, m: ',YE_FRC(NVX,NTX),
     &        '; Z, m: ',ZE_FRC(NVX,NTX)
            WRITE(36,'(6X,A,I1,3(A,1PE11.4))') 'Vertex # ',NVX,
     &        '; X, m: ',XE_FRC(NVX,NTX),
     &        '; Y, m: ',YE_FRC(NVX,NTX),
     &        '; Z, m: ',ZE_FRC(NVX,NTX)
!
!---      Dynamic Sneddon-Barton-Bandis read Sneddon radius
!         and fraction  ---
!
          ELSEIF( IJM_FRC(NFX).EQ.11 ) THEN
            VARB = 'Sneddon Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,JRC_FRC(NTX))
            VARB = 'Sneddon Radius Units'
            UNTS = 'm'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IUNM = 1
            IDFLT = 1
            WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',JRC_FRC(NTX)
            WRITE(36,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',JRC_FRC(NTX)
            INDX = 0
            CALL RDUNIT(UNTS,JRC_FRC(NTX),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',JRC_FRC(NTX),', m)'
            WRITE(36,'(A,1PE11.4,A)') ' (',JRC_FRC(NTX),', m)'
            VARB = 'Sneddon Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,JCS_FRC(NTX))
            WRITE(IWR,'(6X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',JCS_FRC(NTX)
            WRITE(36,'(6X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',JCS_FRC(NTX)
          ENDIF
        ENDDO
!
!---    Reset fracture type to constant aperture  ---
!
        IF( IJM_FRC(NFX).EQ.10 ) IJM_FRC(NFX) = 0
!
!---    Check that triangles of a fracture are connected and not
!       duplicates  ---
!
        DO MT1X = 1,NTP_FRC(NFX)
          NT1X = MT1X + IP_FRC(1,NFX) - 1
          IFLG1 = 0
          DO MT2X = 1,NTP_FRC(NFX)
            NT2X = MT2X + IP_FRC(1,NFX) - 1
            IFLG2 = 0
            IF( NT1X.EQ.NT2X ) CYCLE
            DO NV1X = 1,3
              DO NV2X = 1,3
                DXE = ABS(XE_FRC(NV1X,NT1X)-XE_FRC(NV2X,NT2X))
                DYE = ABS(YE_FRC(NV1X,NT1X)-YE_FRC(NV2X,NT2X))
                DZE = ABS(ZE_FRC(NV1X,NT1X)-ZE_FRC(NV2X,NT2X))
                IF( (DXE+DYE+DZE)/EPSL.LT.EPSL ) IFLG2 = IFLG2 + 1
              ENDDO
            ENDDO
!
!---        Duplicate triangles  ---
!
            IF( IFLG2.EQ.3 .AND. IFLG2.EQ.3 ) THEN
              INDX = 4
              WRITE(FORM5(5:5),'(I1)') ICOUNT(NF_FRC)
              WRITE(FORM5(10:10),'(I1)') ICOUNT(NT1X)
              WRITE(FORM5(15:15),'(I1)') ICOUNT(NT2X)
              NC = ICOUNT(NF_FRC) + ICOUNT(NT1X) + ICOUNT(NT2X) + 18
              WRITE(BDUM(1:NC),FORM5) '(FRAC',NF_FRC,': TRIA',NT1X,
     &          ': TRIA',NT2X,')'
              CHMSG = 'Duplicate Triangle: ' // BDUM(1:NC)
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Connected triangles  ---
!
            IF( IFLG2.GE.2 .AND. NTP_FRC(NFX).GT.1 ) IFLG1 = 1
          ENDDO
!
!---      Unconnected triangle  ---
!
          IF( IFLG1.EQ.0 .AND. NTP_FRC(NFX).GT.1 ) THEN
            INDX = 4
            WRITE(FORM4(5:5),'(I1)') ICOUNT(NF_FRC)
            WRITE(FORM4(10:10),'(I1)') ICOUNT(NT1X)
            NC = ICOUNT(NF_FRC) + ICOUNT(NT1X) + 12
            WRITE(BDUM(1:NC),FORM4) '(FRAC',NF_FRC,': TRIA',NT1X,')'
            CHMSG = 'Unconnected Triangle: ' // BDUM(1:NC)
            CALL WRMSGS( INDX )
          ENDIF
        ENDDO
!
!---    Determine triangle areas and centroids  ---
!
        VOL_FRCX = 0.D+0
        DO MTX = 1,NTP_FRC(NFX)
          NTX = MTX + IP_FRC(1,NFX) - 1
          CALL TRGAREA( XE_FRC(1,NTX),YE_FRC(1,NTX),ZE_FRC(1,NTX),
     &      XE_FRC(2,NTX),YE_FRC(2,NTX),ZE_FRC(2,NTX),
     &      XE_FRC(3,NTX),YE_FRC(3,NTX),ZE_FRC(3,NTX),AF_FRC(NTX) )
          CALL PGCNTRD( 3,XE_FRC(1,NTX),YE_FRC(1,NTX),ZE_FRC(1,NTX),
     &      XP_FRC(NTX),YP_FRC(NTX),ZP_FRC(NTX) )
          VOL_FRCX = VOL_FRCX + AF_FRC(NTX)*JRC_FRC(NTX)
        ENDDO
        IP_FRC(2,NFX) = NTX
        IF( JNDX.EQ.1 ) THEN
          WRITE(36,'(/,4X,A,1PE12.5)')'Fracture Volume, m^3: ',VOL_FRCX
        ELSEIF( JNDX.EQ.2 ) THEN
          WRITE(36,'(/,4X,A,1PE12.5)') 'Fault Volume, m^3: ',VOL_FRCX
        ENDIF     
        IF( IFILE.EQ.1 ) THEN
          IRD = IRDX
          IFILE = 0
          CLOSE(UNIT=35)
        ENDIF
      ENDDO L1
      NT_FRAC = NTX
      IF( NT_FRAC.GT.LT_FRC ) THEN
        INDX = 5
        IF( JNDX.EQ.1 ) THEN
          CHMSG = 'Number of Fracture Triangles > Parameter LT_FRC'
        ELSEIF( JNDX.EQ.2 ) THEN
          CHMSG = 'Number of Fault Triangles > Parameter LT_FRC'
        ENDIF     
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Surface normals of the fracture triangles  ---
!
      IF( JNDX.EQ.1 ) THEN
        WRITE(36,'(/,A,I9)') 'Fracture to Node Connections'
      ELSEIF( JNDX.EQ.2 ) THEN
        WRITE(36,'(/,A,I9)') 'Fault to Node Connections'
      ENDIF     
      NFNCX = 0
      NAT_FRC = 0
      DO NFX = 1,NF_FRC
        IF( JNDX.EQ.1 ) THEN
          WRITE(36,'(/,2X,A,I9)') 'Fracture #: ',NFX
        ELSEIF( JNDX.EQ.2 ) THEN
          WRITE(36,'(/,2X,A,I9)') 'Fault #: ',NFX
        ENDIF     
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
          AFN_FRCX = 0.D+0
          IPN_FRC(1,NTX) = NFNCX + 1
!
!---      x,y,z coordinates of the fracture triangle vertices ---
!
          DO NVX = 1,3
            V1X(1,NVX) = XE_FRC(NVX,NTX)
            V1X(2,NVX) = YE_FRC(NVX,NTX)
            V1X(3,NVX) = ZE_FRC(NVX,NTX)
          ENDDO
!
!---      Surface parallel vectors of the fracture triangles ---
!
          VAX(1) = XE_FRC(2,NTX) - XE_FRC(1,NTX)
          VBX(1) = XE_FRC(3,NTX) - XE_FRC(1,NTX)
          VAX(2) = YE_FRC(2,NTX) - YE_FRC(1,NTX)
          VBX(2) = YE_FRC(3,NTX) - YE_FRC(1,NTX)
          VAX(3) = ZE_FRC(2,NTX) - ZE_FRC(1,NTX)
          VBX(3) = ZE_FRC(3,NTX) - ZE_FRC(1,NTX)
!
!---      Surface normal vector of the fracture triangle ---
!
          SIGN_FRC(1,NTX) = VAX(2)*VBX(3) - VAX(3)*VBX(2)
          SIGN_FRC(2,NTX) = VAX(3)*VBX(1) - VAX(1)*VBX(3)
          SIGN_FRC(3,NTX) = VAX(1)*VBX(2) - VAX(2)*VBX(1)
!
!---      Unit surface normal vector of the fracture triangle ---
!
          SNMX = SQRT( SIGN_FRC(1,NTX)**2 + 
     &      SIGN_FRC(2,NTX)**2 + SIGN_FRC(3,NTX)**2)
          SIGN_FRC(1,NTX) = SIGN_FRC(1,NTX)/SNMX
          SIGN_FRC(2,NTX) = SIGN_FRC(2,NTX)/SNMX
          SIGN_FRC(3,NTX) = SIGN_FRC(3,NTX)/SNMX
          PX(1) = SIGN_FRC(1,NTX)
          PX(2) = SIGN_FRC(2,NTX)
          PX(3) = SIGN_FRC(3,NTX)
!
!---      Rotate surface normal vector from global x,y,z coordinates
!         to local geomechanics x',y',z' coordinates
!
          CALL MATMUL( ROTMAT,PX,SIGN_FRC(1,NTX),3,3,1 )
!
!---      Unit surface parallel vector of the fracture triangle ---
!
          VAMX(1) = SQRT( VAX(1)**2 +  VAX(2)**2 + VAX(3)**2)
          VAX(1) = VAX(1)/VAMX(1)
          VAX(2) = VAX(2)/VAMX(1)
          VAX(3) = VAX(3)/VAMX(1)
!
!---      Unit orthogonal vector to surface normal and 
!         surface parallel vectors ---
!
          VDX(1) = VAX(2)*SIGN_FRC(3,NTX) - VAX(3)*SIGN_FRC(2,NTX)
          VDX(2) = VAX(3)*SIGN_FRC(1,NTX) - VAX(1)*SIGN_FRC(3,NTX)
          VDX(3) = VAX(1)*SIGN_FRC(2,NTX) - VAX(2)*SIGN_FRC(1,NTX)
!
!---      Normal stress (Pa) on fracture triangle---
!
!          CALL MATMUL( SIGN_FRC(1,NTX),SIGX,VCX,1,3,3 )
!          CALL MATMUL( VCX,SIGN_FRC(1,NTX),SIGNX,1,3,1 )
!
!---      Shear stress (Pa) on fracture triangle in direction of VAX ---
!
          CALL MATMUL( VAX,SIGX,VCX,1,3,3 )
          CALL MATMUL( VCX,SIGN_FRC(1,NTX),VAMX,1,3,1 )
!
!---      Shear stress (Pa) on fracture triangle in direction of VDX ---
!
          CALL MATMUL( VDX,SIGX,VCX,1,3,3 )
          CALL MATMUL( VCX,SIGN_FRC(1,NTX),VDMX(1),1,3,1 )
!
!---      Magnitude weighted normal and shear stress (Pa) on 
!         fracture triangle ---
!
          DO NVX = 1,3
!            SIGN_FRC(NVX,NTX) = SIGN_FRC(NVX,NTX)*SIGNX(1)
            SIGS_FRC(NVX,NTX) = VAX(NVX)*VAMX(1) + VDX(NVX)*VDMX(1)
          ENDDO
!
!---      For MPI preprocessing find the node with the closest
!         centroid to the fracture triangle centroid, and then
!         associate just that node with the fracture triangle. One
!         field node can be associated with multiple fracture triangles
!         but each fracture triangle is only associate with one
!         field node  ---
!
          IF( ISLC(67).EQ.1 ) THEN
            DISTX = 1.D+20
            ICHKX = 0
            NX = 0
            DO N = 1,NFLD
              IF( IXP(N).EQ.0 ) CYCLE
              IFLGX = 0
              XMNX = XE(1,N)
              XMAX = XE(1,N)
              YMNX = YE(1,N)
              YMAX = YE(1,N)
              ZMNX = ZE(1,N)
              ZMAX = ZE(1,N)
              DO M = 2,8
                XMNX = MIN( XMNX,XE(M,N) )
                XMAX = MAX( XMAX,XE(M,N) )
                YMNX = MIN( YMNX,YE(M,N) )
                YMAX = MAX( YMAX,YE(M,N) )
                ZMNX = MIN( ZMNX,ZE(M,N) )
                ZMAX = MAX( ZMAX,ZE(M,N) )
              ENDDO
              XMNX = XMNX - 1.D-6
              XMAX = XMAX + 1.D-6
              YMNX = YMNX - 1.D-6
              YMAX = YMAX + 1.D-6
              ZMNX = ZMNX - 1.D-6
              ZMAX = ZMAX + 1.D-6
!
!---          Check for triangle vertices within domain  ---
!
              DO NVX = 1,3
                IF( XE_FRC(NVX,NTX).LT.XMNX .OR. 
     &            XE_FRC(NVX,NTX).GT.XMAX .OR. 
     &            YE_FRC(NVX,NTX).LT.YMNX .OR. 
     &            YE_FRC(NVX,NTX).GT.YMAX .OR. 
     &            ZE_FRC(NVX,NTX).LT.ZMNX .OR. 
     &            ZE_FRC(NVX,NTX).GT.ZMAX ) CYCLE
                CALL WITH_IN( XE_FRC(NVX,NTX),YE_FRC(NVX,NTX),
     &            ZE_FRC(NVX,NTX),IFLGX,N )
                ICHKX = MAX( ICHKX,IFLGX )
              ENDDO
!
!---          Check for triangle centroid within domain  ---
!
              IF( XP_FRC(NTX).GE.XMNX .AND. 
     &          XP_FRC(NTX).LE.XMAX .AND. 
     &          YP_FRC(NTX).GE.YMNX .AND. 
     &          YP_FRC(NTX).LE.YMAX .AND. 
     &          ZP_FRC(NTX).GE.ZMNX .AND. 
     &          ZP_FRC(NTX).LE.ZMAX ) THEN
                CALL WITH_IN( XP_FRC(NTX),YP_FRC(NTX),
     &            ZP_FRC(NTX),IFLGX,N )
                ICHKX = MAX( ICHKX,IFLGX )
              ENDIF
              VARX = SQRT((XP_FRC(NTX)-XP(N))**2 + 
     &          (YP_FRC(NTX)-YP(N))**2 + (ZP_FRC(NTX)-ZP(N))**2)
              IF( VARX.LT.DISTX ) THEN
                DISTX = VARX
                NX = N
              ENDIF
            ENDDO
            IF( ICHKX.EQ.1 ) THEN
              DO NVZ = 1,3
                ND_FRC(NVZ,NTX) = NX
              ENDDO
!
!---        All three fracture triangle vertices and centroid outside
!           of computational domain, make triangle inactive  ---
!
            ELSE
              NXP_FRC = NXP_FRC + 1
              CYCLE
            ENDIF
            IF( JNDX.EQ.1 ) THEN
              WRITE(36,'(2X,A,I9)') 'Fracture Triangle #: ',NTX
            ELSEIF( JNDX.EQ.2 ) THEN
              WRITE(36,'(2X,A,I9)') 'Fault Triangle #: ',NTX
            ENDIF     
            CALL TRGAREA( XE_FRC(1,NTX),YE_FRC(1,NTX),
     &        ZE_FRC(1,NTX),XE_FRC(2,NTX),YE_FRC(2,NTX),
     &        ZE_FRC(2,NTX),XE_FRC(3,NTX),YE_FRC(3,NTX),
     &        ZE_FRC(3,NTX),APX )
            NFNCX = NFNCX + 1
            IF( NFNCX.GT.LNC_FRC ) THEN
              IF( JNDX.EQ.1 ) THEN
                CHMSG = 'Number of Fracture Triangle to ' // 
     &            'Node Connections > Parameter LNC_FRC'
              ELSEIF( JNDX.EQ.2 ) THEN
                CHMSG = 'Number of Fault Triangle to ' // 
     &            'Node Connections > Parameter LNC_FRC'
              ENDIF     
              INDX = 5
              CALL WRMSGS( INDX )
            ENDIF
            VARX = 1.D+0
            INDX = 4
            IUNM = 1
            CALL RDUNIT(UNLN,VARX,INDX)
            INCM_FRC(NFNCX) = NX
            AFN_FRC(NFNCX) = APX
            AFN_FRCX = AFN_FRCX + APX
            CALL PGCNTRD( 3,XE_FRC(1,NTX),YE_FRC(1,NTX),ZE_FRC(1,NTX),
     &        XCNTX,YCNTX,ZCNTX )
            XPFN_FRC(NFNCX) = XCNTX
            YPFN_FRC(NFNCX) = YCNTX
            ZPFN_FRC(NFNCX) = ZCNTX
!
!---        Divide grid-cell into 512 sub-blocks and compute
!           average distance to fracture surface  ---
!
            DFN_FRC(NFNCX) = 5.D-1*(VOL(NX)**(1.D+0/3.D+0))
            ADMFX = 0.D+0
            DO KX = 1,8
            DO JX = 1,8
            DO IX = 1,8
              XDX = (REAL(IX)-5.D-1)/8.D+0
              YDX = (REAL(JX)-5.D-1)/8.D+0
              ZDX = (REAL(KX)-5.D-1)/8.D+0
              CALL TRILINEAR( XDX,YDX,ZDX,XE(1,NX),PX(1) ) 
              CALL TRILINEAR( XDX,YDX,ZDX,YE(1,NX),PX(2) ) 
              CALL TRILINEAR( XDX,YDX,ZDX,ZE(1,NX),PX(3) )
              CALL DIST_PT( V1X(1,1),V1X(1,2),V1X(1,3),PX,DISTX )
              ADMFX = ADMFX + DISTX
            ENDDO
            ENDDO
            ENDDO
            DFN_FRC(NFNCX) = ADMFX/5.12D+2
          ELSE
!
!---        Resident node number for fracture triangle vertices  ---
!
            DO NVX = 1,3
              DO N = 1,NFLD
                IF( IXP(N).EQ.0 ) CYCLE
                XMNX = XE(1,N)
                XMAX = XE(1,N)
                YMNX = YE(1,N)
                YMAX = YE(1,N)
                ZMNX = ZE(1,N)
                ZMAX = ZE(1,N)
                DO M = 2,8
                  XMNX = MIN( XMNX,XE(M,N) )
                  XMAX = MAX( XMAX,XE(M,N) )
                  YMNX = MIN( YMNX,YE(M,N) )
                  YMAX = MAX( YMAX,YE(M,N) )
                  ZMNX = MIN( ZMNX,ZE(M,N) )
                  ZMAX = MAX( ZMAX,ZE(M,N) )
                ENDDO
                XMNX = XMNX - 1.D-6
                XMAX = XMAX + 1.D-6
                YMNX = YMNX - 1.D-6
                YMAX = YMAX + 1.D-6
                ZMNX = ZMNX - 1.D-6
                ZMAX = ZMAX + 1.D-6
                IF( XE_FRC(NVX,NTX).LT.XMNX .OR. 
     &            XE_FRC(NVX,NTX).GT.XMAX .OR. 
     &            YE_FRC(NVX,NTX).LT.YMNX .OR. 
     &            YE_FRC(NVX,NTX).GT.YMAX .OR. 
     &            ZE_FRC(NVX,NTX).LT.ZMNX .OR. 
     &            ZE_FRC(NVX,NTX).GT.ZMAX ) CYCLE
                CALL WITHIN( XE_FRC(NVX,NTX),YE_FRC(NVX,NTX),
     &            ZE_FRC(NVX,NTX),IFLGX,N )
                IF( IFLGX.NE.0 ) THEN
                  ND_FRC(NVX,NTX) = N
                  EXIT
                ENDIF
              ENDDO
            ENDDO
!
!---        For fracture triangle vertices in multiple grid cells,
!           check whether a single grid cell is inclusive  ---
!
            IF( (ND_FRC(1,NTX).NE.ND_FRC(2,NTX) .OR. 
     &        ND_FRC(2,NTX).NE.ND_FRC(3,NTX)) .AND.
     &        ND_FRC(1,NTX)*ND_FRC(2,NTX)*ND_FRC(3,NTX).GT.0 ) THEN
              DO NVX = 1,3
                ICHKX = 1
                N = ND_FRC(NVX,NTX)
                DO NVZ = 1,3
                  CALL WITHIN( XE_FRC(NVZ,NTX),YE_FRC(NVZ,NTX),
     &              ZE_FRC(NVZ,NTX),IFLGX,N )
                  IF( IFLGX.EQ.0 ) ICHKX = 0
                ENDDO
                IF( ICHKX.EQ.1 ) THEN
                  DO NVZ = 1,3
                    ND_FRC(NVZ,NTX) = N
                  ENDDO
                  EXIT
                ENDIF
              ENDDO
            ENDIF
!
!---        Any fracture triangle vertex outside of 
!           computational domain, make triangle inactive  ---
!
            IF( ND_FRC(1,NTX).EQ.0 .OR. ND_FRC(2,NTX).EQ.0 .OR.
     &        ND_FRC(3,NTX).EQ.0 ) THEN
              NXP_FRC = NXP_FRC + 1
              CYCLE
!!
!!---        Fracture triangle vertex outside of computational domain  ---
!!
!            IF( ND_FRC(1,NTX).EQ.0 .OR. ND_FRC(2,NTX).EQ.0 .OR.
!       &      ND_FRC(3,NTX).EQ.0 ) THEN
!              INDX = 4
!              WRITE(FORM4(5:5),'(I1)') ICOUNT(NFX)
!              WRITE(FORM4(10:10),'(I1)') ICOUNT(NTX)
!              NC = ICOUNT(NFX) + ICOUNT(NTX) + 
!       &        ICOUNT(N) + 25
!              IF( JNDX.EQ.1 ) THEN
!                WRITE(BDUM(1:NC),FORM4) '(Fracture #:',NFX,' Triangle #:',
!       &          NTX,')'
!                CHMSG = 'Fracture Triangle Vertex Outside of Domain: ' 
!       &          // BDUM(1:NC)
!              ELSEIF( JNDX.EQ.2 ) THEN
!                WRITE(BDUM(1:NC),FORM4) '(Fault #:',NFX,' Triangle #:',
!       &          NTX,')'
!                CHMSG = 'Fault Triangle Vertex Outside of Domain: ' 
!       &          // BDUM(1:NC)
!              ENDIF     
!              CALL WRMSGS( INDX )
!
!---        Grid cell to fracture triangle contacts  ---
!
            ELSEIF( ND_FRC(1,NTX).NE.ND_FRC(2,NTX) .OR.
     &        ND_FRC(2,NTX).NE.ND_FRC(3,NTX) ) THEN
              IF( JNDX.EQ.1 ) THEN
                WRITE(36,'(2X,A,I9)') 'Fracture Triangle #: ',NTX
              ELSEIF( JNDX.EQ.2 ) THEN
                WRITE(36,'(2X,A,I9)') 'Fault Triangle #: ',NTX
              ENDIF     
!
!---          Fracture triangle contacts multiple nodes  ---
!
              IHI = 1
              ILO = IFLD
              JHI = 1
              JLO = JFLD
              KHI = 1
              KLO = KFLD
              DO NVX = 1,3
                IF( ND_FRC(NVX,NTX).EQ.0 ) CYCLE
                ILO = MIN( ILO,ID(ND_FRC(NVX,NTX)) )
                IHI = MAX( IHI,ID(ND_FRC(NVX,NTX)) )
                JLO = MIN( JLO,JD(ND_FRC(NVX,NTX)) )
                JHI = MAX( JHI,JD(ND_FRC(NVX,NTX)) )
                KLO = MIN( KLO,KD(ND_FRC(NVX,NTX)) )
                KHI = MAX( KHI,KD(ND_FRC(NVX,NTX)) )
              ENDDO
!
!---          Loop over the range of nodes  ---
!
              DO K = KLO,KHI
              DO J = JLO,JHI
              DO I = ILO,IHI
                NPTX = 0
                N = ND(I,J,K)
                IF( IXP(N).EQ.0 ) CYCLE
                DO NVX = 1,3
                  IF( ND_FRC(NVX,NTX).EQ.N ) THEN
                    NPTX = NPTX+1
                    PTX(1,NPTX) = XE_FRC(NVX,NTX)
                    PTX(2,NPTX) = YE_FRC(NVX,NTX)
                    PTX(3,NPTX) = ZE_FRC(NVX,NTX)
                  ENDIF
                ENDDO
!
!---            Loop over the surfaces of the node and the triangles
!               of the node surface, checking for intersections of 
!               the fracture triangle and node surface triangles  ---
!
                ICPX = 0
                DO NS1X = 1,6
                  DO NS2X = 1,4
                    MX = MSX(NS2X,NS1X)
!
!---                Cylindrical coordinates---
!
                    IF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
                      XPX(NS2X) = XE(MX,N)*COS(YE(MX,N))
                      YPX(NS2X) = XE(MX,N)*SIN(YE(MX,N))
                      ZPX(NS2X) = ZE(MX,N)
                    ELSE
                      XPX(NS2X) = XE(MX,N)
                      YPX(NS2X) = YE(MX,N)
                      ZPX(NS2X) = ZE(MX,N)
                    ENDIF                
                  ENDDO
                  NPX = 4
                  CALL PGCNTRD( NPX,XPX(1),YPX(1),ZPX(1),
     &              XPX(5),YPX(5),ZPX(5) )            
                  V2X(1,3) = XPX(5)
                  V2X(2,3) = YPX(5)
                  V2X(3,3) = ZPX(5)
                  DO NS2X = 1,4
!
!---                x,y,z coordinates of the node surface triangle 
!                   vertices ---
!
                    V2X(1,1) = XPX(N1X(NS2X))
                    V2X(2,1) = YPX(N1X(NS2X))
                    V2X(3,1) = ZPX(N1X(NS2X))
                    V2X(1,2) = XPX(N2X(NS2X))
                    V2X(2,2) = YPX(N2X(NS2X))
                    V2X(3,2) = ZPX(N2X(NS2X))
!
!---                Fast triangle-triangle intersection test, with 
!                   intersection point finder when intersection
!                   noted  ---
!
!                    CALL FTTIT( V1X,V2X,P1X,P2X,PN1X,PN2X,ICPX,ITX )
                    CALL TTIT3D(V1X,V2X,P1X,P2X,PN1X,PN2X,CPX,ICPX,ITX)
!
!---                Record intersection points  ---
!
                    IF( ITX.EQ.1 ) THEN
!
!---                  Skip duplicate points  ---
!
                      ISKPX = 0
                      DO M = 1,NPTX
                        IF( ABS(P1X(1)-PTX(1,M)).LT.EPSL .AND. 
     &                    ABS(P1X(2)-PTX(2,M)).LT.EPSL .AND.
     &                    ABS(P1X(3)-PTX(3,M)).LT.EPSL ) THEN
                          ISKPX = 1
                          EXIT
                        ENDIF
                      ENDDO
                      IF( ISKPX.EQ.0 ) THEN
                        NPTX = NPTX + 1
                        PTX(1,NPTX) = P1X(1)
                        PTX(2,NPTX) = P1X(2)
                        PTX(3,NPTX) = P1X(3)
                      ENDIF
!
!---                  Skip duplicate points  ---
!
                      ISKPX = 0
                      DO M = 1,NPTX
                        IF( ABS(P2X(1)-PTX(1,M)).LT.EPSL .AND. 
     &                    ABS(P2X(2)-PTX(2,M)).LT.EPSL .AND.
     &                    ABS(P2X(3)-PTX(3,M)).LT.EPSL ) THEN
                          ISKPX = 1
                          EXIT
                        ENDIF
                      ENDDO
                      IF( ISKPX.EQ.0 ) THEN
                        NPTX = NPTX + 1
                        PTX(1,NPTX) = P2X(1)
                        PTX(2,NPTX) = P2X(2)
                        PTX(3,NPTX) = P2X(3)
                      ENDIF
                    ENDIF
!
!---                Coplanar triangle sets identified, record
!                   intersection points  ---
!
                    IF( ICPX.NE.0 ) THEN
                      DO IC = 1,ICPX
!
!---                    Skip duplicate points  ---
!
                        ISKPX = 0
                        DO M = 1,NPTX
                          IF( ABS(CPX(1,IC)-PTX(1,M)).LT.EPSL .AND. 
     &                      ABS(CPX(2,IC)-PTX(2,M)).LT.EPSL .AND.
     &                      ABS(CPX(3,IC)-PTX(3,M)).LT.EPSL ) THEN
                            ISKPX = 1
                            EXIT
                          ENDIF
                        ENDDO
                        IF( ISKPX.EQ.0 ) THEN
                          NPTX = NPTX + 1
                          PTX(1,NPTX) = CPX(1,IC)
                          PTX(2,NPTX) = CPX(2,IC)
                          PTX(3,NPTX) = CPX(3,IC)
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
!
!---            Surface area of fracture triangle within grid-cell  ---
!
                IF( NPTX.GT.2 ) THEN
                  IF( NPTX.GT.3 ) THEN
                    DO M = 1,NPTX
                      XPTX(M) = PTX(1,M)
                      YPTX(M) = PTX(2,M)
                      ZPTX(M) = PTX(3,M)
                    ENDDO
!
!---                Check for nonzero surface area of polygon  ---
!
                    APX = 0.D+0
                    DO M = 1,NPTX-2
                      CALL TRGAREA( PTX(1,M),PTX(2,M),PTX(3,M),
     &                  PTX(1,M+1),PTX(2,M+1),PTX(3,M+1),
     &                  PTX(1,NPTX),PTX(2,NPTX),PTX(3,NPTX),ATX )
                      APX = APX + ATX
                    ENDDO
!
!---                Ignore relatively small surface areas   ---
!
                    IF( APX.LT.1.D-2*AF_FRC(NTX) ) CYCLE
                    CALL PGCNTRD(NPTX,XPTX,YPTX,ZPTX,XCNTX,YCNTX,ZCNTX)
!
!---                Angles between centroid to vertex 1 vector and
!                   centroid to vertex n vector  ---
!
                    ROTX(1) = 0.D+0
                    VPX(1) = PTX(1,1)-XCNTX
                    VPX(2) = PTX(2,1)-YCNTX
                    VPX(3) = PTX(3,1)-ZCNTX
                    VPLX = SQRT( VPX(1)**2 + VPX(2)**2 + VPX(3)**2 )
                    DO M = 2,NPTX
                      VQX(1) = PTX(1,M)-XCNTX
                      VQX(2) = PTX(2,M)-YCNTX
                      VQX(3) = PTX(3,M)-ZCNTX
                      VQLX = SQRT( VQX(1)**2 + VQX(2)**2 + VQX(3)**2 )
                      DPNX = VDOTP(VPX,VQX)/(VPLX*VQLX)
                      DPNX = MAX( MIN( DPNX,1.D+0 ),-1.D+0 )
                      ROTX(M) = ACOS( DPNX )
                      CALL VCROSSP( VQX,VPX,VNX )
                      ROTX(M) = SIGN( ROTX(M),VDOTP(VNX,PN1X) )
                      IF( ROTX(M).LT.0.D+0 ) ROTX(M) = 2.D+0*GPI+ROTX(M)
                    ENDDO
!
!---                Sort vertices by rotation angles  ---
!
                    DO M = 2,NPTX
                      RTX = ROTX(M)
                      PT1X = PTX(1,M)
                      PT2X = PTX(2,M)
                      PT3X = PTX(3,M)
                      DO L = M-1,1,-1
                        IF( ROTX(L).LE.RTX ) GOTO 10
                        ROTX(L+1) = ROTX(L)
                        PTX(1,L+1) = PTX(1,L)
                        PTX(2,L+1) = PTX(2,L)
                        PTX(3,L+1) = PTX(3,L)
                      ENDDO
                      L = 0
   10                 CONTINUE
                      ROTX(L+1) = RTX
                      PTX(1,L+1) = PT1X
                      PTX(2,L+1) = PT2X
                      PTX(3,L+1) = PT3X
                    ENDDO
                    NPTX = NPTX + 1
                    PTX(1,NPTX) = XCNTX
                    PTX(2,NPTX) = YCNTX
                    PTX(3,NPTX) = ZCNTX
                  ENDIF
                  APX = 0.D+0
                  DO M = 1,NPTX-2
                    CALL TRGAREA( PTX(1,M),PTX(2,M),PTX(3,M),
     &                PTX(1,M+1),PTX(2,M+1),PTX(3,M+1),
     &                PTX(1,NPTX),PTX(2,NPTX),PTX(3,NPTX),ATX )
                    APX = APX + ATX
                  ENDDO
                  IF( NPTX.GT.3 ) THEN
                    CALL TRGAREA( PTX(1,1),PTX(2,1),PTX(3,1),
     &                PTX(1,NPTX-1),PTX(2,NPTX-1),PTX(3,NPTX-1),
     &                PTX(1,NPTX),PTX(2,NPTX),PTX(3,NPTX),ATX )
                    APX = APX + ATX
                  ENDIF
!              IF( NTX.EQ.7 ) THEN
!                PRINT *,'Node = ',N
!                DO M = 1,NPTX
!                  PRINT *,PTX(1,M),PTX(2,M),PTX(3,M)
!                ENDDO
!              ENDIF
!
!---              Ignore relatively small surface areas   ---
!
                  IF( APX.GE.1.D-2*AF_FRC(NTX) ) THEN
                    NFNCX = NFNCX + 1
!                    WRITE(16,'(4(A,I9),2(A,1PE12.5))') 
!       &              '1: Fracture #: ',NFX,
!       &              ' Triangle #: ',NTX,
!       &              ' Node #: ',N,' Connection #: ',NFNCX,
!       &              ' Area,m^2: ',APX,
!       &              ' Triangle Area,m^2: ',AF_FRC(NTX)
                    IF( NFNCX.GT.LNC_FRC ) THEN
                      IF( JNDX.EQ.1 ) THEN
                        CHMSG = 'Number of Fracture Triangle to ' // 
     &                    'Node Connections > Parameter LNC_FRC'
                      ELSEIF( JNDX.EQ.2 ) THEN
                        CHMSG = 'Number of Fault Triangle to ' // 
     &                    'Node Connections > Parameter LNC_FRC'
                      ENDIF     
                      INDX = 5
                      CALL WRMSGS( INDX )
                    ENDIF
                    VARX = 1.D+0
                    INDX = 4
                    IUNM = 1
                    CALL RDUNIT(UNLN,VARX,INDX)
                    NCH = INDEX( UNLN(1:),'  ')-1
                    INCM_FRC(NFNCX) = N
                    IF( ICPX.NE.0 ) APX = 5.D-1*APX
                    AFN_FRC(NFNCX) = APX
                    AFN_FRCX = AFN_FRCX + APX
                    DO M = 1,NPTX
                      XPTX(M) = PTX(1,M)
                      YPTX(M) = PTX(2,M)
                      ZPTX(M) = PTX(3,M)
                    ENDDO
                    CALL PGCNTRD(NPTX,XPTX,YPTX,ZPTX,XCNTX,YCNTX,ZCNTX)
                    XPFN_FRC(NFNCX) = XCNTX
                    YPFN_FRC(NFNCX) = YCNTX
                    ZPFN_FRC(NFNCX) = ZCNTX
!
!---                Divide grid-cell into 512 sub-blocks and compute
!                   average distance to fracture surface  ---
!
                    DFN_FRC(NFNCX) = 5.D-1*(VOL(N)**(1.D+0/3.D+0))
                    ADMFX = 0.D+0
                    DO KX = 1,8
                    DO JX = 1,8
                    DO IX = 1,8
                      XDX = (REAL(IX)-5.D-1)/8.D+0
                      YDX = (REAL(JX)-5.D-1)/8.D+0
                      ZDX = (REAL(KX)-5.D-1)/8.D+0
                      CALL TRILINEAR( XDX,YDX,ZDX,XE(1,N),PX(1) ) 
                      CALL TRILINEAR( XDX,YDX,ZDX,YE(1,N),PX(2) ) 
                      CALL TRILINEAR( XDX,YDX,ZDX,ZE(1,N),PX(3) )
                      CALL DIST_PT( V1X(1,1),V1X(1,2),V1X(1,3),
     &                  PX,DISTX )
                      ADMFX = ADMFX + DISTX
                    ENDDO
                    ENDDO
                    ENDDO
                    DFN_FRC(NFNCX) = ADMFX/5.12D+2
!                    IF( ICPX.NE.0 ) THEN
!                      IF( ISLC(35).EQ. 1) 
!       &                WRITE(IWR,'(4X,A,I9,2(3A,1PE11.4,A,1PE11.4),A)') 
!       &                'Node #: ',N,'; Coplanar Area, ',
!       &                UNLN(1:NCH),'^2 : ',(AFN_FRC(NFNCX)*VARX*VARX),
!       &                ' (',AFN_FRC(NFNCX),', m^2); Distance, ',
!       &                UNLN(1:NCH),' : ',(DFN_FRC(NFNCX)*VARX),
!       &                ' (',DFN_FRC(NFNCX),', m)'
!                    ELSE
!                      IF( ISLC(35).EQ. 1) 
!       &                WRITE(IWR,'(4X,A,I9,2(3A,1PE11.4,A,1PE11.4),A)') 
!       &                'Node #: ',N,'; Area, ',
!       &                UNLN(1:NCH),'^2 : ',(AFN_FRC(NFNCX)*VARX*VARX),
!       &                ' (',AFN_FRC(NFNCX),', m^2); Distance, ',
!       &                UNLN(1:NCH),' : ',(DFN_FRC(NFNCX)*VARX),
!       &                ' (',DFN_FRC(NFNCX),', m)'
!                    ENDIF
!
!---                For coplanar fractures, link to adjacent node  ---
!
                    IF( ICPX.NE.0 ) THEN
                      NFNCX = NFNCX + 1
!                      WRITE(16,'(4(A,I9),2(A,1PE12.5))') 
!       &                '2: Fracture #: ',NFX,
!       &                ' Triangle #: ',NTX,
!       &                ' Node #: ',N,' Connection #: ',NFNCX,
!       &                ' Area,m^2: ',APX,
!       &                ' Triangle Area,m^2: ',AF_FRC(NTX)
                      IF( NFNCX.GT.LNC_FRC ) THEN
                        IF( JNDX.EQ.1 ) THEN
                          CHMSG = 'Number of Fracture Triangle to ' // 
     &                      'Node Connections > Parameter LNC_FRC'
                        ELSEIF( JNDX.EQ.2 ) THEN
                          CHMSG = 'Number of Fault Triangle to ' // 
     &                      'Node Connections > Parameter LNC_FRC'
                        ENDIF     
                        INDX = 5
                        CALL WRMSGS( INDX )
                      ENDIF
                      VARX = 1.D+0
                      INDX = 4
                      IUNM = 1
                      CALL RDUNIT(UNLN,VARX,INDX)
                      NCH = INDEX( UNLN(1:),'  ')-1
                      NAX = ICM(1,NS1X,N)
                      INCM_FRC(NFNCX) = NAX
                      AFN_FRC(NFNCX) = APX
                      AFN_FRCX = AFN_FRCX + APX
                      XPFN_FRC(NFNCX) = XCNTX
                      YPFN_FRC(NFNCX) = YCNTX
                      ZPFN_FRC(NFNCX) = ZCNTX
!
!---                  Divide grid-cell into 512 sub-blocks and compute
!                     average distance to fracture surface  ---
!
                      DFN_FRC(NFNCX) = 5.D-1*(VOL(NAX)**(1.D+0/3.D+0))
                      ADMFX = 0.D+0
                      DO KX = 1,8
                      DO JX = 1,8
                      DO IX = 1,8
                        XDX = (REAL(IX)-5.D-1)/8.D+0
                        YDX = (REAL(JX)-5.D-1)/8.D+0
                        ZDX = (REAL(KX)-5.D-1)/8.D+0
                        CALL TRILINEAR( XDX,YDX,ZDX,XE(1,NAX),PX(1) ) 
                        CALL TRILINEAR( XDX,YDX,ZDX,YE(1,NAX),PX(2) ) 
                        CALL TRILINEAR( XDX,YDX,ZDX,ZE(1,NAX),PX(3) )
                        CALL DIST_PT( V1X(1,1),V1X(1,2),V1X(1,3),
     &                    PX,DISTX )
                        ADMFX = ADMFX + DISTX
                      ENDDO
                      ENDDO
                      ENDDO
                      DFN_FRC(NFNCX) = ADMFX/5.12D+2
!                      IF( ISLC(35).EQ. 1) 
!       &                WRITE(IWR,'(4X,A,I9,2(3A,1PE11.4,A,1PE11.4),A)') 
!       &                'Node #: ',NAX,' Coplanar Area, ',
!       &                UNLN(1:NCH),'^2 : ',(AFN_FRC(NFNCX)*VARX*VARX),
!       &                ' (',AFN_FRC(NFNCX),', m^2); Distance, ',
!       &                UNLN(1:NCH),' : ',(DFN_FRC(NFNCX)*VARX),
!       &                ' (',DFN_FRC(NFNCX),', m)'  
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
              ENDDO
              ENDDO
!
!---        Fracture triangle completely within grid-cell, record
!           points and surface area  ---
!
            ELSE
              IF( JNDX.EQ.1 ) THEN
                WRITE(36,'(2X,A,I9)') 'Fracture Triangle #: ',NTX
              ELSEIF( JNDX.EQ.2 ) THEN
                WRITE(36,'(2X,A,I9)') 'Fault Triangle #: ',NTX
              ENDIF     
              CALL TRGAREA( XE_FRC(1,NTX),YE_FRC(1,NTX),
     &          ZE_FRC(1,NTX),XE_FRC(2,NTX),YE_FRC(2,NTX),
     &          ZE_FRC(2,NTX),XE_FRC(3,NTX),YE_FRC(3,NTX),
     &          ZE_FRC(3,NTX),APX )
              IF( NTX.EQ.7 ) THEN
                PRINT *,'Node = ',N
                DO M = 1,3
                  PRINT *,XE_FRC(M,NTX),YE_FRC(M,NTX),ZE_FRC(M,NTX)
                ENDDO
              ENDIF
!
!---          Ignore relatively small surface areas   ---
!
              IF( APX.GE.1.D-2*AF_FRC(NTX) ) THEN
                NFNCX = NFNCX + 1
!                WRITE(16,'(4(A,I9),2(A,1PE12.5))') 
!       &          '3: Fracture #: ',NFX,
!       &          ' Triangle #: ',NTX,
!       &          ' Node #: ',N,' Connection #: ',NFNCX,
!       &          ' Area,m^2: ',APX,
!       &          ' Triangle Area,m^2: ',AF_FRC(NTX)
                IF( NFNCX.GT.LNC_FRC ) THEN
                  IF( JNDX.EQ.1 ) THEN
                    CHMSG = 'Number of Fracture Triangle to ' // 
     &                'Node Connections > Parameter LNC_FRC'
                  ELSEIF( JNDX.EQ.2 ) THEN
                    CHMSG = 'Number of Fault Triangle to ' // 
     &                'Node Connections > Parameter LNC_FRC'
                  ENDIF     
                  INDX = 5
                  CALL WRMSGS( INDX )
                ENDIF
                VARX = 1.D+0
                INDX = 4
                IUNM = 1
                CALL RDUNIT(UNLN,VARX,INDX)
                INCM_FRC(NFNCX) = N
                AFN_FRC(NFNCX) = APX
                AFN_FRCX = AFN_FRCX + APX
                CALL PGCNTRD( 3,XE_FRC(1,NTX),YE_FRC(1,NTX),
     &            ZE_FRC(1,NTX),XCNTX,YCNTX,ZCNTX )
                XPFN_FRC(NFNCX) = XCNTX
                YPFN_FRC(NFNCX) = YCNTX
                ZPFN_FRC(NFNCX) = ZCNTX
!
!---            Divide grid-cell into 512 sub-blocks and compute
!               average distance to fracture surface  ---
!
                DFN_FRC(NFNCX) = 5.D-1*(VOL(N)**(1.D+0/3.D+0))
                ADMFX = 0.D+0
                DO KX = 1,8
                DO JX = 1,8
                DO IX = 1,8
                  XDX = (REAL(IX)-5.D-1)/8.D+0
                  YDX = (REAL(JX)-5.D-1)/8.D+0
                  ZDX = (REAL(KX)-5.D-1)/8.D+0
                  CALL TRILINEAR( XDX,YDX,ZDX,XE(1,N),PX(1) ) 
                  CALL TRILINEAR( XDX,YDX,ZDX,YE(1,N),PX(2) ) 
                  CALL TRILINEAR( XDX,YDX,ZDX,ZE(1,N),PX(3) )
                  CALL DIST_PT( V1X(1,1),V1X(1,2),V1X(1,3),PX,DISTX )
                  ADMFX = ADMFX + DISTX
                ENDDO
                ENDDO
                ENDDO
                DFN_FRC(NFNCX) = ADMFX/5.12D+2
              ENDIF
            ENDIF
          ENDIF
          NAT_FRC = NAT_FRC + 1
          IXP_FRC(NTX) = NAT_FRC
          IPN_FRC(2,NTX) = NFNCX
          DO MFNCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
            IF( AFN_FRCX.GT.AF_FRC(NTX) ) THEN
              SFX = MIN( 1.D+0,AF_FRC(NTX)/AFN_FRCX )
              AFN_FRC(MFNCX) = AFN_FRC(MFNCX)*SFX
              AFN_FRCX = AF_FRC(NTX)
            ENDIF
            WRITE(36,'(4X,A,I9,2(3A,1PE11.4,A,1PE11.4),A)') 
     &        'Node #: ',INCM_FRC(MFNCX),' Area, ',
     &         UNLN(1:NCH),'^2 : ',(AFN_FRC(MFNCX)*VARX*VARX),
     &         ' (',AFN_FRC(MFNCX),', m^2); Distance, ',
     &         UNLN(1:NCH),' : ',(DFN_FRC(MFNCX)*VARX),
     &         ' (',DFN_FRC(MFNCX),', m)'
          ENDDO
          IF( JNDX.EQ.1 ) THEN
            WRITE(36,'(4X,A,1PE11.4)') 
     &        'Total Fracture Triangle Area, m^2 :',AF_FRC(NTX)
            WRITE(36,'(4X,A,1PE11.4)') 
     &        'Effective Fracture Triangle Area, m^2 :',AFN_FRCX
            AF_FRC(NTX) = AFN_FRCX
          ELSEIF( JNDX.EQ.2 ) THEN
            WRITE(36,'(4X,A,1PE11.4)') 
     &        'Total Fault Triangle Area, m^2 :',AF_FRC(NTX)
            WRITE(36,'(4X,A,1PE11.4)') 
     &        'Effective Fault Triangle Area, m^2 :',AFN_FRCX
            AF_FRC(NTX) = AFN_FRCX
          ENDIF
        ENDDO
      ENDDO
!
!---  Fracture to fracture connections  ---
!
      IF( JNDX.EQ.1 ) THEN
        WRITE(36,'(/,A,I9)') 'Fracture to Fracture Connections'
      ELSEIF( JNDX.EQ.2 ) THEN
        WRITE(36,'(/,A,I9)') 'Fault to Fault Connections'
      ENDIF     
      NFFCX = 0
!
!---  Loop over fractures  ---
!
      DO NF1X = 1,NF_FRC
        IF( JNDX.EQ.1 ) THEN
          WRITE(36,'(/,2X,A,I9)') 'Fracture #: ',NF1X
        ELSEIF( JNDX.EQ.2 ) THEN
          WRITE(36,'(/,2X,A,I9)') 'Fault #: ',NF1X
        ENDIF     
!
!---    Loop over fracture triangles  ---
!
        DO NT1X = IP_FRC(1,NF1X),IP_FRC(2,NF1X)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NT1X).EQ.0 ) CYCLE
          IF( JNDX.EQ.1 ) THEN
            WRITE(36,'(2X,A,I9)') 'Fracture Triangle #: ',NT1X
          ELSEIF( JNDX.EQ.2 ) THEN
            WRITE(36,'(2X,A,I9)') 'Fault Triangle #: ',NT1X
          ENDIF     
          IPF_FRC(1,NT1X) = NFFCX + 1
!
!---      Loop over fracture triangles  ---
!
          DO NT2X = IP_FRC(1,NF1X),IP_FRC(2,NF1X)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            IF( NT1X.EQ.NT2X ) CYCLE
            NVVCX = 0
!
!---        Loop over vertices on first fracture triangle  ---
!
            DO IVP1X = 1,3
!
!---        Loop over vertices on second fracture triangle  ---
!
              DO IVP2X = 1,3
                DVX = ABS(XE_FRC(IVP1X,NT1X) - XE_FRC(IVP2X,NT2X)) +
     &            ABS(YE_FRC(IVP1X,NT1X) - YE_FRC(IVP2X,NT2X)) + 
     &            ABS(ZE_FRC(IVP1X,NT1X) - ZE_FRC(IVP2X,NT2X)) 
!
!---            Matching vertices found  ---
!
                IF( DVX.LE.1.D-9 ) THEN
                  NVVCX = NVVCX + 1
                  IVPX(NVVCX) = IVP1X
!
!---              Overlapping fracture triangles found  ---
!
                  IF( NVVCX.EQ.3 ) THEN
                    INDX = 4
                    WRITE(FORM4(5:5),'(I1)') ICOUNT(NF1X)
                    WRITE(FORM4(10:10),'(I1)') ICOUNT(NT1X)
                    NC = ICOUNT(NF1X) + ICOUNT(NT1X) + 
     &                ICOUNT(N) + 25
                    IF( JNDX.EQ.1 ) THEN
                      WRITE(BDUM(1:NC),FORM4) '(Fracture #:',NF1X,
     &                  ' Triangle #:',NT1X,')'
                      CHMSG = 'Overlapping Fracture Triangle: ' 
     &                  // BDUM(1:NC)
                    ELSEIF( JNDX.EQ.2 ) THEN
                      WRITE(BDUM(1:NC),FORM4) '(Fault #:',NF1X,
     &                  ' Triangle #:',NT1X,')'
                      CHMSG = 'Overlapping Fault Triangle: ' 
     &                  // BDUM(1:NC)
                    ENDIF     
                    CALL WRMSGS( INDX )
                  ENDIF
                  EXIT
                ENDIF
              ENDDO
            ENDDO
!
!---        Connecting fracture triangle found  ---
!
            IF( NVVCX.EQ.2 ) THEN
              NFFCX = NFFCX + 1
              AFF_FRC(NFFCX) = SQRT( 
     &          ((XE_FRC(IVPX(1),NT1X)-XE_FRC(IVPX(2),NT1X))**2) + 
     &          ((YE_FRC(IVPX(1),NT1X)-YE_FRC(IVPX(2),NT1X))**2) +
     &          ((ZE_FRC(IVPX(1),NT1X)-ZE_FRC(IVPX(2),NT1X))**2) )
              P3X(1) = 5.D-1*(XE_FRC(IVPX(1),NT1X)+XE_FRC(IVPX(2),NT1X))
              P3X(2) = 5.D-1*(YE_FRC(IVPX(1),NT1X)+YE_FRC(IVPX(2),NT1X))
              P3X(3) = 5.D-1*(ZE_FRC(IVPX(1),NT1X)+ZE_FRC(IVPX(2),NT1X))
              DFFM_FRC(NFFCX) = SQRT( 
     &          ((XP_FRC(NT1X)-P3X(1))**2) + 
     &          ((YP_FRC(NT1X)-P3X(2))**2) +
     &          ((ZP_FRC(NT1X)-P3X(3))**2) )
              DFF_FRC(NFFCX) = DFFM_FRC(NFFCX) + 
     &          SQRT( (XP_FRC(NT2X)-P3X(1))**2 + 
     &          (YP_FRC(NT2X)-P3X(2))**2 + (ZP_FRC(NT2X)-P3X(3))**2 )
              VARX = SQRT( 
     &          ((XP_FRC(NT1X)-XP_FRC(NT2X))**2) + 
     &          ((YP_FRC(NT1X)-YP_FRC(NT2X))**2) )
!              IF( VARX.LT.EPSL ) THEN
!                GRVF_FRC(NFFCX) = GRAV
!              ELSE
!                GRVF_FRC(NFFCX) = GRAV*SIN(ATAN((ZP_FRC(NT2X) -
!     &            ZP_FRC(NT1X))/VARX))
!              ENDIF
              VARX = 1.D+0
              INDX = 4
              IUNM = 1
              CALL RDUNIT(UNLN,VARX,INDX)
              IF( JNDX.EQ.1 ) THEN
                WRITE(36,'(4X,A,I9,2(3A,1PE11.4,A,1PE11.4),A)') 
     &          'Connecting Fracture Triangle #: ',NT2X,
     &          ' Connection Length, ',UNLN(1:NCH),' : ',
     &          (AFF_FRC(NFFCX)*VARX),' (',AFF_FRC(NFFCX),
     &          ', m); Inter-Fracture Distance, ',UNLN(1:NCH),' : ',
     &          (DFF_FRC(NFFCX)*VARX),' (',DFF_FRC(NFFCX),', m)' 
              ELSEIF( JNDX.EQ.2 ) THEN
                WRITE(36,'(4X,A,I9,2(3A,1PE11.4,A,1PE11.4),A)') 
     &          'Connecting Fault Triangle #: ',NT2X,
     &          ' Connection Length, ',UNLN(1:NCH),' : ',
     &          (AFF_FRC(NFFCX)*VARX),' (',AFF_FRC(NFFCX),
     &          ', m); Inter-Fault Distance, ',UNLN(1:NCH),' : ',
     &          (DFF_FRC(NFFCX)*VARX),' (',DFF_FRC(NFFCX),', m)' 
              ENDIF     
              ITCM_FRC(NFFCX) = NT2X
            ENDIF
          ENDDO
!
!---      Check for fracture triangle intersections with other
!         fractures  ---
!
          DO M = 1,3
            V1X(1,M) = XE_FRC(M,NT1X)
            V1X(2,M) = YE_FRC(M,NT1X)
            V1X(3,M) = ZE_FRC(M,NT1X)
          ENDDO
!
!---      Loop over fractures  ---
!
          DO NF2X = 1,NF_FRC
            IF( NF1X.EQ.NF2X ) CYCLE
!
!---        Loop over fracture triangles  ---
!
            DO NT2X = IP_FRC(1,NF2X),IP_FRC(2,NF2X)
              IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
              IF( NT1X.EQ.NT2X ) CYCLE
              DO M = 1,3
                V2X(1,M) = XE_FRC(M,NT2X)
                V2X(2,M) = YE_FRC(M,NT2X)
                V2X(3,M) = ZE_FRC(M,NT2X)
              ENDDO
!
!---          Fast triangle-triangle intersection test, with 
!             intersection point finder when intersection
!             noted  ---
!
              ICPX = 0
              ITX = 0
!              CALL FTTIT( V1X,V2X,P1X,P2X,PN1X,PN2X,ICPX,ITX )
              CALL TTIT3D( V1X,V2X,P1X,P2X,PN1X,PN2X,CPX,ICPX,ITX )
              IF( ITX.EQ.1 ) THEN
                NFFCX = NFFCX + 1
                AFF_FRC(NFFCX) = SQRT( (P1X(1)-P2X(1))**2 + 
     &            (P1X(2)-P2X(2))**2 + (P1X(3)-P2X(3))**2 )
                DO M = 1,3
                  P3X(M) = 5.D-1*(P1X(M)+P2X(M))
                ENDDO
                DFFM_FRC(NFFCX) = SQRT( 
     &            ((XP_FRC(NT1X)-P3X(1))**2) + 
     &            ((YP_FRC(NT1X)-P3X(2))**2) +
     &            ((ZP_FRC(NT1X)-P3X(3))**2) )
                DFF_FRC(NFFCX) = DFFM_FRC(NFFCX) + 
     &            SQRT( (XP_FRC(NT2X)-P3X(1))**2 + 
     &            (YP_FRC(NT2X)-P3X(2))**2 + (ZP_FRC(NT2X)-P3X(3))**2 )
                VARX = 1.D+0
                INDX = 4
                IUNM = 1
                CALL RDUNIT(UNLN,VARX,INDX)
                IF( JNDX.EQ.1 ) THEN
                  WRITE(36,'(4X,A,I9,2(3A,1PE11.4,A,1PE11.4),A)') 
     &            'Intersecting Fracture Triangle #: ',NT2X,
     &            ' Connection Length, ',UNLN(1:NCH),' : ',
     &            (AFF_FRC(NFFCX)*VARX),' (',AFF_FRC(NFFCX),
     &            ', m); Inter-Fracture Distance, ',UNLN(1:NCH),' : ',
     &            (DFF_FRC(NFFCX)*VARX),' (',DFF_FRC(NFFCX),', m)' 
                ELSEIF( JNDX.EQ.2 ) THEN
                  WRITE(36,'(4X,A,I9,2(3A,1PE11.4,A,1PE11.4),A)') 
     &            'Intersecting Fault Triangle #: ',NT2X,
     &            ' Connection Length, ',UNLN(1:NCH),' : ',
     &            (AFF_FRC(NFFCX)*VARX),' (',AFF_FRC(NFFCX),
     &            ', m); Inter-Fault Distance, ',UNLN(1:NCH),' : ',
     &            (DFF_FRC(NFFCX)*VARX),' (',DFF_FRC(NFFCX),', m)' 
                ENDIF     
                ITCM_FRC(NFFCX) = NT2X
              ELSEIF( ICPX.GT.0 ) THEN
                NFFCX = NFFCX + 1
                DO I = 1,ICPX
                  XPX(I) = CPX(1,I)
                  YPX(I) = CPX(2,I)
                  ZPX(I) = CPX(3,I)
                ENDDO
                CALL PGCNTRD( ICPX,XPX,YPX,ZPX,P3X(1),P3X(2),P3X(3) )
                AFF_FRC(NFFCX) = 0.D+0
                DO I = 1,ICPX-1
                  CALL TRGAREA( XPX(I),YPX(I),ZPX(I),XPX(I+1),YPX(I+1),
     &              ZPX(I+1),P3X(1),P3X(2),P3X(3),AFF_FRCX )
                  AFF_FRC(NFFCX) =  AFF_FRC(NFFCX) + AFF_FRCX
                ENDDO
                CALL TRGAREA( XPX(1),YPX(1),ZPX(1),XPX(ICPX),YPX(ICPX),
     &            ZPX(ICPX),P3X(1),P3X(2),P3X(3),AFF_FRCX )
                AFF_FRC(NFFCX) =  AFF_FRC(NFFCX) + AFF_FRCX
                DFFM_FRC(NFFCX) = SQRT( 
     &            ((XP_FRC(NT1X)-P3X(1))**2) + 
     &            ((YP_FRC(NT1X)-P3X(2))**2) +
     &            ((ZP_FRC(NT1X)-P3X(3))**2) )
                DFF_FRC(NFFCX) = DFFM_FRC(NFFCX) + 
     &            SQRT( (XP_FRC(NT2X)-P3X(1))**2 + 
     &            (YP_FRC(NT2X)-P3X(2))**2 + (ZP_FRC(NT2X)-P3X(3))**2 )
                VARX = 1.D+0
                INDX = 4
                IUNM = 1
                CALL RDUNIT(UNLN,VARX,INDX)
                IF( JNDX.EQ.1 ) THEN
                  WRITE(36,'(4X,A,I9,2(3A,1PE11.4,A,1PE11.4),A)') 
     &            'Intersecting Fracture Triangle #: ',NT2X,
     &            ' Connection Length, ',UNLN(1:NCH),' : ',
     &            (AFF_FRC(NFFCX)*VARX),' (',AFF_FRC(NFFCX),
     &            ', m); Inter-Fracture Distance, ',UNLN(1:NCH),' : ',
     &            (DFF_FRC(NFFCX)*VARX),' (',DFF_FRC(NFFCX),', m)' 
                ELSEIF( JNDX.EQ.2 ) THEN
                  WRITE(36,'(4X,A,I9,2(3A,1PE11.4,A,1PE11.4),A)') 
     &            'Intersecting Fault Triangle #: ',NT2X,
     &            ' Connection Length, ',UNLN(1:NCH),' : ',
     &            (AFF_FRC(NFFCX)*VARX),' (',AFF_FRC(NFFCX),
     &            ', m); Inter-Fault Distance, ',UNLN(1:NCH),' : ',
     &            (DFF_FRC(NFFCX)*VARX),' (',DFF_FRC(NFFCX),', m)' 
                ENDIF     
                ITCM_FRC(NFFCX) = NT2X
              ENDIF
            ENDDO
          ENDDO
          IPF_FRC(2,NT1X) = NFFCX
        ENDDO
      ENDDO
!!
!!---  Write a fracture data file for Tecplot  ---
!!
!      VAR = 1.D+0
!      INDX = 4
!      IUNM = 1
!      CALL RDUNIT(UNLN,VAR,INDX)
!      NCH = INDEX(UNLN(1:),'  ')-1
!      OPEN(UNIT=26, FILE='fracture.dat', STATUS='UNKNOWN', 
!     &  FORM='FORMATTED')
!      CLOSE(UNIT=26,STATUS='DELETE')
!      OPEN(UNIT=26, FILE='fracture.dat', STATUS='NEW', FORM='FORMATTED')
!      WRITE(26,'(A)') 'Title = "Fractures"'
!      WRITE(26,'(23A)') 'VARIABLES = "X, ',UNLN(1:NCH),'" ',
!     &  '"Y, ',UNLN(1:NCH),'" ','"Z, ',UNLN(1:NCH),'" ',
!     &  '"Area, ',UNLN(1:NCH),'^2" ',
!     &  '"Aperture, ',UNLN(1:NCH),'" ',
!     &  '"Triangle Number" ',
!     &  '"Normal Stress, Pa" ',
!     &  '"X-Direction Normal Stress, Pa" ',
!     &  '"Y-Direction Normal Stress, Pa" ',
!     &  '"Z-Direction Normal Stress, Pa" ',
!     &  '"X-Direction Shear Stress, Pa" ',
!     &  '"Y-Direction Shear Stress, Pa" ',
!     &  '"Z-Direction Shear Stress, Pa" '
!      NTX = 0
!      DO NFX = 1,NF_FRC
!        NTX = NTX + IP_FRC(2,NFX) - IP_FRC(1,NFX) + 1
!      ENDDO
!      WRITE(26,'(A,I9,A,I9,A)') 'ZONE T = "Fractures", NODES = ',
!     &  (3*NTX),', ELEMENTS = ',NTX,
!     &  ', DATAPACKING = BLOCK, ZONETYPE = FETRIANGLE'
!      WRITE(26,'(A)')'VARLOCATION=([4,5,6,7,8,9,10,11,12,13]' // 
!     &  '=CELLCENTERED)'
!!
!!---    Write x-dimension vertices  ---
!!
!      DO NFX = 1,NF_FRC
!        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!          WRITE(26,'(3(1PE16.9,1X))') ((VAR*XE_FRC(I,NTX)),I=1,3)
!        ENDDO
!      ENDDO
!!
!!---    Write y-dimension vertices  ---
!!
!      DO NFX = 1,NF_FRC
!        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!          WRITE(26,'(3(1PE16.9,1X))') ((VAR*YE_FRC(I,NTX)),I=1,3)
!        ENDDO
!      ENDDO
!!
!!---    Write z-dimension vertices  ---
!!
!      DO NFX = 1,NF_FRC
!        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!          WRITE(26,'(3(1PE16.9,1X))') ((VAR*ZE_FRC(I,NTX)),I=1,3)
!        ENDDO
!      ENDDO
!!
!!---  Write triangle areas  ---
!!
!      WRITE(26,'(10(1PE16.9,1X))') 
!     &  ((((VAR**2)*AF_FRC(NTX)),NTX=IP_FRC(1,NFX),IP_FRC(2,NFX)),
!     &  NFX=1,NF_FRC)
!!
!!---  Write triangle mechanical apertures  ---
!!
!      WRITE(26,'(10(1PE16.9,1X))') 
!     &  (((VAR*APM_FRC(2,NTX)),NTX=IP_FRC(1,NFX),IP_FRC(2,NFX)),
!     &  NFX=1,NF_FRC)
!!
!!---  Write triangle number  ---
!!
!      WRITE(26,'(10(1PE16.9,1X))') 
!     &  (((REAL(NTX)),NTX=IP_FRC(1,NFX),IP_FRC(2,NFX)),
!     &  NFX=1,NF_FRC)
!!
!!---  Write triangle normal stress, Pa  ---
!!
!      WRITE(26,'(10(1PE16.9,1X))') 
!     &  (((SQRT((SIGN_FRC(1,NTX)**2) + (SIGN_FRC(2,NTX)**2) +
!     &  (SIGN_FRC(3,NTX)**2))),NTX=IP_FRC(1,NFX),IP_FRC(2,NFX)),
!     &  NFX=1,NF_FRC)
!!
!!---  Write x-direction of the triangle normal stress, Pa  ---
!!
!      WRITE(26,'(10(1PE16.9,1X))') 
!     &  (((SIGN_FRC(1,NTX)),NTX=IP_FRC(1,NFX),IP_FRC(2,NFX)),
!     &  NFX=1,NF_FRC)
!!
!!---  Write y-direction of the triangle normal stress, Pa  ---
!!
!      WRITE(26,'(10(1PE16.9,1X))') 
!     &  (((SIGN_FRC(2,NTX)),NTX=IP_FRC(1,NFX),IP_FRC(2,NFX)),
!     &  NFX=1,NF_FRC)
!!
!!---  Write z-direction of the triangle normal stress, Pa  ---
!!
!      WRITE(26,'(10(1PE16.9,1X))') 
!     &  (((SIGN_FRC(3,NTX)),NTX=IP_FRC(1,NFX),IP_FRC(2,NFX)),
!     &  NFX=1,NF_FRC)
!!
!!---  Write x-direction of the triangle shear stress, Pa  ---
!!
!      WRITE(26,'(10(1PE16.9,1X))') 
!     &  (((SIGS_FRC(1,NTX)),NTX=IP_FRC(1,NFX),IP_FRC(2,NFX)),
!     &  NFX=1,NF_FRC)
!!
!!---  Write y-direction of the triangle shear stress, Pa  ---
!!
!      WRITE(26,'(10(1PE16.9,1X))') 
!     &  (((SIGS_FRC(2,NTX)),NTX=IP_FRC(1,NFX),IP_FRC(2,NFX)),
!     &  NFX=1,NF_FRC)
!!
!!---  Write z-direction of the triangle shear stress, Pa  ---
!!
!      WRITE(26,'(10(1PE16.9,1X))') 
!     &  (((SIGS_FRC(3,NTX)),NTX=IP_FRC(1,NFX),IP_FRC(2,NFX)),
!     &  NFX=1,NF_FRC)
!
!---  Surface normal x, y, z components  ---
!
      CALL SURFN_FRC
!
!---  Open 'connect_frc' file for writing  ---
!
      OPEN(UNIT=27, FILE='connect_frc', STATUS='UNKNOWN', 
     &  FORM='FORMATTED')
      CLOSE(UNIT=27,STATUS='DELETE')
      OPEN(UNIT=27, FILE='connect_frc', STATUS='NEW', 
     &  FORM='FORMATTED')
!
!---  Node connection map  ---
!
      NC1X = 0
      DO NF1X = 1,NF_FRC
      DO NT1X = IP_FRC(1,NF1X),IP_FRC(2,NF1X)
      IF( IXP_FRC(NT1X).EQ.0 ) THEN
        NC1X = NC1X + 3
        CYCLE
      ENDIF
      DO NV1X = 1,3
        NC1X = NC1X + 1
        NODX(NV1X) = NC1X
        NC2X = 0
        DO NF2X = 1,NF_FRC
        DO NT2X = IP_FRC(1,NF2X),IP_FRC(2,NF2X)
        DO NV2X = 1,3
          NC2X = NC2X + 1
          IF( NF1X.EQ.NF2X .AND. NT1X.EQ.NT2X .AND. NV1X.EQ.NV2X ) CYCLE
          DXE = ABS(XE_FRC(NV1X,NT1X)-XE_FRC(NV2X,NT2X))
          DYE = ABS(YE_FRC(NV1X,NT1X)-YE_FRC(NV2X,NT2X))
          DZE = ABS(ZE_FRC(NV1X,NT1X)-ZE_FRC(NV2X,NT2X))
          IF( (DXE+DYE+DZE)/EPSL.LT.EPSL ) THEN
            NODX(NV1X) = MIN( NODX(NV1X),NC2X )
          ENDIF
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      WRITE(27,'(3(I9,1X))') (NODX(I),I=1,3)
      ENDDO
      ENDDO      
!
!---  Dynamic memory deallocation  ---
!
      IF( ALLOCATED(XCTX) ) THEN
        DEALLOCATE( XCTX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: XCTX'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IF( ALLOCATED(YCTX) ) THEN
        DEALLOCATE( YCTX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: YCTX'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IF( ALLOCATED(ZCTX) ) THEN
        DEALLOCATE( ZCTX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: ZCTX'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Close fracture/fault data file  ---
!
      CLOSE(UNIT=25)
      CLOSE(UNIT=27)
      CLOSE(UNIT=36)
!
!---  Open fracture/fault preprocessed file for subsequent 
!     simulations  ---
!
      OPEN(UNIT=36, FILE='ppfrac.dat', STATUS='UNKNOWN', 
     &  FORM='FORMATTED')
      CLOSE(UNIT=36,STATUS='DELETE')
      OPEN(UNIT=36, FILE='ppfrac.dat', STATUS='NEW', 
     &  FORM='FORMATTED')
!
!---  Write fracture parameters  ---
!
      WRITE(36,*) LF_FRC,LFC_FRC,LNC_FRC,LT_FRC,LTC_FRC,LXP_FRC
!
!---  Write number of fractures and number of inactive
!     fracture triangles  ---
!
      WRITE(36,*) NF_FRC,NXP_FRC
!
!---  Loop over fractures/faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Write fracture/fault preprocessed file  ---
!
        CALL WRPREP_FRC( NFX )
      ENDDO     
!
!---  Close fracture/fault preprocessed file  ---
!
      CLOSE(UNIT=36)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDGEOM_FRC group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINAC_FRC( JNDX )
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
!     Read input file for inactive fracture triangle information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 20 February 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE GEOM_FRC
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINAC_FRC'
!
!---  Write card information to ouput file  ---
!
      IF( JNDX.EQ.1 ) THEN
        CARD = 'Inactive Fracture Triangle Card'
      ELSE
        CARD = 'Inactive Fault Triangle Card'
      ENDIF
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of lines of inactive triangle entries   ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Inactive Triangle Entries'
      CALL RDINT( ISTART,ICOMMA,CHDUM,NLINES )
      IF( JNDX.EQ.1 ) THEN
        WRITE(IWR,'(/,A)' ) 'Declared Inactive Fracture Triangles'
      ELSE
        WRITE(IWR,'(/,A)' ) 'Declared Inactive Fault Triangles'
      ENDIF
!
!---  Loop over number of inactive triangle entries  ---
!
      DO NL = 1,NLINES
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        IF( JNDX.EQ.1 ) THEN
          VARB = 'Fracture Number'
        ELSE
          VARB = 'Fault Number'
        ENDIF
        CALL RDINT( ISTART,ICOMMA,CHDUM,NFX )
        VARB = 'Triangle Number'
        CALL RDINT( ISTART,ICOMMA,CHDUM,NTLX )
!
!---    Check fracture and triangle numbers  ---
!
        IF( NTLX.LT.1 .OR. NTLX.GT.IP_FRC(2,NFX) ) THEN
          INDX = 7
          IF( JNDX.EQ.1 ) THEN
            CHMSG = 'Out-of-Range Fault Triangle'
          ELSE
            CHMSG = 'Out-of-Range Fracture Triangle'
          ENDIF
          IMSG = NTLX
          CALL WRMSGS( INDX )
        ELSE
          NTX = IP_FRC(1,NFX) + NTLX - 1
          IF( IXP_FRC(NTX).NE.0 ) THEN
            IXP_FRC(NTX) = 0
            NXP_FRC = NXP_FRC + 1
          ENDIF
!
!---      Record fracture/fault triangle as inactive  ---
!
          IF( JNDX.EQ.1 ) THEN
            WRITE(IWR,'(4X,A,I6,A,I6)') 'Fracture = ',NFX,
     &        ' Local Triangle = ',NTLX,' Global Triangle = ',NTX
          ELSE
            WRITE(IWR,'(4X,A,I6,A,I6)') 'Fault = ',NFX,
     &        ' Local Triangle = ',NTLX,' Global Triangle = ',NTX
          ENDIF
        ENDIF
      ENDDO
!
!---  Renumber fracture triangle equation indexing for
!     inactive fracture triangles  ---
!
      NC = 0
      DO NFX = 1,NF_FRC
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
          NC = NC + 1
          IXP_FRC(NTX) = NC
        ENDDO
      ENDDO
!
!---  Search for additional definitions  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINAC_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RD4850D_SIG
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
!     Read input file for stress around the 4850 West Access Drift.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 15 Nov. 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE GEO_MECH
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      CHARACTER*512 CHDUM
      REAL*8, DIMENSION(3,3) :: SIGX
      REAL*8, DIMENSION(3,3) :: V1X,V2X,V3X
      REAL*8, DIMENSION(4) :: SHMINX,SHMAXX,SVX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RD4850D_SIG'
!
!---  Write card information to ouput file  ---
!
      CARD = '4850 Drift Stress Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      WRITE (IWR,'(/,A)') 'Shmin Parameters'
!
!---  Read the Shmin Hill equation parameters  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Shmin Base'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SHMINX(1))
      VARB = 'Shmin Base Units'
      UNTS = 'pa'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IUNM = -1
      IUNKG = 1
      IUNS = -2
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',SHMINX(1)
      INDX = 0
      CALL RDUNIT(UNTS,SHMINX(1),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SHMINX(1),', Pa)'
      VARB = 'Shmin Max'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SHMINX(2))
      VARB = 'Shmin Max Units'
      UNTS = 'pa'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IUNM = -1
      IUNKG = 1
      IUNS = -2
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',SHMINX(2)
      INDX = 0
      CALL RDUNIT(UNTS,SHMINX(2),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SHMINX(2),', Pa)'
      VARB = 'Shmin Xhalf'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SHMINX(3))
      VARB = 'Shmin Xhalf Units'
      UNTS = 'm'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IUNM = 1
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',SHMINX(3)
      INDX = 0
      CALL RDUNIT(UNTS,SHMINX(3),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SHMINX(3),', m)'
      VARB = 'Shmin Rate'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SHMINX(4))
      WRITE(IWR,'(6X,2A,1PE11.4)') VARB(1:IVR),': ',SHMINX(4)
!
!---  Read the Shmax Hill equation parameters  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Shmax Base'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SHMAXX(1))
      VARB = 'Shmax Base Units'
      UNTS = 'pa'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IUNM = -1
      IUNKG = 1
      IUNS = -2
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',SHMAXX(1)
      INDX = 0
      CALL RDUNIT(UNTS,SHMAXX(1),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SHMAXX(1),', Pa)'
      VARB = 'Shmax Max'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SHMAXX(2))
      VARB = 'Shmax Max Units'
      UNTS = 'pa'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IUNM = -1
      IUNKG = 1
      IUNS = -2
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',SHMAXX(2)
      INDX = 0
      CALL RDUNIT(UNTS,SHMAXX(2),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SHMAXX(2),', Pa)'
      VARB = 'Shmax Xhalf'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SHMAXX(3))
      VARB = 'Shmax Xhalf Units'
      UNTS = 'm'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IUNM = 1
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',SHMAXX(3)
      INDX = 0
      CALL RDUNIT(UNTS,SHMAXX(3),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SHMAXX(3),', m)'
      VARB = 'Shmax Rate'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SHMAXX(4))
      WRITE(IWR,'(6X,2A,1PE11.4)') VARB(1:IVR),': ',SHMAXX(4)
!
!---  Read the Sv Hill equation parameters  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Sv Base'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SVX(1))
      VARB = 'Sv Base Units'
      UNTS = 'pa'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IUNM = -1
      IUNKG = 1
      IUNS = -2
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',SVX(1)
      INDX = 0
      CALL RDUNIT(UNTS,SVX(1),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SVX(1),', Pa)'
      VARB = 'Sv Max'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SVX(2))
      VARB = 'Sv Max Units'
      UNTS = 'pa'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IUNM = -1
      IUNKG = 1
      IUNS = -2
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',SVX(2)
      INDX = 0
      CALL RDUNIT(UNTS,SVX(2),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SVX(2),', Pa)'
      VARB = 'Sv Xhalf'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SVX(3))
      VARB = 'Sv Xhalf Units'
      UNTS = 'm'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IUNM = 1
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',SVX(3)
      INDX = 0
      CALL RDUNIT(UNTS,SVX(3),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SVX(3),', m)'
      VARB = 'Sv Rate'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SVX(4))
      WRITE(IWR,'(6X,2A,1PE11.4)') VARB(1:IVR),': ',SVX(4)
!
!---  Read the Euler angles  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Euler Angle Phi'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,PHIX)
      VARB = 'Euler Angle Phi Units'
      UNTS = 'deg'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',PHIX
      INDX = 0
      CALL RDUNIT(UNTS,PHIX,INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',PHIX,', rad)'
      VARB = 'Euler Angle Theta'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,THETAX)
      VARB = 'Euler Angle Theta Units'
      UNTS = 'deg'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',THETAX
      INDX = 0
      CALL RDUNIT(UNTS,THETAX,INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',THETAX,', rad)'
      VARB = 'Euler Angle Psi'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,PSIX)
      VARB = 'Euler Angle Psi Units'
      UNTS = 'deg'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      IDFLT = 1
      WRITE(IWR,'(6X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &  UNTS(1:NCH),': ',PSIX
      INDX = 0
      CALL RDUNIT(UNTS,PSIX,INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',PSIX,', rad)'
!
!---  Rotation matrix from Euler angles  ---
!
      V1X(1,1) = COS(PSIX)
      V1X(1,2) = SIN(PSIX)
      V1X(1,3) = 0.D+0
      V1X(2,1) = -SIN(PSIX)
      V1X(2,2) = COS(PSIX)
      V1X(2,3) = 0.D+0
      V1X(3,1) = 0.D+0
      V1X(3,2) = 1.D+0
      V2X(1,1) = 1.D+0
      V2X(1,2) = 0.D+0
      V2X(1,3) = 0.D+0
      V2X(2,1) = 0.D+0
      V2X(2,2) = COS(THETAX)
      V2X(2,3) = SIN(THETAX)
      V2X(3,1) = 0.D+0
      V2X(3,2) = -SIN(THETAX)
      V2X(3,3) = COS(THETAX)
      CALL MATMUL( V1X,V2X,V3X,3,3,3 )
      V2X(1,1) = COS(PHIX)
      V2X(1,2) = SIN(PHIX)
      V2X(1,3) = 0.D+0
      V2X(2,1) = -SIN(PHIX)
      V2X(2,2) = COS(PHIX)
      V2X(2,3) = 0.D+0
      V2X(3,1) = 0.D+0
      V2X(3,2) = 0.D+0
      V2X(3,3) = 1.D+0
      CALL MATMUL( V3X,V2X,V1X,3,3,3 )
      CALL MATTRP( V1X,V2X,3,3 )
!
!---  Loop over nodes, computing the principal stresses, based
!     on radial distance from the drift wall  ---
!
      DO N = 1,NFLD
        RADX = SQRT((YP(N)**2)+(ZP(N)**2)) - 1.8D+0
        SIGX(1,1) = SHMINX(1) + (SHMINX(2)-SHMINX(1))/(1.D+0 + 
     &    (SHMINX(3)/RADX)**SHMINX(4))
        SIGX(1,2) = 0.D+0
        SIGX(1,3) = 0.D+0
        SIGX(2,1) = 0.D+0
        SIGX(2,2) = SHMAXX(1) + (SHMAXX(2)-SHMAXX(1))/(1.D+0 + 
     &    (SHMAXX(3)/RADX)**SHMAXX(4))
        SIGX(2,3) = 0.D+0
        SIGX(3,1) = 0.D+0
        SIGX(3,2) = 0.D+0
        SIGX(3,3) = SVX(1) + (SVX(2)-SVX(1))/(1.D+0 + 
     &    (SVX(3)/RADX)**SVX(4))
!
!---    Rotate principal stresses, creating the stress tensor  ---
!
        CALL MATMUL( V1X,SIGX,V3X,3,3,3 )
        CALL MATMUL( V3X,V2X,SIGX,3,3,3 )
!
!---    Assign nodal stress tensor sigma-xx  ---
!
        SIG_GM(1,N) = SIGX(1,1)
!
!---    Assign nodal stress tensor sigma-yy  ---
!
        SIG_GM(2,N) = SIGX(2,2)
!
!---    Assign nodal stress tensor sigma-zz  ---
!
        SIG_GM(3,N) = SIGX(3,3)
!
!---    Assign nodal stress tensor sigma-yz  ---
!
        SIG_GM(4,N) = SIGX(2,3)
!
!---    Assign nodal stress tensor sigma-xz  ---
!
        SIG_GM(5,N) = SIGX(1,3)
!
!---    Assign nodal stress tensor sigma-xy  ---
!
        SIG_GM(6,N) = SIGX(1,2)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RD4850D_SIG group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT_VAR( AX,BX )
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
!     Sort two real variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 August 2017.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( AX.GT.BX ) THEN
        CX = AX
        AX = BX
        BX = CX
      ENDIF
!
!---  End of SORT_VAR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SURFN_FRC
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
!     Surface normals for fractures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 March 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GEOM_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER N1(3),N2(3)
      REAL*8 P1X(3),P2X(3),P3X(3)
!
!----------------------Data Statements---------------------------------!
!
      DATA N1 /1,2,3/
      DATA N2 /2,3,1/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SURFN_FRC'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fracture triangles  ---
!
        DO NT1X = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NT1X).EQ.0 ) CYCLE
          P1X(1) = XE_FRC(1,NT1X)
          P1X(2) = YE_FRC(1,NT1X)
          P1X(3) = ZE_FRC(1,NT1X)
          P2X(1) = XE_FRC(2,NT1X)
          P2X(2) = YE_FRC(2,NT1X)
          P2X(3) = ZE_FRC(2,NT1X)
          P3X(1) = XE_FRC(3,NT1X)
          P3X(2) = YE_FRC(3,NT1X)
          P3X(3) = ZE_FRC(3,NT1X)
          DO M = 1,3
            P2X(M) = P2X(M)-P1X(M)
            P3X(M) = P3X(M)-P1X(M)
          ENDDO
          CALL VCROSSP( P2X,P3X,P1X )
          SFNT_FRC(1,NT1X) = P1X(1)/SQRT(P1X(1)**2+P1X(2)**2+P1X(3)**2)
          SFNT_FRC(2,NT1X) = P1X(2)/SQRT(P1X(1)**2+P1X(2)**2+P1X(3)**2)
          SFNT_FRC(3,NT1X) = P1X(3)/SQRT(P1X(1)**2+P1X(2)**2+P1X(3)**2)
!
!---      Loop over fracture triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DXX = XP_FRC(NT2X)-XP_FRC(NT1X)
            DYX = YP_FRC(NT2X)-YP_FRC(NT1X)
            DZX = ZP_FRC(NT2X)-ZP_FRC(NT1X)
            DISTX = SQRT(DXX**2 + DYX**2 + DZX**2)
            SFNC_FRC(1,NCX) = DXX/DISTX
            SFNC_FRC(2,NCX) = DYX/DISTX
            SFNC_FRC(3,NCX) = DZX/DISTX
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SURFN_FRC group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PNT_TRI_CHK( VX,UX,I0X,I1X,ITX )
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
!     Check to see if triangle V is totally within triangle U
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2017.
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
      REAL*8 VX(3),UX(3,3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PNT_TRI_CHK'
      ITX = 0
      AX = UX(I1X,2) - UX(I1X,1)
      BX = -(UX(I0X,2) - UX(I0X,1))
      CX = -AX*UX(I0X,1) - BX*UX(I1X,1)
      D0X = AX*VX(I0X) + BX*VX(I1X) + CX
      AX = UX(I1X,3) - UX(I1X,2)
      BX = -(UX(I0X,3) - UX(I0X,2))
      CX = -AX*UX(I0X,2) - BX*UX(I1X,2)
      D1X = AX*VX(I0X) + BX*VX(I1X) + CX
      AX = UX(I1X,1) - UX(I1X,3)
      BX = -(UX(I0X,1) - UX(I0X,3))
      CX = -AX*UX(I0X,3) - BX*UX(I1X,3)
      D2X = AX*VX(I0X) + BX*VX(I1X) + CX
      IF( D0X*D1X .GT. 0.D+0 ) THEN
        IF( D0X*D2X .GT. 0.D+0 ) THEN
          ITX = 1
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PNT_TRI_CHK group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VSUB( V1X,V2X,V3X )
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
!     Subtraction of two three-dimensional vectors
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 August 2017.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 V1X(3),V2X(3),V3X(3)
!
!----------------------Executable Lines--------------------------------!
!
      V3X(1) = V1X(1) - V2X(1)
      V3X(2) = V1X(2) - V2X(2)
      V3X(3) = V1X(3) - V2X(3)
!
!---  End of VSUB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VSCALAR( V1X,VARX,V2X )
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
!     Multiple a vector by a scalar
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 March 2021.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 V1X(3),V2X(3)
!
!----------------------Executable Lines--------------------------------!
!
      V2X(1) = VARX*V1X(1)
      V2X(2) = VARX*V1X(2)
      V2X(3) = VARX*V1X(3)
!
!---  End of VSUB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VERTICES_BH( RBX,VX,XEX,YEX,ZEX )
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
!     Find the vertices of a borehole node.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 17 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XEX(4),YEX(4),ZEX(4)
      REAL*8 RM(3,3)
      REAL*8 UX(3),VX(3),WX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VERTICES_BH'
!
!---  Borehole generally headed in a global z-direction  ---
!
      IF( ABS(VX(3)).GE.ABS(VX(1)) .AND. 
     &  ABS(VX(3)).GE.ABS(VX(2)) ) THEN
        IF( VX(3).LE.0.D+0 ) THEN
          UX(1) = 0.D+0
          UX(2) = 0.D+0
          UX(3) = -1.D+0
        ELSE
          UX(1) = 0.D+0
          UX(2) = 0.D+0
          UX(3) = 1.D+0
        ENDIF
!
!---  Borehole generally headed in a global y-direction  ---
!
      ELSEIF( ABS(VX(2)).GE.ABS(VX(1)) .AND. 
     &  ABS(VX(2)).GE.ABS(VX(3)) ) THEN
        IF( VX(2).LE.0.D+0 ) THEN
          UX(1) = 0.D+0
          UX(2) = -1.D+0
          UX(3) = 0.D+0
        ELSE
          UX(1) = 0.D+0
          UX(2) = 1.D+0
          UX(3) = 0.D+0
        ENDIF
!
!---  Borehole generally headed in a global x-direction  ---
!
      ELSE
        IF( VX(1).LE.0.D+0 ) THEN
          UX(1) = -1.D+0
          UX(2) = 0.D+0
          UX(3) = 0.D+0
        ELSE
          UX(1) = 1.D+0
          UX(2) = 0.D+0
          UX(3) = 0.D+0
        ENDIF
      ENDIF
      CALL VCROSSP( UX,VX,WX )
      WNX = SQRT( (WX(1)**2) + WX(2)**2 + WX(3)**2 )
      IF( WNX.GT.EPSL ) THEN
        DO M = 1,3
          WX(M) = WX(M)/WNX
        ENDDO
      ELSE
        DO M = 1,3
          WX(M) = 0.D+0
        ENDDO
      ENDIF
      ALPHA = ACOS(VDOTP(UX,VX))
      RM(1,1) = WX(1)*WX(1)*(1.D+0-COS(ALPHA)) + COS(ALPHA)
      RM(1,2) = WX(1)*WX(2)*(1.D+0-COS(ALPHA)) - SIN(ALPHA)*WX(3)
      RM(1,3) = WX(1)*WX(3)*(1.D+0-COS(ALPHA)) + SIN(ALPHA)*WX(2)
      RM(2,1) = WX(1)*WX(2)*(1.D+0-COS(ALPHA)) + SIN(ALPHA)*WX(3)
      RM(2,2) = WX(2)*WX(2)*(1.D+0-COS(ALPHA)) + COS(ALPHA)
      RM(2,3) = WX(2)*WX(3)*(1.D+0-COS(ALPHA)) - SIN(ALPHA)*WX(1)
      RM(3,1) = WX(1)*WX(3)*(1.D+0-COS(ALPHA)) - SIN(ALPHA)*WX(2)
      RM(3,2) = WX(2)*WX(3)*(1.D+0-COS(ALPHA)) + SIN(ALPHA)*WX(1)
      RM(3,3) = WX(3)*WX(3)*(1.D+0-COS(ALPHA)) + COS(ALPHA)
!
!---  Borehole generally headed in a global z-direction  ---
!
      IF( ABS(VX(3)).GE.ABS(VX(1)) .AND. 
     &  ABS(VX(3)).GE.ABS(VX(2)) ) THEN
!
!---    Negative z-direction  ---
!
        IF( VX(3).LE.0.D+0 ) THEN
          UX(1) = -0.7071068*RBX
          UX(2) = -0.7071068*RBX
          UX(3) = 0.D+0
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(1) = VX(1)
          YEX(1) = VX(2)
          ZEX(1) = VX(3)
          UX(1) = 0.7071068*RBX
          UX(2) = -0.7071068*RBX
          UX(3) = 0.D+0
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(2) = VX(1)
          YEX(2) = VX(2)
          ZEX(2) = VX(3)
          UX(1) = -0.7071068*RBX
          UX(2) = 0.7071068*RBX
          UX(3) = 0.D+0
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(3) = VX(1)
          YEX(3) = VX(2)
          ZEX(3) = VX(3)
          UX(1) = 0.7071068*RBX
          UX(2) = 0.7071068*RBX
          UX(3) = 0.D+0
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(4) = VX(1)
          YEX(4) = VX(2)
          ZEX(4) = VX(3)
!
!---    Positive z-direction  ---
!
        ELSE
          UX(1) = -0.7071068*RBX
          UX(2) = -0.7071068*RBX
          UX(3) = 0.D+0
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(1) = VX(1)
          YEX(1) = VX(2)
          ZEX(1) = VX(3)
          UX(1) = -0.7071068*RBX
          UX(2) = 0.7071068*RBX
          UX(3) = 0.D+0
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(2) = VX(1)
          YEX(2) = VX(2)
          ZEX(2) = VX(3)
          UX(1) = 0.7071068*RBX
          UX(2) = -0.7071068*RBX
          UX(3) = 0.D+0
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(3) = VX(1)
          YEX(3) = VX(2)
          ZEX(3) = VX(3)
          UX(1) = 0.7071068*RBX
          UX(2) = 0.7071068*RBX
          UX(3) = 0.D+0
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(4) = VX(1)
          YEX(4) = VX(2)
          ZEX(4) = VX(3)
        ENDIF
!
!---  Borehole generally headed in a global y-direction  ---
!
      ELSEIF( ABS(VX(2)).GE.ABS(VX(1)) .AND. 
     &  ABS(VX(2)).GE.ABS(VX(3)) ) THEN
!
!---    Negative y-direction  ---
!
        IF( VX(2).LE.0.D+0 ) THEN
          UX(1) = -0.7071068*RBX
          UX(2) = 0.D+0
          UX(3) = -0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(1) = VX(1)
          YEX(1) = VX(2)
          ZEX(1) = VX(3)
          UX(1) = -0.7071068*RBX
          UX(2) = 0.D+0
          UX(3) = 0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(2) = VX(1)
          YEX(2) = VX(2)
          ZEX(2) = VX(3)
          UX(1) = 0.7071068*RBX
          UX(2) = 0.D+0
          UX(3) = -0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(3) = VX(1)
          YEX(3) = VX(2)
          ZEX(3) = VX(3)
          UX(1) = 0.7071068*RBX
          UX(2) = 0.D+0
          UX(3) = 0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(4) = VX(1)
          YEX(4) = VX(2)
          ZEX(4) = VX(3)
!
!---    Positive y-direction  ---
!
        ELSE
          UX(1) = -0.7071068*RBX
          UX(2) = 0.D+0
          UX(3) = -0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(1) = VX(1)
          YEX(1) = VX(2)
          ZEX(1) = VX(3)
          UX(1) = 0.7071068*RBX
          UX(2) = 0.D+0
          UX(3) = -0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(2) = VX(1)
          YEX(2) = VX(2)
          ZEX(2) = VX(3)
          UX(1) = -0.7071068*RBX
          UX(2) = 0.D+0
          UX(3) = 0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(3) = VX(1)
          YEX(3) = VX(2)
          ZEX(3) = VX(3)
          UX(1) = 0.7071068*RBX
          UX(2) = 0.D+0
          UX(3) = 0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(4) = VX(1)
          YEX(4) = VX(2)
          ZEX(4) = VX(3)
        ENDIF
!
!---  Borehole generally headed in a global x-direction  ---
!
      ELSE
!
!---    Negative x-direction  ---
!
        IF( VX(1).LE.0.D+0 ) THEN
          UX(1) = 0.D+0
          UX(2) = -0.7071068*RBX
          UX(3) = -0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(1) = VX(1)
          YEX(1) = VX(2)
          ZEX(1) = VX(3)
          UX(1) = 0.D+0
          UX(2) = 0.7071068*RBX
          UX(3) = -0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(2) = VX(1)
          YEX(2) = VX(2)
          ZEX(2) = VX(3)
          UX(1) = 0.D+0
          UX(2) = -0.7071068*RBX
          UX(3) = 0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(3) = VX(1)
          YEX(3) = VX(2)
          ZEX(3) = VX(3)
          UX(1) = 0.D+0
          UX(2) = 0.7071068*RBX
          UX(3) = 0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(4) = VX(1)
          YEX(4) = VX(2)
          ZEX(4) = VX(3)
!
!---    Positive x-direction  ---
!
        ELSE
          UX(1) = 0.D+0
          UX(2) = -0.7071068*RBX
          UX(3) = -0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(1) = VX(1)
          YEX(1) = VX(2)
          ZEX(1) = VX(3)
          UX(1) = 0.D+0
          UX(2) = -0.7071068*RBX
          UX(3) = 0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(2) = VX(1)
          YEX(2) = VX(2)
          ZEX(2) = VX(3)
          UX(1) = 0.D+0
          UX(2) = 0.7071068*RBX
          UX(3) = -0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(3) = VX(1)
          YEX(3) = VX(2)
          ZEX(3) = VX(3)
          UX(1) = 0.D+0
          UX(2) = 0.7071068*RBX
          UX(3) = 0.7071068*RBX
          CALL MATMUL( RM,UX,VX,3,3,1 )
          XEX(4) = VX(1)
          YEX(4) = VX(2)
          ZEX(4) = VX(3)
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VERTICES_BH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WR_BOREHOLE
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
!
!     Write borehole.dat file.
!     
!     ID_BH(1,LN_BH) - starting borehole interval index for borehole
!     ID_BH(2,LN_BH) - ending borehole interval index for borehole
!     ID_BH(3,LN_BH) - starting borehole node index for borehole
!     ID_BH(4,LN_BH) - ending borehole node index for borehole
!     ID_BH(5,LN_BH) - principal borehole node index for borehole
!
!     IT_BH(1,LN_BH) - type index for starting borehole node
!     IT_BH(2,LN_BH) - type index for ending borehole node
!     Borehole types
!     1 - gas mass injection w/ zero water
!     2 - gas mass injection w/ water mass frac.
!     3 - gas mass injection w/ water rel. humidity
!     4 - aqueous mass injection w/ zero air and salt
!     5 - aqueous mass injection w/ air and salt mass frac.
!     6 - aqueous mass injection w/ air and salt rel. sat.
!     7 - fixed pressure, state condition #1 (saturated)
!     8 - fixed pressure, state condition #2 (partially saturated)
!     9 - fixed pressure, state condition #3 (unsaturated)
!    11 - gas volumetric injection w/ zero water
!    12 - gas volumetric injection w/ water mass frac.
!    13 - gas volumetric injection w/ water rel. humidity
!    14 - aqueous volumetric injection w/ zero air and salt
!    15 - aqueous volumetric injection w/ air and salt mass frac.
!    16 - aqueous volumetric injection w/ air and salt rel. sat.
!    21 - Dirichlet energy
!    22 - Neumann energy
!
!     N_BH - number of coupled boreholes
!     XTP_BH(2,LI_BH) - x-transition points for borehole interval, m
!     YTP_BH(2,LI_BH) - y-transition points for borehole interval, m
!     ZTP_BH(2,LI_BH) - z-transition points for borehole interval, m
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 July 2014.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE OUTPU
      USE GEOM_BH
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*6 FORM1
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WR_BOREHOLE'
!
!---  Open 'borehole.dat' file for writing  ---
!
      OPEN(UNIT=26, FILE='borehole.dat', STATUS='UNKNOWN', 
     &  FORM='FORMATTED')
      CLOSE(UNIT=26,STATUS='DELETE')
      OPEN(UNIT=26, FILE='borehole.dat', STATUS='NEW', FORM='FORMATTED')
!
!---  Length conversion for output ---
!
      VARX = 1.D+0
      INDX = 4
      IUNM = 1
      CALL RDUNIT(UNLN,VARX,INDX)
!
!---  Loop over number of coupled boreholes  ---
!
      DO 400 NBH = 1,N_BH
!
!---    Write a purple borehole name text tag and geometry data to  
!       "borehole.dat" for an injection borehole  ---
!
        IF( IT_BH(1,NBH).GT.0 ) THEN
          NCH = INDEX( NM_BH(NBH),'  ')-1
          NICW = ID_BH(1,NBH)
          WRITE(26,'(2A,3(1PE12.5,A),2A)') 
     &      'TEXT C=PURPLE HU=POINT H=10 CS=GRID3D',
     &      ', X=',VARX*XTP_BH(1,NICW),
     &      ', Y=',VARX*YTP_BH(1,NICW),
     &      ', Z=',VARX*ZTP_BH(1,NICW),
     &      ', T="',NM_BH(NBH)(1:NCH),'"'
          WRITE(26,'(A)') 'GEOMETRY T=LINE3D, C=PURPLE LT=0.2'
          WRITE(26,'(A)') '1'
          NC = 0
          DO 110 NICW = ID_BH(1,NBH),ID_BH(2,NBH)
            IF( NICW.GT.ID_BH(1,NBH) ) THEN
              IF( ABS(XTP_BH(1,NICW)-XTP_BH(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(YTP_BH(1,NICW)-YTP_BH(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(ZTP_BH(1,NICW)-ZTP_BH(2,NICW-1)).LT.1.D-3 ) THEN
                NC = NC + 1
              ELSE
                NC = NC + 2
              ENDIF
            ELSE
              NC = NC + 1
            ENDIF
            IF( NICW.EQ.ID_BH(2,NBH) ) THEN
              NC = NC + 1
            ENDIF
  110     CONTINUE
          WRITE(FORM1(3:3),'(I1)') ICOUNT(NC)
          WRITE(26,FORM1) NC
          DO 112 NICW = ID_BH(1,NBH),ID_BH(2,NBH)
            IF( NICW.EQ.ID_BH(1,NBH) ) THEN
              WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_BH(1,NICW),
     &          VARX*YTP_BH(1,NICW),VARX*ZTP_BH(1,NICW)
            ENDIF
            IF( NICW.GT.ID_BH(1,NBH) ) THEN
              IF( ABS(XTP_BH(1,NICW)-XTP_BH(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(YTP_BH(1,NICW)-YTP_BH(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(ZTP_BH(1,NICW)-ZTP_BH(2,NICW-1)).LT.1.D-3 ) THEN
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_BH(1,NICW),
     &            VARX*YTP_BH(1,NICW),VARX*ZTP_BH(1,NICW)
              ELSE
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)') 
     &            VARX*XTP_BH(2,NICW-1),
     &            VARX*YTP_BH(2,NICW-1),VARX*ZTP_BH(2,NICW-1)
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_BH(1,NICW),
     &            VARX*YTP_BH(1,NICW),VARX*ZTP_BH(1,NICW)
              ENDIF
            ENDIF
            IF( NICW.EQ.ID_BH(2,NBH) ) THEN
              WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_BH(2,NICW),
     &          VARX*YTP_BH(2,NICW),VARX*ZTP_BH(2,NICW)
            ENDIF
  112     CONTINUE
!
!---    Write a cyan borehole name text tag and geometry data to  
!       "borehole.dat" for a production borehole  ---
!
        ELSEIF( IT_BH(1,NBH).LT.0 ) THEN
          NCH = INDEX( NM_BH(NBH),'  ')-1
          NICW = ID_BH(1,NBH)
          WRITE(26,'(2A,3(1PE12.5,A),2A)') 
     &      'TEXT C=CYAN HU=POINT H=10 CS=GRID3D',
     &      ', X=',VARX*XTP_BH(1,NICW),
     &      ', Y=',VARX*YTP_BH(1,NICW),
     &      ', Z=',VARX*ZTP_BH(1,NICW),
     &      ', T="',NM_BH(NBH)(1:NCH),'"'
          WRITE(26,'(A)')  'GEOMETRY T=LINE3D, C=CYAN LT=0.2'
          WRITE(26,'(A)') '1'
          NC = 0
          DO 120 NICW = ID_BH(1,NBH),ID_BH(2,NBH)
            IF( NICW.GT.ID_BH(1,NBH) ) THEN
              IF( ABS(XTP_BH(1,NICW)-XTP_BH(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(YTP_BH(1,NICW)-YTP_BH(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(ZTP_BH(1,NICW)-ZTP_BH(2,NICW-1)).LT.1.D-3 ) THEN
                NC = NC + 1
              ELSE
                NC = NC + 2
              ENDIF
            ELSE
              NC = NC + 1
            ENDIF
            IF( NICW.EQ.ID_BH(2,NBH) ) THEN
              NC = NC + 1
            ENDIF
  120     CONTINUE
          WRITE(FORM1(3:3),'(I1)') ICOUNT(NC)
          WRITE(26,FORM1) NC
          DO 122 NICW = ID_BH(1,NBH),ID_BH(2,NBH)
            IF( NICW.EQ.ID_BH(1,NBH) ) THEN
              WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_BH(1,NICW),
     &          VARX*YTP_BH(1,NICW),VARX*ZTP_BH(1,NICW)
            ENDIF
            IF( NICW.GT.ID_BH(1,NBH) ) THEN
              IF( ABS(XTP_BH(1,NICW)-XTP_BH(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(YTP_BH(1,NICW)-YTP_BH(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(ZTP_BH(1,NICW)-ZTP_BH(2,NICW-1)).LT.1.D-3 ) THEN
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_BH(1,NICW),
     &            VARX*YTP_BH(1,NICW),VARX*ZTP_BH(1,NICW)
              ELSE
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)')
     &            VARX*XTP_BH(2,NICW-1),
     &            VARX*YTP_BH(2,NICW-1),VARX*ZTP_BH(2,NICW-1)
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_BH(1,NICW),
     &            VARX*YTP_BH(1,NICW),VARX*ZTP_BH(1,NICW)
              ENDIF
            ENDIF
            IF( NICW.EQ.ID_BH(2,NBH) ) THEN
              WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_BH(2,NICW),
     &          VARX*YTP_BH(2,NICW),VARX*ZTP_BH(2,NICW)
            ENDIF
  122     CONTINUE
        ENDIF
  400 CONTINUE
!
!---  Close 'borehole.dat' file  ---
!
      CLOSE(UNIT=26)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WR_BOREHOLE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDPREP_FRC( NFX )
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
!     Read preprocessed input file for fractures/faults
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 9 February 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE OUTPU
      USE GRID
      USE GEOM_FRC
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
!
!----------------------Data Statements---------------------------------!
!
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDPREP_FRC'
!
!---  Read number of fracture triangles  ---
!
      READ(35,*) NTP_FRC(NFX),IP_FRC(1,NFX),IP_FRC(2,NFX)
      NT_FRC = NT_FRC + NTP_FRC(NFX)
      IF( NT_FRC.GT.LT_FRC ) THEN
        INDX = 5
        CHMSG = 'Number of Fracture/Fault Triangles > Parameter LT_FRC'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Read joint model parameters  ---
!
      READ(35,*) IJM_FRC(NFX),(DSBB_FRC(M,NFX),M=1,7),FRCX,BMAX,
     &  SKF_FRC(NFX),NM_FRC(NFX)
!
!---  Loop over the number of triangles in the fracture  ---
!
      DO MTX = 1,NTP_FRC(NFX)
        NTX = MTX + IP_FRC(1,NFX) - 1
!
!---    Fracture triangle index  ---
!
        READ(35,*) IXP_FRC(NTX)
!
!---    Read the fracture triangle vertices  ---
!
        READ(35,*) (XE_FRC(I,NTX),I=1,3),(YE_FRC(J,NTX),J=1,3),
     &    (ZE_FRC(K,NTX),K=1,3)
!
!---    Read the fracture triangle centroids and area  ---
!
        READ(35,*) XP_FRC(NTX),YP_FRC(NTX),ZP_FRC(NTX),AF_FRC(NTX)
!
!---    Surface normal vector of the fracture triangle ---
!
        READ(35,*) (SIGN_FRC(I,NTX),I=1,3),(SIGS_FRC(J,NTX),J=1,3),
     &    (SFNT_FRC(K,NTX),K=1,3)
!
!---    Read triangle joint parameters  ---
!
        READ(35,*) JRC_FRC(NTX),JCS_FRC(NTX),UCS_FRC(NTX)
!
!---    Skip for inactive triangles  ---
!
        IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---    Resident node number for fracture triangle vertices  ---
!
        READ(35,*) (ND_FRC(I,NTX),I=1,3)
!
!---    Number of node to fracture and fracture to fracture 
!       connections ---
!
        READ(35,*) (IPN_FRC(I,NTX),I=1,2),(IPF_FRC(J,NTX),J=1,2)
!
!---    Loop over fracture triangle to triangle connections ---
!
        DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
          READ(35,*) ITCM_FRC(NCX),AFF_FRC(NCX),DFF_FRC(NCX),
     &      DFFM_FRC(NCX),(SFNC_FRC(M,NCX),M=1,3)
        ENDDO
!
!---    Loop over fracture triangle to grid-cell connections ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
          READ(35,*) INCM_FRC(NCX),AFN_FRC(NCX),DFN_FRC(NCX),
     &      DFNM_FRC(NCX),XPFN_FRC(NCX),YPFN_FRC(NCX),ZPFN_FRC(NCX)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDPREP_FRC group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRPREP_FRC( NFX )
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
!     Write preprocessed input file for fractures/faults
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 9 February 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE OUTPU
      USE GRID
      USE GEOM_FRC
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
!
!----------------------Data Statements---------------------------------!
!
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRPREP_FRC'
!
!---  Write number of fracture triangles  ---
!
      WRITE(36,*) NTP_FRC(NFX),IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---  Write joint model parameters  ---
!
      WRITE(36,*) IJM_FRC(NFX),(DSBB_FRC(M,NFX),M=1,7),FRCX,BMAX,
     &  SKF_FRC(NFX),NM_FRC(NFX)
!
!---  Loop over the number of triangles in the fracture  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Fracture triangle index  ---
!
        WRITE(36,*) IXP_FRC(NTX)
!
!---    Write the fracture triangle vertices  ---
!
        WRITE(36,*) (XE_FRC(I,NTX),I=1,3),(YE_FRC(J,NTX),J=1,3),
     &    (ZE_FRC(K,NTX),K=1,3)
!
!---    Write the fracture triangle centroids and area  ---
!
        WRITE(36,*) XP_FRC(NTX),YP_FRC(NTX),ZP_FRC(NTX),AF_FRC(NTX)
!
!---    Surface normal vector of the fracture triangle ---
!
        WRITE(36,*) (SIGN_FRC(I,NTX),I=1,3),(SIGS_FRC(J,NTX),J=1,3),
     &    (SFNT_FRC(K,NTX),K=1,3)
!
!---    Write triangle joint parameters  ---
!
        WRITE(36,*) JRC_FRC(NTX),JCS_FRC(NTX),UCS_FRC(NTX)
!
!---    Skip for inactive triangles  ---
!
        IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---    Resident node number for fracture triangle vertices  ---
!
        WRITE(36,*) (ND_FRC(I,NTX),I=1,3)
!
!---    Number of node to fracture and fracture to fracture 
!       connections ---
!
        WRITE(36,*) (IPN_FRC(I,NTX),I=1,2),(IPF_FRC(J,NTX),J=1,2)
!
!---    Loop over fracture triangle to triangle connections ---
!
        DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
          WRITE(36,*) ITCM_FRC(NCX),AFF_FRC(NCX),DFF_FRC(NCX),
     &      DFFM_FRC(NCX),(SFNC_FRC(M,NCX),M=1,3)
        ENDDO
!
!---    Loop over fracture triangle to grid-cell connections ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
          WRITE(36,*) INCM_FRC(NCX),AFN_FRC(NCX),DFN_FRC(NCX),
     &      DFNM_FRC(NCX),XPFN_FRC(NCX),YPFN_FRC(NCX),ZPFN_FRC(NCX)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRPREP_FRC group ---
!
      RETURN
      END



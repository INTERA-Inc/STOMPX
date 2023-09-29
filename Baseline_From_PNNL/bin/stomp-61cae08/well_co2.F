!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_COUP_WELL
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
!     STOMP-CO2
!
!     Define well nodes, determine trajectory points, and 
!     check for well trajectories within node surface planes
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 31 March 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
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
      REAL*8 XPX(5),YPX(5),ZPX(5)
      REAL*8 XIX(2),YIX(2),ZIX(2)
      REAL*8 AJ(3,3),BJ(3)
      INTEGER MSX(4,6),IJ(3)
      INTEGER N1X(4),N2X(4)
      CHARACTER*9 FORM1
      CHARACTER*8 FORM2
      CHARACTER*23 FORM5
!
!----------------------Data Statements---------------------------------!
!
      DATA MSX / 1,2,4,3,1,5,6,2,1,3,7,5,2,6,8,4,3,4,8,7,5,7,8,6 /
      DATA N1X / 2,3,4,1 /
      DATA N2X / 1,2,3,4 /
      DATA FORM1 /'(2X,I1,$)'/
      DATA FORM2 /'(2X,I1)'/
      DATA FORM5 /'(A,I3,A,I3,A,I3,A,I8,A)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHK_COUP_WELL'
      EPSLX = 1.D-12
!
!---  Loop over coupled wells ---
!
      DO 600 NCW = 1,N_CW
        ID_CW(3,NCW) = NWN_CW+1
        ID_CW(5,NCW) = NWF_CW+1
!
!---    Loop over number of well intervals  ---
!
        DO 490 NICW = ID_CW(1,NCW),ID_CW(2,NCW)
!
!---    Loop over active nodes to find well nodes and well
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
                  DZPX1 = ZTP_CW(NPT,NICW)-ZPX(1)
                  DZPX2 = ZPX(2)-ZTP_CW(NPT,NICW)
                  IF( ABS(DZPX1).LT.EPSLX ) DZPX1 = 0.D+0
                  IF( ABS(DZPX2).LT.EPSLX ) DZPX2 = 0.D+0
!
!---              Transition point within vertical limits of node  ---
!
                  IF( DZPX1.GE.0.D+0 .AND. DZPX2.GE.0.D+0 ) THEN
                    NC = NC+1
                    XIX(NC) = 0.D+0
                    YIX(NC) = 0.D+0
                    ZIX(NC) = ZTP_CW(NPT,NICW)
                  ENDIF
                ENDIF
                GOTO 190
              ENDIF
!
!---          Check for point with hexahedron  ---
!
              CALL WITHIN( XTP_CW(NPT,NICW),YTP_CW(NPT,NICW),
     &          ZTP_CW(NPT,NICW),ICWX,N )
!
!---          Opposing rotations found, point outside hexahedron  ---
!
              IF( ICWX.EQ.0 ) GOTO 190
!
!---          No opposing rotations found, point inside hexahedron  ---
!
              NC = NC+1
              XIX(NC) = XTP_CW(NPT,NICW)
              YIX(NC) = YTP_CW(NPT,NICW)
              ZIX(NC) = ZTP_CW(NPT,NICW)
  190       CONTINUE
!
!---        Both transition points inside hexahedron, skip
!           search for well path crossing hexahedron surfaces  ---
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
              DZPX1 = ZPX(1)-ZTP_CW(1,NICW)
              IF( ABS(DZPX1).LT.EPSLX ) DZPX1 = 0.D+0
              DZPX2 = ZPX(1)-ZTP_CW(2,NICW)
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
              DZPX1 = ZPX(2)-ZTP_CW(1,NICW)
              IF( ABS(DZPX1).LT.EPSLX ) DZPX1 = 0.D+0
              DZPX2 = ZPX(2)-ZTP_CW(2,NICW)
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
                AJ(1,1) = XTP_CW(1,NICW)-XTP_CW(2,NICW)
                AJ(2,1) = YTP_CW(1,NICW)-YTP_CW(2,NICW)
                AJ(3,1) = ZTP_CW(1,NICW)-ZTP_CW(2,NICW)
                AJ(1,2) = XPX(N1X(NT))-XPX(5)
                AJ(2,2) = YPX(N1X(NT))-YPX(5)
                AJ(3,2) = ZPX(N1X(NT))-ZPX(5)
                AJ(1,3) = XPX(N2X(NT))-XPX(5)
                AJ(2,3) = YPX(N2X(NT))-YPX(5)
                AJ(3,3) = ZPX(N2X(NT))-ZPX(5)
                BJ(1) = XTP_CW(1,NICW)-XPX(5)
                BJ(2) = YTP_CW(1,NICW)-YPX(5)
                BJ(3) = ZTP_CW(1,NICW)-ZPX(5)
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
                  XTX = XTP_CW(1,NICW) 
     &              + (XTP_CW(2,NICW)-XTP_CW(1,NICW))*TX
                  YTX = YTP_CW(1,NICW)
     &              + (YTP_CW(2,NICW)-YTP_CW(1,NICW))*TX
                  ZTX = ZTP_CW(1,NICW)
     &              + (ZTP_CW(2,NICW)-ZTP_CW(1,NICW))*TX
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
!---        Check that two well points are distinct  ---
!
            IF( NC.EQ.2 ) THEN
              DPX = SQRT( ((XIX(2)-XIX(1))**2) + ((YIX(2)-YIX(1))**2)
     &          + ((ZIX(2)-ZIX(1))**2) )
            ENDIF
!
!---        Two distinct well points within node define or update
!           well node trajectory points  ---
!
            IF( NC.EQ.2 .AND. DPX.GT.EPSLX ) THEN
!
!---          Cylindrical coordinates with azimuthal symmetry  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) GOTO 362
!
!---          Check if line between well points is contained in
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
!---          Loop over well nodes in current well  ---
!
              DO 370 NWN = ID_CW(3,NCW),NWN_CW
!
!---            Well node previously counted  ---
!
                IF( IWN_CW(NWN).EQ.N ) THEN
                  GOTO 372
                ENDIF
  370         CONTINUE
!
!---          Increment the number of field nodes with
!             coupled-well nodes  ---
!
              NWF_CW = NWF_CW + 1
              IF( NWF_CW.GT.LWF_CW ) THEN
               INDX = 5
               CHMSG = 'Number of Field Nodes with Coupled Wells ' //
     &           '> Parameter LWF_CW'
               CALL WRMSGS( INDX )
              ENDIF
  372         CONTINUE
!
!---          Define a new well node and set the trajectory
!             points  ---
!
              NWN_CW = NWN_CW + 1
              IF( NWN_CW.GT.LWN_CW ) THEN
               INDX = 5
               CHMSG = 'Number of Coupled-Well Nodes ' //
     &           '> Parameter LWN_CW'
               CALL WRMSGS( INDX )
              ENDIF
              XP_CW(1,NWN_CW) = XIX(1)
              YP_CW(1,NWN_CW) = YIX(1)
              ZP_CW(1,NWN_CW) = ZIX(1)
              XP_CW(2,NWN_CW) = XIX(2)
              YP_CW(2,NWN_CW) = YIX(2)
              ZP_CW(2,NWN_CW) = ZIX(2)
              IWN_CW(NWN_CW) = N
              INV_CW(NWN_CW) = NICW
!
!---          Well projections on local grid  ---
!
              NWN = NWN_CW
              CALL PROJ_WELL( PLX_CW(NWN),PLY_CW(NWN),PLZ_CW(NWN),
     &          XP_CW(1,NWN),YP_CW(1,NWN),ZP_CW(1,NWN),N )
            ENDIF
  480     CONTINUE
  490   CONTINUE
        ID_CW(4,NCW) = NWN_CW
        ID_CW(6,NCW) = NWF_CW
!
!---    Sequence well nodes according to their distance from
!       the starting point of the well  ---
!
        I1X = ID_CW(1,NCW)
        I3X = ID_CW(3,NCW)
        I4X = ID_CW(4,NCW)
  500   CONTINUE
        DMNX = 1.D+20
        DO 510 KCW = I3X,I4X
          XWPX = 5.D-1*(XP_CW(1,KCW)+XP_CW(2,KCW))
          YWPX = 5.D-1*(YP_CW(1,KCW)+YP_CW(2,KCW))
          ZWPX = 5.D-1*(ZP_CW(1,KCW)+ZP_CW(2,KCW))
          DISTX = SQRT( (XTP_CW(1,I1X)-XWPX)**2 +
     &      (YTP_CW(1,I1X)-YWPX)**2 + (ZTP_CW(1,I1X)-ZWPX)**2 )
          IF( DISTX.LT.DMNX ) THEN
            DMNX = DISTX
            IMNX = KCW
            IWNX = IWN_CW(KCW)
            INVX = INV_CW(KCW)
            XIX(1) = XP_CW(1,KCW)
            XIX(2) = XP_CW(2,KCW)
            YIX(1) = YP_CW(1,KCW)
            YIX(2) = YP_CW(2,KCW)
            ZIX(1) = ZP_CW(1,KCW)
            ZIX(2) = ZP_CW(2,KCW)
            PLX_CWX = PLX_CW(KCW)
            PLY_CWX = PLY_CW(KCW)
            PLZ_CWX = PLZ_CW(KCW)
          ENDIF
  510   CONTINUE
        DO 520 JCW = I4X,I3X,-1
          IF( JCW.LT.IMNX ) THEN
            IWN_CW(JCW+1) = IWN_CW(JCW)
            INV_CW(JCW+1) = INV_CW(JCW)
            XP_CW(1,JCW+1) = XP_CW(1,JCW)
            XP_CW(2,JCW+1) = XP_CW(2,JCW)
            YP_CW(1,JCW+1) = YP_CW(1,JCW)
            YP_CW(2,JCW+1) = YP_CW(2,JCW)
            ZP_CW(1,JCW+1) = ZP_CW(1,JCW)
            ZP_CW(2,JCW+1) = ZP_CW(2,JCW)
            PLX_CW(JCW+1) = PLX_CW(JCW)
            PLY_CW(JCW+1) = PLY_CW(JCW)
            PLZ_CW(JCW+1) = PLZ_CW(JCW)
          ENDIF
  520   CONTINUE
        IWN_CW(I3X) = IWNX
        INV_CW(I3X) = INVX
        XP_CW(1,I3X) = XIX(1)
        XP_CW(2,I3X) = XIX(2)
        YP_CW(1,I3X) = YIX(1)
        YP_CW(2,I3X) = YIX(2)
        ZP_CW(1,I3X) = ZIX(1)
        ZP_CW(2,I3X) = ZIX(2)
        PLX_CW(I3X) = PLX_CWX
        PLY_CW(I3X) = PLY_CWX
        PLZ_CW(I3X) = PLZ_CWX
        I3X = I3X+1
        IF( I3X.LT.I4X ) GOTO 500
!
!---    Sequence well node points according to their distance from
!       the starting point of the well  ---
!
        I1X = ID_CW(1,NCW)
        I3X = ID_CW(3,NCW)
        I4X = ID_CW(4,NCW)
        DO 530 KCW = I3X,I4X
          DIST1X = SQRT( (XTP_CW(1,I1X)-XP_CW(1,KCW))**2 +
     &      (YTP_CW(1,I1X)-YP_CW(1,KCW))**2 + 
     &      (ZTP_CW(1,I1X)-ZP_CW(1,KCW))**2 )
          DIST2X = SQRT( (XTP_CW(1,I1X)-XP_CW(2,KCW))**2 +
     &      (YTP_CW(1,I1X)-YP_CW(2,KCW))**2 + 
     &      (ZTP_CW(1,I1X)-ZP_CW(2,KCW))**2 )
          IF( DIST1X.GT.DIST2X ) THEN
            XIX(2) = XP_CW(1,KCW)
            YIX(2) = YP_CW(1,KCW)
            ZIX(2) = ZP_CW(1,KCW)
            XP_CW(1,KCW) = XP_CW(2,KCW)
            YP_CW(1,KCW) = YP_CW(2,KCW)
            ZP_CW(1,KCW) = ZP_CW(2,KCW)
            XP_CW(2,KCW) = XIX(2)
            YP_CW(2,KCW) = YIX(2)
            ZP_CW(2,KCW) = ZIX(2)
          ENDIF
  530   CONTINUE
  600 CONTINUE
!
!---  Load array for field nodes that contain coupled-well nodes  ---
!
      NC = 0
      DO 640 NCW = 1,N_CW
        NUK_CW = (ISVC*(ABS(ID_CW(6,NCW)-ID_CW(5,NCW))+1))
        IF( NUK_CW.GT.LUK_CW ) THEN
          INDX = 5
          CHMSG = 'Number of Well Equation Unknowns ' //
     &      '> Parameter LUK_CW'
          CALL WRMSGS( INDX )
        ENDIF
        MC = 0
        DO 630 KCW = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(KCW)
          DO 610 JCW = ID_CW(5,NCW),ID_CW(5,NCW)-1+MC
            IF( IWF_CW(JCW).EQ.N ) THEN
              IWP_CW(KCW) = JCW
              GOTO 620
            ENDIF
  610     CONTINUE
          NC = NC + 1
          MC = MC + 1
          IWF_CW(NC) = N
          IWP_CW(KCW) = NC
  620     CONTINUE
  630   CONTINUE
  640 CONTINUE
!
!---  Initialize coupled well pressure to be in hydrostatic
!     equilibrium with local pressure  ---
!
      DO 700 NCW = 1,N_CW
!
!---    Coupled-well pressure previously initialized  ---
!
        IF( P_CW(2,NCW).GT.-(HDOD*RHORL*GRAV) ) GOTO 700
        CALL EQUIL_COUP_WELL( NCW )
  700 CONTINUE
!
!---  Record coupled-well nodes to output file ---
!
      WRITE(IWR,'(//,A,/)') ' --- Coupled-Well Node Record  ---'
      DO 720 NCW = 1,N_CW
        NCH = INDEX( WNM_CW(NCW),'  ') - 1
        WRITE(IWR,'(2X,A,I4,2A)') 'Nodes Containing Coupled-Well ' //
     &    'Number: ',NCW,' Coupled-Well Name: ',WNM_CW(NCW)(1:NCH)
        DO 710 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
          NX = IWF_CW(NWF)
          IX = ID(NX)
          JX = JD(NX)
          KX = KD(NX)
          WRITE(IWR,FORM5) ' (',IX,',',JX,',',KX,':',NX,')'
  710   CONTINUE
  720 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_LEAK_WELL
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
!     Define well nodes, determine trajectory points, and 
!     check for well trajectories within node surface planes
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 31 March 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE OUTPU
      USE LEAK_WELL
      USE GRID
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
      REAL*8 AJ(3,3),BJ(3)
      INTEGER MSX(4,6),IJ(3)
      INTEGER N1X(4),N2X(4)
      CHARACTER*9 FORM1
      CHARACTER*8 FORM2
      CHARACTER*23 FORM5
!
!----------------------Data Statements---------------------------------!
!
      DATA MSX / 1,2,4,3,1,5,6,2,1,3,7,5,2,6,8,4,3,4,8,7,5,7,8,6 /
      DATA N1X / 2,3,4,1 /
      DATA N2X / 1,2,3,4 /
      DATA FORM1 /'(2X,I1,$)'/
      DATA FORM2 /'(2X,I1)'/
      DATA FORM5 /'(A,I3,A,I3,A,I3,A,I8,A)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHK_LEAK_WELL'
      EPSLX = 1.D-12
!
!---  Loop over leaky wells ---
!
      DO 600 NLW = 1,N_LW
        ID_LW(3,NLW) = NWN_LW+1
        ID_LW(5,NLW) = NWF_LW+1
!
!---    Loop over number of well intervals  ---
!
        DO 490 NILW = ID_LW(1,NLW),ID_LW(2,NLW)
!
!---    Loop over active nodes to find well nodes and well
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
                  DZPX1 = ZTP_LW(NPT,NILW)-ZPX(1)
                  DZPX2 = ZPX(2)-ZTP_LW(NPT,NILW)
                  IF( ABS(DZPX1).LT.EPSLX ) DZPX1 = 0.D+0
                  IF( ABS(DZPX2).LT.EPSLX ) DZPX2 = 0.D+0
!
!---              Transition point within vertical limits of node  ---
!
                  IF( DZPX1.GE.0.D+0 .AND. DZPX2.GE.0.D+0 ) THEN
                    NC = NC+1
                    XIX(NC) = 0.D+0
                    YIX(NC) = 0.D+0
                    ZIX(NC) = ZTP_LW(NPT,NILW)
                  ENDIF
                ENDIF
                GOTO 190
              ENDIF
!
!---          Check for point with hexahedron  ---
!
              CALL WITHIN( XTP_LW(NPT,NILW),YTP_LW(NPT,NILW),
     &          ZTP_LW(NPT,NILW),ILWX,N )
!
!---          Opposing rotations found, point outside hexahedron  ---
!
              IF( ILWX.EQ.0 ) GOTO 190
!
!---          No opposing rotations found, point inside hexahedron  ---
!
              NC = NC+1
              XIX(NC) = XTP_LW(NPT,NILW)
              YIX(NC) = YTP_LW(NPT,NILW)
              ZIX(NC) = ZTP_LW(NPT,NILW)
  190       CONTINUE
!
!---        Both transition points inside hexahedron, skip
!           search for well path crossing hexahedron surfaces  ---
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
              DZPX1 = ZPX(1)-ZTP_LW(1,NILW)
              IF( ABS(DZPX1).LT.EPSLX ) DZPX1 = 0.D+0
              DZPX2 = ZPX(1)-ZTP_LW(2,NILW)
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
              DZPX1 = ZPX(2)-ZTP_LW(1,NILW)
              IF( ABS(DZPX1).LT.EPSLX ) DZPX1 = 0.D+0
              DZPX2 = ZPX(2)-ZTP_LW(2,NILW)
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
                AJ(1,1) = XTP_LW(1,NILW)-XTP_LW(2,NILW)
                AJ(2,1) = YTP_LW(1,NILW)-YTP_LW(2,NILW)
                AJ(3,1) = ZTP_LW(1,NILW)-ZTP_LW(2,NILW)
                AJ(1,2) = XPX(N1X(NT))-XPX(5)
                AJ(2,2) = YPX(N1X(NT))-YPX(5)
                AJ(3,2) = ZPX(N1X(NT))-ZPX(5)
                AJ(1,3) = XPX(N2X(NT))-XPX(5)
                AJ(2,3) = YPX(N2X(NT))-YPX(5)
                AJ(3,3) = ZPX(N2X(NT))-ZPX(5)
                BJ(1) = XTP_LW(1,NILW)-XPX(5)
                BJ(2) = YTP_LW(1,NILW)-YPX(5)
                BJ(3) = ZTP_LW(1,NILW)-ZPX(5)
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
                  XTX = XTP_LW(1,NILW) 
     &              + (XTP_LW(2,NILW)-XTP_LW(1,NILW))*TX
                  YTX = YTP_LW(1,NILW)
     &              + (YTP_LW(2,NILW)-YTP_LW(1,NILW))*TX
                  ZTX = ZTP_LW(1,NILW)
     &              + (ZTP_LW(2,NILW)-ZTP_LW(1,NILW))*TX
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
                    CHMSG = 'Three Distinct Leaky Well Points'
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
!---        Check that two well points are distinct  ---
!
            IF( NC.EQ.2 ) THEN
              DPX = SQRT( ((XIX(2)-XIX(1))**2) + ((YIX(2)-YIX(1))**2)
     &          + ((ZIX(2)-ZIX(1))**2) )
            ENDIF
!
!---        Two distinct well points within node define or update
!           well node trajectory points  ---
!
            IF( NC.EQ.2 .AND. DPX.GT.EPSLX ) THEN
!
!---          Cylindrical coordinates with azimuthal symmetry  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) GOTO 362
!
!---          Check if line between well points is contained in
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
!---          Loop over well nodes in current well  ---
!
              DO 370 NWN = ID_LW(3,NLW),NWN_LW
!
!---            Well node previously counted  ---
!
                IF( NF_LW(NWN).EQ.N ) THEN
                  GOTO 372
                ENDIF
  370         CONTINUE
!
!---          Increment the number of field nodes with
!             leaky well nodes  ---
!
              NWF_LW = NWF_LW + 1
              IF( NWF_LW.GT.LWF_LW ) THEN
               INDX = 5
               CHMSG = 'Number of Field Nodes with Leaky Wells ' //
     &           '> Parameter LWF_LW'
               CALL WRMSGS( INDX )
              ENDIF
  372         CONTINUE
!
!---          Define a new well node and set the trajectory
!             points  ---
!
              NWN_LW = NWN_LW + 1
              IF( NWN_LW.GT.LWN_LW ) THEN
               INDX = 5
               CHMSG = 'Number of Leaky Well Nodes ' //
     &           '> Parameter LWN_LW'
               CALL WRMSGS( INDX )
              ENDIF
              XP_LW(1,NWN_LW) = XIX(1)
              YP_LW(1,NWN_LW) = YIX(1)
              ZP_LW(1,NWN_LW) = ZIX(1)
              XP_LW(2,NWN_LW) = XIX(2)
              YP_LW(2,NWN_LW) = YIX(2)
              ZP_LW(2,NWN_LW) = ZIX(2)
              NF_LW(NWN_LW) = N
              INV_LW(NWN_LW) = NILW
!
!---          Well projections on local grid  ---
!
              NWN = NWN_LW
              CALL PROJ_WELL( PLX_LW(NWN),PLY_LW(NWN),PLZ_LW(NWN),
     &          XP_LW(1,NWN),YP_LW(1,NWN),ZP_LW(1,NWN),N )
            ENDIF
  480     CONTINUE
  490   CONTINUE
        ID_LW(4,NLW) = NWN_LW
        ID_LW(6,NLW) = NWF_LW
!
!---    Sequence well nodes according to their distance from
!       the starting point of the well  ---
!
        I1X = ID_LW(1,NLW)
        I3X = ID_LW(3,NLW)
        I4X = ID_LW(4,NLW)
  500   CONTINUE
        DMNX = 1.D+20
        DO 510 KLW = I3X,I4X
          XWPX = 5.D-1*(XP_LW(1,KLW)+XP_LW(2,KLW))
          YWPX = 5.D-1*(YP_LW(1,KLW)+YP_LW(2,KLW))
          ZWPX = 5.D-1*(ZP_LW(1,KLW)+ZP_LW(2,KLW))
          DISTX = SQRT( (XTP_LW(1,I1X)-XWPX)**2 +
     &      (YTP_LW(1,I1X)-YWPX)**2 + (ZTP_LW(1,I1X)-ZWPX)**2 )
          IF( DISTX.LT.DMNX ) THEN
            DMNX = DISTX
            IMNX = KLW
            IWNX = NF_LW(KLW)
            INVX = INV_LW(KLW)
            XIX(1) = XP_LW(1,KLW)
            XIX(2) = XP_LW(2,KLW)
            YIX(1) = YP_LW(1,KLW)
            YIX(2) = YP_LW(2,KLW)
            ZIX(1) = ZP_LW(1,KLW)
            ZIX(2) = ZP_LW(2,KLW)
            PLX_LWX = PLX_LW(KLW)
            PLY_LWX = PLY_LW(KLW)
            PLZ_LWX = PLZ_LW(KLW)
          ENDIF
  510   CONTINUE
        DO 520 JLW = I4X,I3X,-1
          IF( JLW.LT.IMNX ) THEN
            NF_LW(JLW+1) = NF_LW(JLW)
            INV_LW(JLW+1) = INV_LW(JLW)
            XP_LW(1,JLW+1) = XP_LW(1,JLW)
            XP_LW(2,JLW+1) = XP_LW(2,JLW)
            YP_LW(1,JLW+1) = YP_LW(1,JLW)
            YP_LW(2,JLW+1) = YP_LW(2,JLW)
            ZP_LW(1,JLW+1) = ZP_LW(1,JLW)
            ZP_LW(2,JLW+1) = ZP_LW(2,JLW)
            PLX_LW(JLW+1) = PLX_LW(JLW)
            PLY_LW(JLW+1) = PLY_LW(JLW)
            PLZ_LW(JLW+1) = PLZ_LW(JLW)
          ENDIF
  520   CONTINUE
        NF_LW(I3X) = IWNX
        INV_LW(I3X) = INVX
        XP_LW(1,I3X) = XIX(1)
        XP_LW(2,I3X) = XIX(2)
        YP_LW(1,I3X) = YIX(1)
        YP_LW(2,I3X) = YIX(2)
        ZP_LW(1,I3X) = ZIX(1)
        ZP_LW(2,I3X) = ZIX(2)
        PLX_LW(I3X) = PLX_LWX
        PLY_LW(I3X) = PLY_LWX
        PLZ_LW(I3X) = PLZ_LWX
        I3X = I3X+1
        IF( I3X.LT.I4X ) GOTO 500
!
!---    Sequence well node points according to their distance from
!       the starting point of the well  ---
!
        I1X = ID_LW(1,NLW)
        I3X = ID_LW(3,NLW)
        I4X = ID_LW(4,NLW)
        DO 530 KLW = I3X,I4X
          DIST1X = SQRT( (XTP_LW(1,I1X)-XP_LW(1,KLW))**2 +
     &      (YTP_LW(1,I1X)-YP_LW(1,KLW))**2 + 
     &      (ZTP_LW(1,I1X)-ZP_LW(1,KLW))**2 )
          DIST2X = SQRT( (XTP_LW(1,I1X)-XP_LW(2,KLW))**2 +
     &      (YTP_LW(1,I1X)-YP_LW(2,KLW))**2 + 
     &      (ZTP_LW(1,I1X)-ZP_LW(2,KLW))**2 )
          IF( DIST1X.GT.DIST2X ) THEN
            XIX(2) = XP_LW(1,KLW)
            YIX(2) = YP_LW(1,KLW)
            ZIX(2) = ZP_LW(1,KLW)
            XP_LW(1,KLW) = XP_LW(2,KLW)
            YP_LW(1,KLW) = YP_LW(2,KLW)
            ZP_LW(1,KLW) = ZP_LW(2,KLW)
            XP_LW(2,KLW) = XIX(2)
            YP_LW(2,KLW) = YIX(2)
            ZP_LW(2,KLW) = ZIX(2)
          ENDIF
  530   CONTINUE
  600 CONTINUE
!
!---  Load array for field nodes that contain leaky well nodes  ---
!
      NC = 0
      DO 640 NLW = 1,N_LW
        MC = 0
        DO 630 KLW = ID_LW(3,NLW),ID_LW(4,NLW)
          N = NF_LW(KLW)
          DO 610 JLW = ID_LW(5,NLW),ID_LW(5,NLW)-1+MC
            IF( IWF_LW(JLW).EQ.N ) THEN
              IWP_LW(KLW) = JLW
              GOTO 620
            ENDIF
  610     CONTINUE
          NC = NC + 1
          MC = MC + 1
          IWF_LW(NC) = N
          IWP_LW(KLW) = NC
  620     CONTINUE
  630   CONTINUE
  640 CONTINUE
!
!---  Record leaky well nodes to output file ---
!
      WRITE(IWR,'(//,A,/)') ' --- Leaky Well to Field Node ' //
     &  'Connection Record  ---'
      DO 720 NLW = 1,N_LW
        NCH = INDEX( WNM_LW(NLW),'  ') - 1
        WRITE(IWR,'(2X,A,I4,2A)') 'Field Nodes Connected to Leaky ' //
     &    'Well Number: ',NLW,' Leaky Well Name: ',WNM_LW(NLW)(1:NCH)
        DO 710 NWF = ID_LW(5,NLW),ID_LW(6,NLW)
          NX = IWF_LW(NWF)
          IX = ID(NX)
          JX = JD(NX)
          KX = KD(NX)
          WRITE(IWR,FORM5) ' (',IX,',',JX,',',KX,':',NX,')'
  710   CONTINUE
  720 CONTINUE
!
!---  Create intermediate leaky well nodes, loop over leaky wells ---
!
      DO NLW = 1,N_LW
!
!---    Loop over well nodes  ---
!
        DO KLW = ID_LW(3,NLW)+1,ID_LW(4,NLW)
!
!---    Identify discontinuous well node sections  ---
!
          DLX = SQRT( (XP_LW(1,KLW)-XP_LW(2,KLW-1))**2 + 
     &      (YP_LW(1,KLW)-YP_LW(2,KLW-1))**2 + 
     &      (ZP_LW(1,KLW)-ZP_LW(2,KLW-1))**2 )
          IF( DLX.GT.1.D-9 ) THEN
!
!---        Determine the well intervals of the two points  ---
!
            INV1X = 0
            INV2X = 0
            DO NILW = ID_LW(1,NLW),ID_LW(2,NLW)
              XMAX = MAX(XTP_LW(1,NILW),XTP_LW(2,NILW))
              XMIN = MIN(XTP_LW(1,NILW),XTP_LW(2,NILW))
              YMAX = MAX(YTP_LW(1,NILW),YTP_LW(2,NILW))
              YMIN = MIN(YTP_LW(1,NILW),YTP_LW(2,NILW))
              ZMAX = MAX(ZTP_LW(1,NILW),ZTP_LW(2,NILW))
              ZMIN = MIN(ZTP_LW(1,NILW),ZTP_LW(2,NILW))
              IF( XP_LW(1,KLW).GE.XMIN .AND. 
     &          XP_LW(1,KLW).LE.XMAX .AND.
     &          YP_LW(1,KLW).GE.YMIN .AND. 
     &          YP_LW(1,KLW).LE.YMAX .AND.
     &          ZP_LW(1,KLW).GE.ZMIN .AND. 
     &          ZP_LW(1,KLW).LE.ZMAX ) INV1X = NILW
              IF( XP_LW(2,KLW-1).GE.XMIN .AND. 
     &          XP_LW(2,KLW-1).LE.XMAX .AND.
     &          YP_LW(2,KLW-1).GE.YMIN .AND. 
     &          YP_LW(2,KLW-1).LE.YMAX .AND.
     &          ZP_LW(2,KLW-1).GE.ZMIN .AND. 
     &          ZP_LW(2,KLW-1).LE.ZMAX ) INV2X = NILW
            ENDDO
!
!---        Discontinuity found within an interval  ---
!
            IF( INV1X.EQ.INV2X ) THEN
              LWN_LW = LWN_LW + INT(DLX/PAR_LW(3,INV1X)) + 1
!
!---        Discontinuity between intervals  ---
!
            ELSE
              DL1X = SQRT( (XTP_LW(2,INV2X)-XP_LW(2,KLW-1))**2 + 
     &          (YTP_LW(2,INV2X)-YP_LW(2,KLW-1))**2 + 
     &          (ZTP_LW(2,INV2X)-ZP_LW(2,KLW-1))**2 )
              NWNX = 0
              IF( DL1X.GT.1.D-9 ) NWNX = INT(DL1X/PAR_LW(3,INV2X)) + 1
              DO NX = 1,NWNX
!
!---             Define a new well node and set the trajectory
!                points  ---
!
                 NWN_LW = NWN_LW + 1
                 IF( NWN_LW.GT.LWN_LW ) THEN
                   INDX = 5
                   CHMSG = 'Number of Leaky Well Nodes ' //
     &               '> Parameter LWN_LW'
                   CALL WRMSGS( INDX )
                  ENDIF
                  XP_LW(1,NWN_LW) = REAL(NX-1)*
     &              (XTP_LW(2,INV2X)-XP_LW(2,KLW-1))/
     &              REAL(NWNX) + XP_LW(2,KLW-1)
                  YP_LW(1,NWN_LW) = REAL(NX-1)*
     &              (YTP_LW(2,INV2X)-YP_LW(2,KLW-1))/
     &              REAL(NWNX) + YP_LW(2,KLW-1)
                  ZP_LW(1,NWN_LW) = REAL(NX-1)*
     &              (ZTP_LW(2,INV2X)-ZP_LW(2,KLW-1))/
     &              REAL(NWNX) + ZP_LW(2,KLW-1)
                  XP_LW(2,NWN_LW) = REAL(NX)*
     &              (XTP_LW(2,INV2X)-XP_LW(2,KLW-1))/
     &              REAL(NWNX) + XP_LW(2,KLW-1)
                  YP_LW(2,NWN_LW) = REAL(NX)*
     &              (YTP_LW(2,INV2X)-YP_LW(2,KLW-1))/
     &              REAL(NWNX) + YP_LW(2,KLW-1)
                  ZP_LW(2,NWN_LW) = REAL(NX)*
     &              (ZTP_LW(2,INV2X)-ZP_LW(2,KLW-1))/
     &              REAL(NWNX) + ZP_LW(2,KLW-1)
                  NF_LW(NWN_LW) = 0
                  INV_LW(NWN_LW) = INV2X
              ENDDO
              DL2X = SQRT( (XTP_LW(1,INV1X)-XP_LW(1,KLW))**2 + 
     &          (YTP_LW(1,INV1X)-YP_LW(1,KLW))**2 + 
     &          (ZTP_LW(1,INV1X)-ZP_LW(1,KLW))**2 )
              NWNX = 0
              IF( DL2X.GT.1.D-9 ) NWNX = INT(DL2X/PAR_LW(3,INV1X)) + 1
              DO NX = 1,NWNX
!
!---             Define a new well node and set the trajectory
!                points  ---
!
                 NWN_LW = NWN_LW + 1
                 IF( NWN_LW.GT.LWN_LW ) THEN
                   INDX = 5
                   CHMSG = 'Number of Leaky Well Nodes ' //
     &               '> Parameter LWN_LW'
                   CALL WRMSGS( INDX )
                  ENDIF
                  XP_LW(1,NWN_LW) = REAL(NX-1)*
     &              (XP_LW(1,KLW)-XTP_LW(1,INV1X))/
     &              REAL(NWNX) + XTP_LW(1,INV1X)
                  YP_LW(1,NWN_LW) = REAL(NX-1)*
     &              (YP_LW(1,KLW)-YTP_LW(1,INV1X))/
     &              REAL(NWNX) + YTP_LW(1,INV1X)
                  ZP_LW(1,NWN_LW) = REAL(NX-1)*
     &              (ZP_LW(1,KLW)-ZTP_LW(1,INV1X))/
     &              REAL(NWNX) + ZTP_LW(1,INV1X)
                  XP_LW(2,NWN_LW) = REAL(NX)*
     &              (XP_LW(1,KLW)-XTP_LW(1,INV1X))/
     &              REAL(NWNX) + XTP_LW(1,INV1X)
                  YP_LW(2,NWN_LW) = REAL(NX)*
     &              (YP_LW(1,KLW)-YTP_LW(1,INV1X))/
     &              REAL(NWNX) + YTP_LW(1,INV1X)
                  ZP_LW(2,NWN_LW) = REAL(NX)*
     &              (ZP_LW(1,KLW)-ZTP_LW(1,INV1X))/
     &              REAL(NWNX) + ZTP_LW(1,INV1X)
                  NF_LW(NWN_LW) = 0
                  INV_LW(NWN_LW) = INV1X
              ENDDO
            ENDIF
          ENDIF
        ENDDO
        ID_LW(4,NLW) = NWN_LW
        ID_LW(6,NLW) = NWF_LW
!
!---    Sequence well nodes according to their distance from
!       the starting point of the well  ---
!
        I1X = ID_LW(1,NLW)
        I3X = ID_LW(3,NLW)
        I4X = ID_LW(4,NLW)
  800   CONTINUE
        DMNX = 1.D+20
        DO 810 KLW = I3X,I4X
          XWPX = 5.D-1*(XP_LW(1,KLW)+XP_LW(2,KLW))
          YWPX = 5.D-1*(YP_LW(1,KLW)+YP_LW(2,KLW))
          ZWPX = 5.D-1*(ZP_LW(1,KLW)+ZP_LW(2,KLW))
          DISTX = SQRT( (XTP_LW(1,I1X)-XWPX)**2 +
     &      (YTP_LW(1,I1X)-YWPX)**2 + (ZTP_LW(1,I1X)-ZWPX)**2 )
          IF( DISTX.LT.DMNX ) THEN
            DMNX = DISTX
            IMNX = KLW
            IWNX = NF_LW(KLW)
            INVX = INV_LW(KLW)
            XIX(1) = XP_LW(1,KLW)
            XIX(2) = XP_LW(2,KLW)
            YIX(1) = YP_LW(1,KLW)
            YIX(2) = YP_LW(2,KLW)
            ZIX(1) = ZP_LW(1,KLW)
            ZIX(2) = ZP_LW(2,KLW)
            PLX_LWX = PLX_LW(KLW)
            PLY_LWX = PLY_LW(KLW)
            PLZ_LWX = PLZ_LW(KLW)
          ENDIF
  810   CONTINUE
        DO 820 JLW = I4X,I3X,-1
          IF( JLW.LT.IMNX ) THEN
            NF_LW(JLW+1) = NF_LW(JLW)
            INV_LW(JLW+1) = INV_LW(JLW)
            XP_LW(1,JLW+1) = XP_LW(1,JLW)
            XP_LW(2,JLW+1) = XP_LW(2,JLW)
            YP_LW(1,JLW+1) = YP_LW(1,JLW)
            YP_LW(2,JLW+1) = YP_LW(2,JLW)
            ZP_LW(1,JLW+1) = ZP_LW(1,JLW)
            ZP_LW(2,JLW+1) = ZP_LW(2,JLW)
            PLX_LW(JLW+1) = PLX_LW(JLW)
            PLY_LW(JLW+1) = PLY_LW(JLW)
            PLZ_LW(JLW+1) = PLZ_LW(JLW)
          ENDIF
  820   CONTINUE
        NF_LW(I3X) = IWNX
        INV_LW(I3X) = INVX
        XP_LW(1,I3X) = XIX(1)
        XP_LW(2,I3X) = XIX(2)
        YP_LW(1,I3X) = YIX(1)
        YP_LW(2,I3X) = YIX(2)
        ZP_LW(1,I3X) = ZIX(1)
        ZP_LW(2,I3X) = ZIX(2)
        PLX_LW(I3X) = PLX_LWX
        PLY_LW(I3X) = PLY_LWX
        PLZ_LW(I3X) = PLZ_LWX
        I3X = I3X+1
        IF( I3X.LT.I4X ) GOTO 800
!
!---    Sequence well node points according to their distance from
!       the starting point of the well  ---
!
        I1X = ID_LW(1,NLW)
        I3X = ID_LW(3,NLW)
        I4X = ID_LW(4,NLW)
        DO 830 KLW = I3X,I4X
          DIST1X = SQRT( (XTP_LW(1,I1X)-XP_LW(1,KLW))**2 +
     &      (YTP_LW(1,I1X)-YP_LW(1,KLW))**2 + 
     &      (ZTP_LW(1,I1X)-ZP_LW(1,KLW))**2 )
          DIST2X = SQRT( (XTP_LW(1,I1X)-XP_LW(2,KLW))**2 +
     &      (YTP_LW(1,I1X)-YP_LW(2,KLW))**2 + 
     &      (ZTP_LW(1,I1X)-ZP_LW(2,KLW))**2 )
          IF( DIST1X.GT.DIST2X ) THEN
            XIX(2) = XP_LW(1,KLW)
            YIX(2) = YP_LW(1,KLW)
            ZIX(2) = ZP_LW(1,KLW)
            XP_LW(1,KLW) = XP_LW(2,KLW)
            YP_LW(1,KLW) = YP_LW(2,KLW)
            ZP_LW(1,KLW) = ZP_LW(2,KLW)
            XP_LW(2,KLW) = XIX(2)
            YP_LW(2,KLW) = YIX(2)
            ZP_LW(2,KLW) = ZIX(2)
          ENDIF
  830   CONTINUE
      ENDDO
!
!---  Geometric factors for leaky well nodes, setting centroids  ---
!
      N = NFLD
!
!---  Loop over number of leaky wells  ---
!
      DO NLW = 1,N_LW
!
!---    Loop over number of leaky well nodes in the leaky well  ---
!
        DO NWN = ID_LW(3,NLW),ID_LW(4,NLW)
          N = N + 1
          ND_LW(NWN) = N
!
!---      Leaky well node centroid  ---
!
          XP(N) = 5.D-1*(XP_LW(2,NWN)+XP_LW(1,NWN))
          YP(N) = 5.D-1*(YP_LW(2,NWN)+YP_LW(1,NWN))
          ZP(N) = 5.D-1*(ZP_LW(2,NWN)+ZP_LW(1,NWN))
        ENDDO
      ENDDO
!
!---  Geometric factors for leaky well nodes, setting surfaces
!     and vertices  ---
!
      NPX = NSX(NFLD) + 1
      NPY = NSY(NFLD) + IFLD
      NPZ = NSZ(NFLD) + IJFLD
      NPZ = NPZ + 1
!
!---  Loop over number of leaky wells  ---
!
      DO NLW = 1,N_LW
!
!---    Loop over number of leaky well nodes in the leaky well  ---
!
        IFINDX = 0
        DO NWN = ID_LW(3,NLW),ID_LW(4,NLW)
          N = ND_LW(NWN)
          ID(N) = 0
          JD(N) = 0
          KD(N) = 0
          NPX = NPX + 1
          NPY = NPY + 1
          NPZ = NPZ + 1
          NSX(N) = NPX
          NSY(N) = NPY
          NSZ(N) = NPZ     
          NQX = NPX + 1
          NQY = NPY + 1
          INVX = INV_LW(NWN)
          IZ(N) = IZ_LW(INVX)
          IZN = IZ(N)
          PERMX = MAX( PERM(1,IZN),PERM(2,IZN),PERM(3,IZN) )
          PERM(1,IZN) = PERMX
          PERM(2,IZN) = PERMX
          PERM(3,IZN) = PERMX
          IXP(N) = N
!
!---      Connect leaky well nodes  ---
!
          IF( NF_LW(NWN).NE.0 .AND. IS_LW(INVX).EQ.1 ) IFINDX=1
          IF( NWN.GT.ID_LW(3,NLW) ) ICM(1,6,N) = N-1
          IF( NWN.LT.ID_LW(4,NLW) ) ICM(1,1,N) = N+1
!
!---      Block refinement indices  ---
!
          IBR(1,N) = 0
          IBR(2,N) = 0
          IBR(3,N) = 0
          IBR(4,N) = N
          IBR(5,N) = N
!
!---      Split west, east, north, and south surfaces on leaky
!         well nodes from domain  ---
!
          INBS(2,N) = NPY
          INBS(3,N) = NPX
          INBS(4,N) = NQX
          INBS(5,N) = NQY
!
!---      Length of leaky well node  ---
!
          DISTX = SQRT( (XP_LW(1,NWN)-XP_LW(2,NWN))**2 +
     &      (YP_LW(1,NWN)-YP_LW(2,NWN))**2 + 
     &      (ZP_LW(1,NWN)-ZP_LW(2,NWN))**2 )
!
!---      Leaky well vertices  ---
!
          RX = SQRT(GPI*(PAR_LW(2,INVX)**2))
          XE(1,N) = XP_LW(2,NWN) - 5.D-1*RX
          XE(2,N) = XP_LW(2,NWN) + 5.D-1*RX
          XE(3,N) = XP_LW(2,NWN) - 5.D-1*RX
          XE(4,N) = XP_LW(2,NWN) + 5.D-1*RX
          XE(5,N) = XP_LW(1,NWN) - 5.D-1*RX
          XE(6,N) = XP_LW(1,NWN) + 5.D-1*RX
          XE(7,N) = XP_LW(1,NWN) - 5.D-1*RX
          XE(8,N) = XP_LW(1,NWN) + 5.D-1*RX
          YE(1,N) = YP_LW(2,NWN) - 5.D-1*RX
          YE(2,N) = YP_LW(2,NWN) + 5.D-1*RX
          YE(3,N) = YP_LW(2,NWN) - 5.D-1*RX
          YE(4,N) = YP_LW(2,NWN) + 5.D-1*RX
          YE(5,N) = YP_LW(1,NWN) - 5.D-1*RX
          YE(6,N) = YP_LW(1,NWN) + 5.D-1*RX
          YE(7,N) = YP_LW(1,NWN) - 5.D-1*RX
          YE(8,N) = YP_LW(1,NWN) + 5.D-1*RX
          ZE(1,N) = ZP_LW(2,NWN)
          ZE(2,N) = ZP_LW(2,NWN)
          ZE(3,N) = ZP_LW(2,NWN)
          ZE(4,N) = ZP_LW(2,NWN)
          ZE(5,N) = ZP_LW(1,NWN)
          ZE(6,N) = ZP_LW(1,NWN)
          ZE(7,N) = ZP_LW(1,NWN)
          ZE(8,N) = ZP_LW(1,NWN)
!
!---      Leaky well volume and surface areas  ---
!
          VOL(N) = GPI*(PAR_LW(2,INVX)**2)*DISTX
          AFX(NPX) = GPI*(2.D+0*PAR_LW(2,INVX))*DISTX
          AFX(NQX) = GPI*(2.D+0*PAR_LW(2,INVX))*DISTX
          AFY(NPY) = GPI*(2.D+0*PAR_LW(2,INVX))*DISTX
          AFY(NQY) = GPI*(2.D+0*PAR_LW(2,INVX))*DISTX
          AFZ(NPZ) = GPI*(PAR_LW(2,INVX)**2)
          DXGF(N) = 2.D+0*PAR_LW(2,INVX)
          DYGF(N) = 2.D+0*PAR_LW(2,INVX)
          DZGF(N) = DISTX
          GRVX(NPX) = 0.D+0
          GRVX(NQX) = 0.D+0
          GRVY(NPY) = 0.D+0
          GRVY(NQY) = 0.D+0
          DXGP(NPX) = 5.D-1*DXGF(N)
          DXGP(NQX) = 5.D-1*DXGF(N)
          DYGP(NPY) = 5.D-1*DYGF(N)
          DYGP(NQY) = 5.D-1*DYGF(N)
          IF( NWN.EQ.ID_LW(4,NLW) ) THEN
            INBS(1,N) = NPZ
            DZGP(NPZ) = SQRT( (XP(N)-XP_LW(2,NWN))**2 +
     &        (YP(N)-YP_LW(2,NWN))**2 + 
     &        (ZP(N)-ZP_LW(2,NWN))**2 )
            VARX = SQRT((XP(N)-XP_LW(2,NWN))**2 + 
     &        (YP(N)-YP_LW(2,NWN))**2)
            IF( VARX.LT.EPSL ) THEN
              GRVZ(NPZ) = GRAV
            ELSE
              GRVZ(NPZ) = GRAV*SIN(ATAN((ZP(N)-ZP_LW(2,NWN))/VARX))
            ENDIF
          ELSE
            DZGP(NPZ) = SQRT( (XP(N)-XP(N+1))**2 + (YP(N)-YP(N+1))**2 + 
     &        (ZP(N)-ZP(N+1))**2 )
            VARX = SQRT((XP(N)-XP(N+1))**2 + (YP(N)-YP(N+1))**2)
            IF( VARX.LT.EPSL ) THEN
              GRVZ(NPZ) = GRAV
            ELSE
              GRVZ(NPZ) = GRAV*SIN(ATAN((ZP(N)-ZP(N+1))/VARX))
            ENDIF
          ENDIF
          IF( NWN.EQ.ID_LW(3,NLW) ) THEN
            NQZ = NPZ - 1
            INBS(6,N) = NQZ
            AFZ(NQZ) = GPI*(PAR_LW(2,INVX)**2)
            DZGP(NQZ) = SQRT( (XP_LW(1,NWN)-XP(N))**2 +
     &        (YP_LW(1,NWN)-YP(N))**2 + 
     &        (ZP_LW(1,NWN)-ZP(N))**2 )
            VARX = SQRT((XP_LW(1,NWN)-XP(N))**2 + 
     &        (YP_LW(1,NWN)-YP(N))**2)
            IF( VARX.LT.EPSL ) THEN
              GRVZ(NQZ) = GRAV
            ELSE
              GRVZ(NQZ) = GRAV*SIN(ATAN((ZP_LW(1,NWN)-ZP(N))/VARX))
            ENDIF
          ENDIF
          NPX = NPX + 1
          NPY = NPY + 1
        ENDDO
!
!---    No connection found between leaky well and field nodes  ---
!
        IF( IFINDX.EQ.0 ) THEN
          INDX = 7
          CARD = 'Leaky Well Card Card'
          IMSG = NLW
          CHMSG = 'No Connection Found Between Leaky ' // 
     &      'Well and Field Nodes: Leaky Well Number'
          CALL WRMSGS( INDX )
        ENDIF
      ENDDO
!
!---  Check for unrecognized leaky well reference nodes  ---
!
      DO N = 1,NREF
        IF( NDREF(N).LT.0 ) THEN
          NWN = -NDREF(N)
          IF( NWN.LT.1 .OR. NWN.GT.NWN_LW ) THEN
            CARD = 'Output Control Card'
            INDX = 7
            CHMSG = 'Unrecognized Leaky Well Reference Node Number'
            IMSG = NWN
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
      ENDDO
!
!---  Record leaky well nodes to output file ---
!
      WRITE(IWR,'(//,A,/)') ' --- Leaky Well Node Record  ---'
      WRITE(IWR,'(2X,A,I9)')'Total Number of Leaky Well Nodes = ',NWN_LW
!
!---  Loop over number of leaky wells  ---
!
      DO NLW = 1,N_LW
        NCH = INDEX( WNM_LW(NLW),'  ') - 1
        WRITE(IWR,'(2X,A,I4,2A)') 'Leaky Well Nodes for Leaky Well ' //
     &    'Number: ',NLW,' Leaky Well Name: ',WNM_LW(NLW)(1:NCH)
!
!---    Loop over number of leaky well nodes in the leaky well  ---
!
        DO NWN = ID_LW(3,NLW),ID_LW(4,NLW)
          N = ND_LW(NWN)
          WRITE(IWR,'(2X,A,I9,A,I9)') 'Leaky Well Node ',NWN,
     &      ' <=> Field Node ',N
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_LEAK_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE EQUIL_COUP_WELL( NCW )
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
!     STOMP-CO2
!
!     Equilibrate coupled-well pressure with formation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 May 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYST
      USE GRID
      USE FDVP
      USE COUP_WELL
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
      SUB_LOG(ISUB_LOG) = '/EQUIL_COUP_WELL'
!
!---  Injection well, equilibrate with first well node  ---
!
      IF( IT_CW(NCW).GT.0 ) THEN
        N = IWN_CW(ID_CW(3,NCW))
        XCWX = XTP_CW(1,ID_CW(1,NCW))
        YCWX = YTP_CW(1,ID_CW(1,NCW))
        ZCWX = ZTP_CW(1,ID_CW(1,NCW))
!
!---  Withdrawl well, equilibrate with last well node  ---
!
      ELSE
        N = IWN_CW(ID_CW(4,NCW))
        XCWX = XTP_CW(2,ID_CW(2,NCW))
        YCWX = YTP_CW(2,ID_CW(2,NCW))
        ZCWX = ZTP_CW(2,ID_CW(2,NCW))
      ENDIF
!
!---  Aqueous unsaturated conditions  ---
!
      IF( SG(2,N)-SGT(2,N).GT.EPSL ) THEN
        P_CW(2,NCW) = PG(2,N) - (ZCWX-ZP(N))*GRAV*RHOG(2,N)
!
!---  Aqueous saturated conditions  ---
!
      ELSE
        P_CW(2,NCW) = PL(2,N) - (ZCWX-ZP(N))*GRAV*RHOL(2,N)
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of EQUIL_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLUX_COUP_WELL
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
!     STOMP-CO2
!
!     Mass flux of water and CO2 between coupled-well nodes and 
!     field nodes.
!
!     CO2 mass balance residuals for injection type coupled wells.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 19 April 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE HYST
      USE GRID
      USE FDVP
      USE FDVG
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
      REAL*8 XPX(2),YPX(2),ZPX(2)
      REAL*8 VAR_CWX(6)
      INTEGER, DIMENSION(:), ALLOCATABLE :: MCW,MFD
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLUX_COUP_WELL'
!
!---  Dynamic memory allocation  ---
!
      ALLOCATE( MCW(1:(LUK+2)),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: MCW'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( MFD(1:(LUK+2)),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: MFD'
        CALL WRMSGS( INDX )
      ENDIF
      DO 10 M = 1,ISVC+2
        IF( M.NE.ISVC+2 ) THEN
          MCW(M) = 2
        ELSE
          MCW(M) = 3
        ENDIF
        IF( M.NE.ISVC+2 ) THEN
          MFD(M) = M+1
        ELSE
          MFD(M) = 2
        ENDIF
   10 CONTINUE        
!
!---  Initialize coupled-well parameters ---
!
      DO 40 NCW = 1,N_CW
!
!---    Zero coupled-well fluxes ---
!
        QM_CW(1,NCW) = 0.D+0
        QM_CW(3,NCW) = 0.D+0
        QM_CW(7,NCW) = 0.D+0
!
!---    Flow controlled well ---
!
        ID_CW(8,NCW) = 0
!
!---    Loop over coupled-well nodes  ---
!
        DO 30 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---      Zero volumetric injection well fluxes
!
!         Q_CW(1,NWN) - total volumetric flux, m^3/s
!         Q_CW(2,NWN) - aqueous volumetric flux, m^3/s
!         Q_CW(3,NWN) - gas volumetric flux, m^3/s
!
          DO M = 1,3
            Q_CW(M,NWN) = 0.D+0
          ENDDO
!
!---      Loop over increment indices  ---
!
          DO 20 M = 1,ISVC+2
            FXA_CW(M,NWN) = 0.D+0
            FXS_CW(M,NWN) = 0.D+0
            FXW_CW(M,NWN) = 0.D+0
   20     CONTINUE
   30  CONTINUE
   40 CONTINUE        
!
!---  Loop over coupled wells ---
!
      DO 500 NCW = 1,N_CW
        ICHK_CWX = 0
        DQ_CWX = 1.D-6
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        NCT = 0
        DO M = 1,IM_CW(NCW)
          NCT = NCT + IMP_CW(M,NCW)
        ENDDO
        IF( ICC_CW(NCW).EQ.1 ) TMZ = MOD( TM,VAR_CW(1,NCT,NCW) )
!
!---    Coupled well is inactive set well pressure to be in 
!       equilibrium with formation  ---
!
        IF( TMZ.LE.VAR_CW(1,1,NCW) .OR. 
     &    QM_CW(2,NCW).GE.TML_CW(NCW) ) THEN
          CALL EQUIL_COUP_WELL( NCW )
          ID_CW(8,NCW) = 1
          GOTO 500
        ENDIF
        IF( NCT.EQ.1 ) THEN
          DO 80 N = 2,6
            VAR_CWX(N) = VAR_CW(N,1,NCW)
   80     CONTINUE
!
!---      Limit injection rate by total injected mass  ---
!
          VAR_CWX(2) = MIN( VAR_CWX(2),
     &      ((TML_CW(NCW)-QM_CW(2,NCW))*DTI) )
!
!---      Set well state  ---
!
          IT_CWX = IT_CW(NCW)
          IF( IT_CW(NCW).EQ.100 ) IT_CWX = ITS_CW(1,NCW)
        ELSE
          DO 100 M = 2,NCT
            IF( TMZ.LE.VAR_CW(1,M,NCW) ) THEN
              TD_CW = VAR_CW(1,M,NCW)-VAR_CW(1,M-1,NCW)
              DT_CW = MIN( VAR_CW(1,M,NCW)-TMZ,DT )
              TF_CW = (TMZ-VAR_CW(1,M-1,NCW))/TD_CW
              DO 90 N = 2,6
                VAR_CWX(N) = VAR_CW(N,M-1,NCW) + 
     &            TF_CW*(VAR_CW(N,M,NCW)-VAR_CW(N,M-1,NCW))
   90         CONTINUE
!
!---          Limit injection rate by total injected mass  ---
!
              VAR_CWX(2) = MIN( VAR_CWX(2),
     &          ((TML_CW(NCW)-QM_CW(2,NCW))*DTI) )
!
!---          Set well state  ---
!
              IT_CWX = IT_CW(NCW)
              IF( IT_CW(NCW).EQ.100 ) THEN
                NC = 0
                DO N = 1,IM_CW(NCW)
                  NC = NC + IMP_CW(N,NCW)
                  IF( NC.GE.M ) THEN
                    IT_CWX = ITS_CW(N,NCW)
                    EXIT
                  ENDIF
                ENDDO
              ENDIF
              GOTO 110
            ENDIF
  100     CONTINUE
!
!---      Coupled well is inactive set well pressure to be in 
!         equilibrium with formation  ---
!
          CALL EQUIL_COUP_WELL( NCW )
          ID_CW(8,NCW) = 1
          GOTO 500
        ENDIF
  110   CONTINUE
!
!---    Load CO2 mass flux for use in RSDL_COUP_WELL ---
!
        FX_CW(NCW) = VAR_CWX(2)
!
!---    Load pressure limit for use in UPDT_COUP_WELL ---
!
        PL_CW(NCW) = VAR_CWX(3) - PATM
!
!---    Pressure controlled well ---
!
        IF( PL_CW(NCW)-P_CW(2,NCW).LT.EPSL ) THEN
          ID_CW(8,NCW) = 1
        ENDIF
!
!---    Excessive flow rate, pressure controlled well ---
!
        IF( VAR_CWX(2).GT.1.D+5 ) THEN
          ID_CW(8,NCW) = 1
          P_CW(2,NCW) = PL_CW(NCW)
        ENDIF
!
!---    Loop over increment indices ---
!
        DO 300 M = 1,ISVC+2
          MW = MCW(M)
          MF = MFD(M)
!
!---      Injection well  ---
!
          IF( IT_CW(NCW).GT.0 ) THEN
            N = IWN_CW(ID_CW(3,NCW))
            P_CWX = P_CW(MW,NCW)
            T_CWX = T(2,N)
!
!---        CO2 injection with zero solvated water  ---
!
            IF( IT_CWX.EQ.1 ) THEN
              PVW_CWX = 0.D+0
              PVA_CWX = MAX( P_CWX+PATM-PVW_CWX,0.D+0 )
              XMGA_CWX = 1.D+0
              XMGW_CWX = 0.D+0
              XGA_CWX = 1.D+0
              XGW_CWX = 0.D+0
              CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
              ISRX = 2
              CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
              RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---        CO2 injection with solvated water concentration 
!           declared as mass fraction  ---
!
            ELSEIF( IT_CWX.EQ.2 ) THEN
              XGW_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
              XGA_CWX = 1.D+0-XGW_CWX
              XMGA_CWX = (XGA_CWX/WTMA)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
              XMGW_CWX = (XGW_CWX/WTMW)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
              PVA_CWX = XMGA_CWX*(P_CWX+PATM)
              PVW_CWX = XMGW_CWX*(P_CWX+PATM)
              CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
              ISRX = 2
              CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
              RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---        CO2 injection with solvated water concentration 
!           declared as relative humidity  ---
!
            ELSEIF( IT_CWX.EQ.3 ) THEN
              CALL SP_W( T_CWX,PSW_CWX )
              PVW_CWX = PSW_CWX*VAR_CWX(4)
              PVA_CWX = MAX( P_CWX+PATM-PVW_CWX,0.D+0 )
              XMGA_CWX = PVA_CWX/(P_CWX+PATM)
              XMGW_CWX = PVW_CWX/(P_CWX+PATM)
              XGA_CWX = (XMGA_CWX*WTMA)/(XMGA_CWX*WTMA + XMGW_CWX*WTMW)
              XGW_CWX = (XMGW_CWX*WTMW)/(XMGA_CWX*WTMA + XMGW_CWX*WTMW)
              CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
              ISRX = 2
              CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
              RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---        Aqueous injection with zero dissolved CO2  ---
!
            ELSEIF( MOD(IT_CWX,10).EQ.4 ) THEN
              PL_CWX = MAX( P_CWX+PATM,0.D+0 )
              XLA_CWX = 0.D+0
!
!---          Aqueous injection with zero dissolved salt  ---
!
              IF( IT_CWX/10.EQ.1 ) THEN
                XLS_CWX = 0.D+0
!
!---          Aqueous injection with dissolved salt declared as 
!             brine mass fraction  ---
!
              ELSEIF( IT_CWX/10.EQ.2 ) THEN
                XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
!
!---          Aqueous injection with dissolved salt declared as 
!             brine relative saturation  ---
!
              ELSEIF( IT_CWX/10.EQ.3 ) THEN
                CALL SOL_LS( T_CWX,XLSMX )
                XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
              ENDIF
              CALL SP_B( T_CWX,XLS_CWX,PSBX )
              PVBX = PSBX
              CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
              PVWX = PVBX
              XLSSX = XLS_CWX
              CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &          XMGWX,XMLAX,XMLSX,XMLWX )
              XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
              XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
              CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
!
!---        Aqueous injection with dissolved CO2 declared as 
!           mass fraction  ---
!
            ELSEIF( MOD(IT_CWX,10).EQ.5 ) THEN
              PL_CWX = MAX( P_CWX+PATM,0.D+0 )
              XLA_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
!
!---          Aqueous injection with zero dissolved salt  ---
!
              IF( IT_CWX/10.EQ.1 ) THEN
                XLS_CWX = 0.D+0
!
!---          Aqueous injection with dissolved salt declared as 
!             brine mass fraction  ---
!
              ELSEIF( IT_CWX/10.EQ.2 ) THEN
                XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
!
!---          Aqueous injection with dissolved salt declared as 
!             brine relative saturation  ---
!
              ELSEIF( IT_CWX/10.EQ.3 ) THEN
                CALL SOL_LS( T_CWX,XLSMX )
                XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
              ENDIF
              CALL SP_B( T_CWX,XLS_CWX,PSBX )
              PVBX = PSBX
              CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
              PVWX = PVBX
              XLSSX = XLS_CWX
              CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &          XMGWX,XMLAX,XMLSX,XMLWX )
              XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
              XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
              CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
!
!---        Aqueous injection with dissolved CO2 declared as 
!           relative saturation  ---
!
            ELSEIF( MOD(IT_CWX,10).EQ.6 ) THEN
              PL_CWX = MAX( P_CWX+PATM,0.D+0 )
              XLA_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
!
!---          Aqueous injection with zero dissolved salt  ---
!
              IF( IT_CWX/10.EQ.1 ) THEN
                XLS_CWX = 0.D+0
!
!---          Aqueous injection with dissolved salt declared as 
!             brine mass fraction  ---
!
              ELSEIF( IT_CWX/10.EQ.2 ) THEN
                XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
 !
!---          Aqueous injection with dissolved salt declared as 
!             brine relative saturation  ---
!
              ELSEIF( IT_CWX/10.EQ.3 ) THEN
                CALL SOL_LS( T_CWX,XLSMX )
                XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
              ENDIF
              CALL SP_B( T_CWX,XLS_CWX,PSBX )
              PVBX = PSBX
              CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
              PVWX = PVBX
              XLSSX = XLS_CWX
              CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &          XMGWX,XMLAX,XMLSX,XMLWX )
              XLA_CWX = XLA_CWX*XLAX
              XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
              XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
              CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
            ENDIF
!
!---        Store top of coupled-well location in previous
!           coupled-well node location  ---
!
            XPX(1) = XTP_CW(1,ID_CW(1,NCW))
            YPX(1) = YTP_CW(1,ID_CW(1,NCW))
            ZPX(1) = ZTP_CW(1,ID_CW(1,NCW))
            ISX = ID_CW(3,NCW)
            IEX = ID_CW(4,NCW)
            IDX = 1
          ENDIF
!
!---      Loop over the nodes in the coupled well ---
!
          DO 200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
            N = IWN_CW(NWN)
            I = ID(N)
            INVX = INV_CW(NWN)
            IZN = IZ(N)
            T_CWX = T(2,N)
!
!---        Coupled-well node centroids and projections ---
!
!            XLX = ABS(XP_CW(2,NWN)-XP_CW(1,NWN))
!            YLX = ABS(YP_CW(2,NWN)-YP_CW(1,NWN))
!            ZLX = ABS(ZP_CW(2,NWN)-ZP_CW(1,NWN))
            XLX = PLX_CW(NWN)
            YLX = PLY_CW(NWN)
            ZLX = PLZ_CW(NWN)
            XPX(2) = 5.D-1*(XP_CW(2,NWN)+XP_CW(1,NWN))
            YPX(2) = 5.D-1*(YP_CW(2,NWN)+YP_CW(1,NWN))
            ZPX(2) = 5.D-1*(ZP_CW(2,NWN)+ZP_CW(1,NWN))
!
!---        Cylindrical coordinates with azimuthal symmetry,
!           centrally located wells  ---
!
            IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
              XPNX = 0.D+0
              YPNX = 0.D+0
              ZPNX = ZP(N)
!
!---        Cylindrical coordinates  ---
!
            ELSEIF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
              XPNX = XP(N)*COS(YP(N))
              YPNX = XP(N)*SIN(YP(N))
              ZPNX = ZP(N)
!
!---        Cartesian or boundary-fitted orthogonal coordinates  ---
!
            ELSE
              XPNX = XP(N)
              YPNX = YP(N)
              ZPNX = ZP(N)
            ENDIF
!
!---        Well pressure using previous coupled-well node density ---
!
            IF( IT_CWX.LE.3 ) THEN
              P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOG_CWX
            ELSE
              P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOL_CWX
            ENDIF
!
!---        Well pressure at the node centroid, used for coupled-well
!           nodal output  ---
!
            IF( M.EQ.1 ) THEN
              NWF = IWP_CW(NWN)
              IF( IT_CWX.LE.3 ) THEN
                PF_CW(NWF) = P_CWX - (ZPNX-ZPX(1))*GRAV*RHOG_CWX
              ELSE
                PF_CW(NWF) = P_CWX - (ZPNX-ZPX(1))*GRAV*RHOL_CWX
              ENDIF
            ENDIF
!
!---        Adjust the formation pressure to the coupled-well node
!           centroid  ---
!
            IF( (SG(MF,N)-SGT(MF,N)).GT.EPSL ) THEN
              PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(MF,N)
            ELSE
              PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            ENDIF
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
!
!---        Update coupled-well node density and viscosity
!           for CO2 injection, by using a fixed solvated
!           water mass fraction  ---
!
            IF( IT_CWX.LE.3 ) THEN
              XMGA_CWX = (XGA_CWX/WTMA)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
              XMGW_CWX = (XGW_CWX/WTMW)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
              PVA_CWX = XMGA_CWX*(P_CWX+PATM)
              PVW_CWX = XMGW_CWX*(P_CWX+PATM)
              CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
              ISRX = 2
              CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
              RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
              CALL VISC_A( T_CWX,RHOGA_CWX,VISGA_CWX )
              CALL VISC_W( T_CWX,PVW_CWX,RHOGW_CWX,VISGW_CWX )
              CALL VISC_G( VISGA_CWX,VISGW_CWX,XMGA_CWX,XMGW_CWX,
     &          VISG_CWX )
!
!---        Update coupled-well node density and viscosity
!           for aqueous injection, by using a fixed CO2 and salt
!           dissolved mass fraction  ---
!
            ELSE
              PL_CWX = MAX( P_CWX+PATM,0.D+0 )
              CALL SP_B( T_CWX,XLS_CWX,PSBX )
              PVAX = MAX( PL_CWX-PSBX,0.D+0 )
              CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
              WTMLX = 1.D+0/(XLA_CWX/WTMA+XLS_CWX/WTMS+XLW_CWX/WTMW)
              XMLA_CWX = XLA_CWX*WTMLX/WTMA
              CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
              ISRX = 1
              CALL DENS_W( T_CWX,PL_CWX,RHOLWX,RHOX,ISRX )
              CALL VISC_W( T_CWX,PL_CWX,RHOLWX,VISLWX )
              CALL VISC_B( T_CWX,XLS_CWX,VISLWX,VISBX )
              CALL DENS_A( T_CWX,PVAX,RHOGAX,I_VX )
              CALL VISC_A( T_CWX,RHOGAX,VISGAX )
              CALL VISC_L( XMLA_CWX,VISBX,VISGAX,VISL_CWX )
            ENDIF
!
!---        Equivalent field node radius components  ---
!
            PERMX = MAX( PERMV(1,N),1.D-20 )
            PERMY = MAX( PERMV(2,N),1.D-20 )
            PERMZ = MAX( PERMV(3,N),1.D-20 )
            RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---        Cylindrical coordinates with azimuthal symmetry,
!           centrally located wells  ---
!
            IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
              ROZ = RP(I)
              RWX = MIN( RWX,9.999D-1*ROZ )
              PERMX = PERMRF(MF,N)*PERMV(1,N)
              WI_CWX = 2.D+0*GPI*PERMX*ZLX/(LOG(ROZ/RWX)+PAR_CW(1,INVX))
            ELSE
              PERMYZ = SQRT(PERMY/PERMZ)
              PERMZY = SQRT(PERMZ/PERMY)
              DXGFX = DXGF(N)/FF_CW(1,NCW)
              DYGFX = DYGF(N)*RP(I)/FF_CW(2,NCW)
              DZGFX = DZGF(N)/FF_CW(3,NCW)
              ROX = 2.8D-1*SQRT(PERMYZ*(DZGFX**2) + PERMZY*(DYGFX**2))
     &        /(SQRT(PERMYZ)+SQRT(PERMZY))
              PERMZX = SQRT(PERMZ/PERMX)
              PERMXZ = SQRT(PERMX/PERMZ)
              ROY = 2.8D-1*SQRT(PERMZX*(DXGFX**2) + PERMXZ*(DZGFX**2))
     &          /(SQRT(PERMZX)+SQRT(PERMXZ))
              PERMYX = SQRT(PERMY/PERMX)
              PERMXY = SQRT(PERMX/PERMY)
              ROZ = 2.8D-1*SQRT(PERMYX*(DXGFX**2) + PERMXY*(DYGFX**2))
     &          /(SQRT(PERMYX)+SQRT(PERMXY))
!
!---          Well index components  ---
!
              PERMX = PERMRF(MF,N)*PERMV(1,N)
              PERMY = PERMRF(MF,N)*PERMV(2,N)
              PERMZ = PERMRF(MF,N)*PERMV(3,N)
              WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/
     &          (LOG(ROX/RWX)+PAR_CW(1,INVX))
              WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/
     &          (LOG(ROY/RWX)+PAR_CW(1,INVX))
              WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &          (LOG(ROZ/RWX)+PAR_CW(1,INVX))
              WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
            ENDIF
!
!---        Mass fluxes, positive into the node  ---
!
            DPGX = MAX( P_CWX-PGFX,0.D+0 )
            DPLX = MAX( P_CWX-MAX(PLFX,PGFX),0.D+0 )
!
!---        CO2 injection, flux from well to formation  ---
!
            IF( IT_CWX.LE.3 ) THEN
              FX_CWX = WI_CWX*RHOG_CWX*DPGX/VISG_CWX
              FXA_CW(M,NWN) = FX_CWX*XGA_CWX
              FXW_CW(M,NWN) = FX_CWX*XGW_CWX
!
!---          Volumetric injection well fluxes  ---
!
              Q_CW(1,NWN) = FX_CWX/RHOG_CWX
              Q_CW(3,NWN) = FX_CWX/RHOG_CWX
!
!---        Aqueous injection, flux from well to formation  ---
!
            ELSE
              FX_CWX = WI_CWX*RHOL_CWX*DPLX/VISL_CWX
              FXA_CW(M,NWN) = FX_CWX*XLA_CWX
              FXS_CW(M,NWN) = FX_CWX*XLS_CWX
              FXW_CW(M,NWN) = FX_CWX*XLW_CWX
!
!---          Volumetric injection well fluxes  ---
!
              Q_CW(1,NWN) = FX_CWX/RHOG_CWX
              Q_CW(2,NWN) = FX_CWX/RHOG_CWX
            ENDIF
!
!---        Store current coupled-well node location in previous
!           coupled-well node location  ---
!
            XPX(1) = XPX(2)
            YPX(1) = YPX(2)
            ZPX(1) = ZPX(2)
!            IF( M.EQ.1 .AND. FXA_CW(1,NWN).GT.EPSL ) 
!     &        PRINT *,'FXA_CW(1,',NWN,') = ',FXA_CW(1,NWN),
!     &        ' P_CWX = ',P_CWX,' PLFX = ',PLFX,
!     &        ' PGFX = ',PGFX,' WI_CWX = ',WI_CWX,
!     &        ' RHOG_CWX = ',RHOG_CWX,' DPGX = ',DPGX,
!     &        ' VISG_CWX = ',VISG_CWX
  200     CONTINUE
  300   CONTINUE
!
!---    CO2 mass balance residuals for injection type coupled well  ---
!
        NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
        MX = (NWFX*ISVC)+2
        RS_CW(1,NCW) = VAR_CWX(2)
        RS_CW(MX,NCW) = VAR_CWX(2)
        QM_CW(1,NCW) = 0.D+0
        QM_CW(3,NCW) = 0.D+0
        QM_CW(7,NCW) = 0.D+0
        QM_CWX = 0.D+0
!
!---    Loop over coupled-well nodes  ---
!
        DO 400 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          RS_CW(1,NCW) = RS_CW(1,NCW) - FXA_CW(1,NWN)
     &      - FXW_CW(1,NWN) - FXS_CW(1,NWN)
          RS_CW(MX,NCW) = RS_CW(MX,NCW) - FXA_CW(ISVC+2,NWN) 
     &      - FXW_CW(ISVC+2,NWN) - FXS_CW(ISVC+2,NWN)
          QM_CW(1,NCW) = QM_CW(1,NCW) + FXA_CW(1,NWN)
          QM_CW(3,NCW) = QM_CW(3,NCW) + FXW_CW(1,NWN)
          QM_CW(7,NCW) = QM_CW(7,NCW) + FXS_CW(1,NWN)
          QM_CWX = QM_CWX + FXA_CW(ISVC+2,NWN)
     &      + FXW_CW(ISVC+2,NWN) + FXS_CW(ISVC+2,NWN)
  400   CONTINUE
!!
!!---    Excessive increment in coupled-well pressure to create
!!       flow from well, decrease increment  ---
!!
!        IF( QM_CW(1,NCW).LT.EPSL .AND. QM_CWX.GT.VAR_CWX(2) .AND.
!     &    VAR_CWX(2).GT.EPSL .AND. ID_CW(8,NCW).EQ.0 ) THEN
!          DNR_CW(NCW) = 8.D-1*DNR_CW(NCW)
!          P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
!!
!!---      Zero coupled-well fluxes ---
!!
!          QM_CW(1,NCW) = 0.D+0
!          QM_CW(3,NCW) = 0.D+0
!!
!!---      Loop over coupled-well nodes  ---
!!
!          DO 412 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!!
!!---        Loop over increment indices  ---
!!
!            DO 410 M = 1,ISVC+2
!              FXA_CW(M,NWN) = 0.D+0
!              FXW_CW(M,NWN) = 0.D+0
!  410       CONTINUE
!  412     CONTINUE
!          GOTO 110
!
!---    Insufficient increment in coupled-well pressure to create
!       flow from well, increase increment and well pressure  ---
!
!        ELSEIF( QM_CWX.LT.EPSL ) THEN
!        VAR6X = MAX( EPSL,(1.D-6*VAR_CWX(2)) )
!        IF( (IT_CWX.LE.3 .AND. QM_CWX.LT.EPSL) .OR. 
!     &    (IT_CWX.GE.4 .AND. QM_CWX.LT.VAR6X) ) THEN
        IF( ABS(QM_CWX/DNR_CW(NCW)).LT.1.D-7 ) THEN
          DNR_CW(NCW) = 1.25D+0*DNR_CW(NCW)
          P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
          IF( P_CW(3,NCW).GE.PL_CW(NCW) ) THEN
!            IF( ICHK_CWX.EQ.0 ) THEN
!              P_CW(3,NCW) = PL_CW(NCW)
!              DNR_CW(NCW) = PL_CW(NCW) - P_CW(2,NCW)
!              ICHK_CWX = 1
!            ELSE
              P_CW(2,NCW) = PL_CW(NCW)
              DNR_CW(NCW) = 1.D-1
              P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
              ID_CW(8,NCW) = 1
              GOTO 500
 !           ENDIF
          ENDIF
!
!---      Zero coupled-well fluxes ---
!
          QM_CW(1,NCW) = 0.D+0
          QM_CW(3,NCW) = 0.D+0
          QM_CW(7,NCW) = 0.D+0
!
!---      Flow controlled well ---
!
          ID_CW(8,NCW) = 0
!
!---      Loop over coupled-well nodes  ---
!
          DO 422 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---        Loop over increment indices  ---
!
            DO 420 M = 1,ISVC+2
              FXA_CW(M,NWN) = 0.D+0
              FXS_CW(M,NWN) = 0.D+0
              FXW_CW(M,NWN) = 0.D+0
  420       CONTINUE
  422     CONTINUE
          GOTO 110
!
!---    Excessive well flow for pressure controlled well,
!       reduce well pressure  ---
!
        ELSEIF( ABS(ID_CW(8,NCW)).EQ.1 .AND. 
     &    (QM_CW(1,NCW)+QM_CW(3,NCW)+QM_CW(7,NCW)).GT.VAR_CWX(2) ) THEN
!          RS_CWX = VAR_CWX(2) - QM_CW(1,NCW)*(1.D+0-DQ_CWX)
!          DQ_CWX = 1.D+1*DQ_CWX
!          DP_CWX = -RS_CWX*DNR_CW(NCW)/
!     &      (RS_CW(MX,NCW)-RS_CW(1,NCW))
!          
!!          DP_CWX = -RS_CW(1,NCW)*DNR_CW(NCW)/
!!     &      (RS_CW(MX,NCW)-RS_CW(1,NCW))
          DQ_CWX = 1.D+1*DQ_CWX
          DP_CWX = -DQ_CWX*PATM
          P_CW(2,NCW) = P_CW(2,NCW) + DP_CWX
          P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
!
!---      Zero coupled-well fluxes ---
!
          QM_CW(1,NCW) = 0.D+0
          QM_CW(3,NCW) = 0.D+0
          QM_CW(7,NCW) = 0.D+0
!
!---      Pressure correction for pressure controlled well ---
!
          ID_CW(8,NCW) = -1
!
!---      Loop over coupled-well nodes  ---
!
          DO 432 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---        Loop over increment indices  ---
!
            DO 430 M = 1,ISVC+2
              FXA_CW(M,NWN) = 0.D+0
              FXS_CW(M,NWN) = 0.D+0
              FXW_CW(M,NWN) = 0.D+0
  430       CONTINUE
  432     CONTINUE
          GOTO 110
!
!---    Acceptable pressure reduction for pressure controlled well, 
!       switch to flow controlled well  ---
!
        ELSEIF( ID_CW(8,NCW).EQ.-1 .AND. 
     &    (QM_CW(1,NCW)+QM_CW(3,NCW)+QM_CW(7,NCW)).LT.VAR_CWX(2) ) THEN
          ID_CW(8,NCW) = 0
        ENDIF
!
!---    Loop over field nodes that contain coupled-well nodes  ---
!
        DO 470 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
          M1 = (NWF-ID_CW(5,NCW))*ISVC + 1
          DO 460 M2 = 1,ISVC
            M3 = M1+M2
            RS_CW(M3,NCW) = VAR_CWX(2)
!
!---        Loop over coupled-well nodes  ---
!
            DO 450 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---          If coupled-well node is within the current field
!             node, use incremented fluxes  ---
!
              IF( IWF_CW(NWF).EQ.IWN_CW(NWN) ) THEN
                RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXA_CW(M2+1,NWN)
     &            - FXW_CW(M2+1,NWN) - FXS_CW(M2+1,NWN)
!
!---          If coupled-well node is outside the current field
!             node, use un-incremented fluxes  ---
!
              ELSE
                RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXA_CW(1,NWN)
     &            - FXW_CW(1,NWN) - FXS_CW(1,NWN)
              ENDIF
  450       CONTINUE
  460     CONTINUE
  470   CONTINUE
  500 CONTINUE
!
!---  Dynamic memory deallocation  ---
!
      IF( ALLOCATED(MCW) ) THEN
      DEALLOCATE( MCW,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: MCW'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      IF( ALLOCATED(MFD) ) THEN
      DEALLOCATE( MFD,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: MFD'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLUX_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLUX_LEAK_WELL
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
!     STOMP-CO2
!
!     Aqueous, gas, and salt flux from leaky well nodes
!     to field nodes.
!
!     For increment indexing the leaky well node is considered to
!     be the lower node and the field node is considered to be the
!     upper node. Flux values are stored in the east surface locations
!     of the leaky well node (i.e., NQX = NSX(N) + 1)
!
!     UL(ISVF,NQX) - aqueous flux, m/s
!     UG(ISVF,NQX) - gas flux, m/s
!     UGW(ISVF,NQX) - water gas flux, kg/m^2 s
!     UDLA(ISVF,NQX) - CO2 aqueous flux, kmol/m^2 s
!     UDS(ISVF,NQX) - salt aqueous diffusive flux, kg/m^2 s
!     US(ISVF,NQX) - salt flux, kg/m^2 s
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 October 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE LEAK_WELL
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
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
      SUB_LOG(ISUB_LOG) = '/FLUX_LEAK_WELL'
!
!---  Initialize leaky well parameters ---
!
!
!---  Loop over number of leaky wells  ---
!
      DO NLW = 1,N_LW
!
!---    Loop over number of leaky well nodes in the leaky well  ---
!
        DO NWN = ID_LW(3,NLW),ID_LW(4,NLW)
          NW = ND_LW(NWN)
          N = NF_LW(NWN)
!
!---      Leaky well node is not connected to a field node  ---
!
          IF( N.EQ.0 ) CYCLE
          I = ID(N)
          INVX = INV_LW(NWN)
          NPX = NSX(NW)
          NQX = NSX(NW) + 1
!
!---      Equivalent field node radius components  ---
!
          PERMX = MAX( PERMV(1,N),1.D-20 )
          PERMY = MAX( PERMV(2,N),1.D-20 )
          PERMZ = MAX( PERMV(3,N),1.D-20 )
          RWX = MAX( PAR_LW(2,INVX),1.D-20 )
!
!---      Leaky well projections ---
!
          XLX = PLX_LW(NWN)
          YLX = PLY_LW(NWN)
          ZLX = PLZ_LW(NWN)
          DISTX = SQRT( (XLX**2) + (YLX**2) + (ZLX**2) )
          AREAX = 2.D+0*GPI*RWX*DISTX
!
!---      Aqueous and gas flux on screened intervals  ---
!
          IF( IS_LW(INVX).EQ.1 ) THEN
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
                ROZ = RP(I)
                RWX = MIN( RWX,9.999D-1*ROZ )
                PERMX = PERMRF(MP,N)*PERMV(1,N)
                WI_LWX = 2.D+0*GPI*PERMX*ZLX/
     &            (LOG(ROZ/RWX)+PAR_LW(1,INVX))
              ELSE
                PERMYZ = SQRT(PERMY/PERMZ)
                PERMZY = SQRT(PERMZ/PERMY)
                DXGFX = DXGF(N)/FF_LW(1,NLW)
                DYGFX = DYGF(N)*RP(I)/FF_LW(2,NLW)
                DZGFX = DZGF(N)/FF_LW(3,NLW)
                ROX = 2.8D-1*SQRT(PERMYZ*(DZGFX**2) + PERMZY*(DYGFX**2))
     &          /(SQRT(PERMYZ)+SQRT(PERMZY))
                PERMZX = SQRT(PERMZ/PERMX)
                PERMXZ = SQRT(PERMX/PERMZ)
                ROY = 2.8D-1*SQRT(PERMZX*(DXGFX**2) + PERMXZ*(DZGFX**2))
     &            /(SQRT(PERMZX)+SQRT(PERMXZ))
                PERMYX = SQRT(PERMY/PERMX)
                PERMXY = SQRT(PERMX/PERMY)
                ROZ = 2.8D-1*SQRT(PERMYX*(DXGFX**2) + PERMXY*(DYGFX**2))
     &            /(SQRT(PERMYX)+SQRT(PERMXY))
!
!---            Well index components  ---
!
                PERMX = PERMRF(MP,N)*PERMV(1,N)
                PERMY = PERMRF(MP,N)*PERMV(2,N)
                PERMZ = PERMRF(MP,N)*PERMV(3,N)
                WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/
     &            (LOG(ROX/RWX)+PAR_LW(1,INVX))
                WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/
     &            (LOG(ROY/RWX)+PAR_LW(1,INVX))
                WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &            (LOG(ROZ/RWX)+PAR_LW(1,INVX))
                WI_LWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
              ENDIF
!
!---          Aqueous flux  ---
!
              RHOLX = (RHOL(MN,NW)*VOL(NW) + RHOL(MP,N)*VOL(N))/
     &          (VOL(NW)+VOL(N))
              HDLX = PL(MN,NW) - PL(MP,N) + (ZP(NW)-ZP(N))*GRAV*RHOLX
              IF( M.EQ.1 ) HDL = HDLX
              PERMWX = PERMRF(MN,NW)*PERMV(3,NW)
              VARX = 2.D+0*GPI*DISTX*PERMWX
              IF( WI_LWX.GT.EPSL .AND. VARX.GT.EPSL ) THEN
                WIX = 1.D+0/((1.D+0/WI_LWX) + (1.D+0/VARX))
              ELSE
                WIX = 0.D+0
              ENDIF
              INDX = 8
              RKLM = DIFMN(RKL(1,MN,NW),RKL(1,MP,N),DXGF(NW),
     &          DXGF(N),HDL,INDX)
              INDX = 5
              VLM = DIFMN(VISL(MN,NW),VISL(MP,N),DXGF(NW),DXGF(N),
     &          HDL,INDX)
              UL(M,NQX) = WIX*RKLM*HDLX/VLM/AREAX
              UL(M,NPX) = UL(M,NQX)
!
!---          Gas flux  ---
!
              RHOGX = (RHOG(MN,NW)*VOL(NW) + RHOG(MP,N)*VOL(N))/
     &          (VOL(NW)+VOL(N))
              HDGX = PSO(MN,NW) - PSO(MP,N) + (ZP(NW)-ZP(N))*GRAV*RHOGX
              IF( M.EQ.1 ) HDG = HDGX
              PERMWX = PERMRF(MN,NW)*PERMV(3,NW)
              VARX = 2.D+0*GPI*DISTX*PERMWX
              IF( WI_LWX.GT.EPSL .AND. VARX.GT.EPSL ) THEN
                WIX = 1.D+0/((1.D+0/WI_LWX) + (1.D+0/VARX))
              ELSE
                WIX = 0.D+0
              ENDIF
              INDX = 9
              RKGM = DIFMN(RKG(MN,NW),RKG(MP,N),DXGF(NW),
     &          DXGF(N),HDG,INDX)
              INDX = 6
              VGM = DIFMN(VISG(MN,NW),VISG(MP,N),DXGF(NW),DXGF(N),
     &          HDG,INDX)
              UG(M,NQX) = WIX*RKGM*HDGX/VGM/AREAX
              UG(M,NPX) = UG(M,NQX)
            ENDDO
          ENDIF
!
!---      Water gas diffusion on screened intervals  ---
!
          IF( IS_LW(INVX).EQ.1 ) THEN
            DXMGW = XMGW(2,NW)-XMGW(2,N)
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &          DFGW(MP,N)
              DFW = TORG(MN,NW)*PORD(MN,NW)*(SG(MN,NW)-SGT(MN,NW))*
     &          DFGW(MN,NW)
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
                ROZ = RP(I)
                RWX = MIN( RWX,9.999D-1*ROZ )
                WI_LWX = 2.D+0*GPI*DFP*ZLX/
     &            (LOG(ROZ/RWX)+PAR_LW(1,INVX))
              ELSE
                DXGFX = DXGF(N)/FF_LW(1,NLW)
                DYGFX = DYGF(N)*RP(I)/FF_LW(2,NLW)
                DZGFX = DZGF(N)/FF_LW(3,NLW)
                ROX = 2.8D-1*SQRT((DZGFX**2) + (DYGFX**2))
                ROY = 2.8D-1*SQRT((DXGFX**2) + (DZGFX**2))
                ROZ = 2.8D-1*SQRT((DXGFX**2) + (DYGFX**2))
!
!---            Well index components  ---
!
                DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &            DFGW(MP,N)
                WIX = 2.D+0*GPI*DFP*XLX/(LOG(ROX/RWX)+PAR_LW(1,INVX))
                WIY = 2.D+0*GPI*DFP*YLX/(LOG(ROY/RWX)+PAR_LW(1,INVX))
                WIZ = 2.D+0*GPI*DFP*ZLX/(LOG(ROZ/RWX)+PAR_LW(1,INVX))
                WI_LWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
              ENDIF
              VARX = 2.D+0*GPI*DISTX*DFW
              IF( WI_LWX.GT.EPSL .AND. VARX.GT.EPSL ) THEN
                WIX = 1.D+0/((1.D+0/WI_LWX) + (1.D+0/VARX))
              ELSE
                WIX = 0.D+0
              ENDIF
              UDGW(M,NQX) = WIX*(XMGW(MN,NW)*RHOMG(MN,NW)
     &          - XMGW(MP,N)*RHOMG(MP,N))
              UDGW(M,NPX) = UDGW(M,NQX)
              FGWP = XGW(MP,N)*RHOG(MP,N)
              FGWW = XGW(MN,NW)*RHOG(MN,NW)
              INDX = 3
              FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,NQX),INDX )
              UGW(M,NQX) = UG(M,NQX)*FGW + WTMW*UDGW(M,NQX)
              UGW(M,NPX) = UGW(M,NQX)
            ENDDO
          ENDIF
!
!---      CO2 aqueous diffusion on screened intervals  ---
!
          IF( IS_LW(INVX).EQ.1 ) THEN
            DXLA = (XMLA(2,NW)-XMLA(2,N))
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &          *RHOML(MP,N)
              DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*DFLA(MN,NW)
     &          *RHOML(MN,NW)
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
                ROZ = RP(I)
                RWX = MIN( RWX,9.999D-1*ROZ )
                WI_LWX = 2.D+0*GPI*DFP*ZLX/
     &            (LOG(ROZ/RWX)+PAR_LW(1,INVX))
              ELSE
                DXGFX = DXGF(N)/FF_LW(1,NLW)
                DYGFX = DYGF(N)*RP(I)/FF_LW(2,NLW)
                DZGFX = DZGF(N)/FF_LW(3,NLW)
                ROX = 2.8D-1*SQRT((DZGFX**2) + (DYGFX**2))
                ROY = 2.8D-1*SQRT((DXGFX**2) + (DZGFX**2))
                ROZ = 2.8D-1*SQRT((DXGFX**2) + (DYGFX**2))
!
!---            Well index components  ---
!
                WIX = 2.D+0*GPI*DFP*XLX/(LOG(ROX/RWX)+PAR_LW(1,INVX))
                WIY = 2.D+0*GPI*DFP*YLX/(LOG(ROY/RWX)+PAR_LW(1,INVX))
                WIZ = 2.D+0*GPI*DFP*ZLX/(LOG(ROZ/RWX)+PAR_LW(1,INVX))
                WI_LWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
              ENDIF
              VARX = 2.D+0*GPI*DISTX*DFW
              IF( WI_LWX.GT.EPSL .AND. VARX.GT.EPSL ) THEN
                WIX = 1.D+0/((1.D+0/WI_LWX) + (1.D+0/VARX))
              ELSE
                WIX = 0.D+0
              ENDIF
              UDGW(M,NQX) = WIX*(XMLA(MN,NW)-XMLA(MP,N))
              UDGW(M,NPX) = UDGW(M,NQX)
            ENDDO
          ENDIF
!
!---      Salt aqueous flux on screened intervals  ---
!
          IF( IS_LW(INVX).EQ.1 ) THEN
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
!
!---          Diffusion coefficients  ---
!
              IF( IEDLS.EQ.1 ) THEN
                TCOR = (T(MP,N)+TABS)/TSPRF
                SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
                DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
                TCOR = (T(MN,NW)+TABS)/TSPRF
                SMDLP = DFLS(MN,NW)*TCOR*(VISRL/VISL(MN,NW))
                DFCLW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*SMDLP
              ELSEIF( IEDLS.EQ.2 ) THEN
                DFCLP = SDCLS(1,IZ(N))*SDCLS(2,IZ(N))*
     &            EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
                DFCLW = SDCLS(1,IZ(NW))*SDCLS(2,IZ(NW))*
     &            EXP(SL(MN,NW)*PORD(MN,NW)*SDCLS(3,IZ(NW)))
              ELSEIF( IEDLS.EQ.3 ) THEN
                DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
                DFCLW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*DFLS(MN,NW)
              ENDIF
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
                ROZ = RP(I)
                RWX = MIN( RWX,9.999D-1*ROZ )
                WI_LWX = 2.D+0*GPI*DFCLP*ZLX/
     &            (LOG(ROZ/RWX)+PAR_LW(1,INVX))
              ELSE
                DXGFX = DXGF(N)/FF_LW(1,NLW)
                DYGFX = DYGF(N)*RP(I)/FF_LW(2,NLW)
                DZGFX = DZGF(N)/FF_LW(3,NLW)
                ROX = 2.8D-1*SQRT((DZGFX**2) + (DYGFX**2))
                ROY = 2.8D-1*SQRT((DXGFX**2) + (DZGFX**2))
                ROZ = 2.8D-1*SQRT((DXGFX**2) + (DYGFX**2))
!
!---            Well index components  ---
!
                WIX = 2.D+0*GPI*DFCLP*XLX/(LOG(ROX/RWX)+PAR_LW(1,INVX))
                WIY = 2.D+0*GPI*DFCLP*YLX/(LOG(ROY/RWX)+PAR_LW(1,INVX))
                WIZ = 2.D+0*GPI*DFCLP*ZLX/(LOG(ROZ/RWX)+PAR_LW(1,INVX))
                WI_LWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
              ENDIF
              VARX = 2.D+0*GPI*DISTX*DFCLW
              IF( WI_LWX.GT.EPSL .AND. VARX.GT.EPSL ) THEN
                WIX = 1.D+0/((1.D+0/WI_LWX) + (1.D+0/VARX))
              ELSE
                WIX = 0.D+0
              ENDIF
              DDLW = WIX/AREAX
!
!---          Patankar salt transport  ---
!
              AL = MAX( UL(M,NQX),ZERO ) +
     &          DDLW*MAX((ONE-(TENTH*ABS(UL(M,NQX))/
     &          (DDLW+SMALL)))**5,ZERO)
              ALP = MAX( -UL(M,NQX),ZERO ) +
     &          DDLW*MAX((ONE-(TENTH*ABS(UL(M,NQX))/
     &          (DDLW+SMALL)))**5,ZERO)
              US(M,NQX) = XLS(MN,NW)*RHOL(MN,NW)*AL -
     &          XLS(MP,N)*RHOL(MP,N)*ALP
              US(M,NPX) = US(M,NQX)
              UDS(M,NQX) = DDLW*(XLS(MN,NW)*RHOL(MN,NW) -
     &          XLS(MP,N)*RHOL(MP,N))
              UDS(M,NPX) = UDS(M,NQX)     
            ENDDO
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLUX_LEAK_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_COUP_WELL
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
!     STOMP-CO2
!
!     Define well nodes, determine trajectory points, and 
!     check for well trajectories within node surface planes
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 14 April 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM_COUP_WELL'
!
!---  Loop over coupled wells ---
!
      DO 100 NCW = 1,N_CW
!
!---    Coupled-well is an injection type well, estimate the
!       well pressure and whether the well is pressure or flow
!       controlled  ---
!
        IF( IT_CW(NCW).GT.0 .AND. ISLC(70).EQ.1 ) THEN
          CALL INJP_COUP_WELL( NCW )
        ENDIF
        DNR_CW(NCW) = 1.D-1
        P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
  100 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INJ_COUP_WELL( NCW )
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
!     STOMP-CO2
!
!     Injection coupled well model
!     
!     Rate controlled or pressure controlled
!
!     Flux of water mass, and CO2 mass
!     from coupled-well nodes to field nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 March 2016.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE HYST
      USE GRID
      USE FDVP
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
      REAL*8 XPX(2),YPX(2),ZPX(2)
      REAL*8 VAR_CWX(6)
      INTEGER, DIMENSION(:), ALLOCATABLE :: MCW,MFD
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INJ_COUP_WELL'
!
!---  Dynamic memory allocation  ---
!
      ALLOCATE( MCW(1:(LUK+2)),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: MCW'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( MFD(1:(LUK+2)),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: MFD'
        CALL WRMSGS( INDX )
      ENDIF
      DO 10 M = 1,ISVC+2
        IF( M.NE.ISVC+2 ) THEN
          MCW(M) = 2
        ELSE
          MCW(M) = 3
        ENDIF
        IF( M.NE.ISVC+2 ) THEN
          MFD(M) = M+1
        ELSE
          MFD(M) = 2
        ENDIF
   10 CONTINUE        
!
!---  Zero injection well fluxes ---
!
      QM_CW(1,NCW) = 0.D+0
      QM_CW(3,NCW) = 0.D+0
      QM_CW(7,NCW) = 0.D+0
!
!---  Loop over coupled-well nodes  ---
!
      DO 30 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---    Zero volumetric injection well fluxes
!
!       Q_CW(1,NWN) - total volumetric flux, m^3/s
!       Q_CW(2,NWN) - aqueous volumetric flux, m^3/s
!       Q_CW(3,NWN) - gas volumetric flux, m^3/s
!
        DO M = 1,3
          Q_CW(M,NWN) = 0.D+0
        ENDDO
!
!---    Loop over increment indices  ---
!
        DO 20 M = 1,ISVC+2
          FXA_CW(M,NWN) = 0.D+0
          FXS_CW(M,NWN) = 0.D+0
          FXW_CW(M,NWN) = 0.D+0
   20   CONTINUE
   30 CONTINUE
!
!---  Injection well time interval ---
!
      DQ_CWX = 1.D-6
      TMZ = TM
      IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
      NCT = 0
      DO M = 1,IM_CW(NCW)
        NCT = NCT + IMP_CW(M,NCW)
      ENDDO
      IF( ICC_CW(NCW).EQ.1 ) TMZ = MOD( TM,VAR_CW(1,NCT,NCW) )
!
!---  Coupled well is inactive set well pressure to be in 
!     equilibrium with formation  ---
!
      IF( TMZ.LE.VAR_CW(1,1,NCW) .OR. 
     &  QM_CW(2,NCW).GE.TML_CW(NCW) ) THEN
        CALL EQUIL_COUP_WELL( NCW )
        ID_CW(8,NCW) = 1
        GOTO 500
      ENDIF
      IF( NCT.EQ.1 ) THEN
        DO 80 N = 2,6
          VAR_CWX(N) = VAR_CW(N,1,NCW)
   80   CONTINUE
!
!---    Limit injection rate by total injected mass  ---
!
        VAR_CWX(2) = MIN( VAR_CWX(2),
     &    ((TML_CW(NCW)-QM_CW(2,NCW))*DTI) )
!
!---    Set well state  ---
!
        IT_CWX = IT_CW(NCW)
        IF( IT_CW(NCW).EQ.100 ) IT_CWX = ITS_CW(1,NCW)
      ELSE
        DO 100 M = 2,NCT
          IF( TMZ.LE.VAR_CW(1,M,NCW) ) THEN
            TD_CW = VAR_CW(1,M,NCW)-VAR_CW(1,M-1,NCW)
            DT_CW = MIN( VAR_CW(1,M,NCW)-TMZ,DT )
            TF_CW = (TMZ-VAR_CW(1,M-1,NCW))/TD_CW
            DO 90 N = 2,6
              VAR_CWX(N) = VAR_CW(N,M-1,NCW) + 
     &          TF_CW*(VAR_CW(N,M,NCW)-VAR_CW(N,M-1,NCW))
   90       CONTINUE
!
!---        Limit injection rate by total injected mass  ---
!
            VAR_CWX(2) = MIN( VAR_CWX(2),
     &        ((TML_CW(NCW)-QM_CW(2,NCW))*DTI) )
!
!---        Set well state  ---
!
            IT_CWX = IT_CW(NCW)
            IF( IT_CW(NCW).EQ.100 ) THEN
              NC = 0
              DO N = 1,IM_CW(NCW)
                NC = NC + IMP_CW(N,NCW)
                IF( NC.GE.M ) THEN
                  IT_CWX = ITS_CW(N,NCW)
                  EXIT
                ENDIF
              ENDDO
            ENDIF
            GOTO 122
          ENDIF
  100   CONTINUE
      ENDIF
!
!---  Injection well is inactive set well pressure to be in 
!     equilibrium with reservoir  ---
!
      CALL EQUIL_COUP_WELL( NCW )
      ID_CW(8,NCW) = 1
      GOTO 500
  122 CONTINUE
!
!---  Load pressure limit ---
!
      PL_CW(NCW) = VAR_CWX(3) - PATM
!
!---  Loop over increment indices ---
!
      DO 300 M = 1,ISVC+2
        MW = MCW(M)
        MF = MFD(M)
        N = IWN_CW(ID_CW(3,NCW))
        P_CWX = P_CW(MW,NCW)
        T_CWX = T(2,N)
        XLS_CWX = 0.D+0
        CALL SP_W( T_CWX,PSW_CWX )
!
!---    Zero solvated water  ---
!
        IF( IT_CWX.EQ.1 ) THEN
          PVW_CWX = 0.D+0
          PVA_CWX = MAX( P_CWX+PATM-PVW_CWX,0.D+0 )
          XMGA_CWX = 1.D+0
          XMGW_CWX = 0.D+0
          XGA_CWX = 1.D+0
          XGW_CWX = 0.D+0
          CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
          ISRX = 2
          CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
          RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---    Solvated water concentration declared as mass fraction  ---
!
        ELSEIF( IT_CWX.EQ.2 ) THEN
          XGW_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
          XGA_CWX = 1.D+0-XGW_CWX
          XMGA_CWX = (XGA_CWX/WTMA)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
          XMGW_CWX = (XGW_CWX/WTMW)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
          PVA_CWX = XMGA_CWX*(P_CWX+PATM)
          PVW_CWX = XMGW_CWX*(P_CWX+PATM)
          CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
          ISRX = 2
          CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
          RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---    Solvated water concentration declared as relative
!       humidity---
!
        ELSEIF( IT_CWX.EQ.3 ) THEN
          PVW_CWX = PSW_CWX*VAR_CWX(4)
          PVA_CWX = MAX( P_CWX+PATM-PVW_CWX,0.D+0 )
          XMGA_CWX = PVA_CWX/(P_CWX+PATM)
          XMGW_CWX = PVW_CWX/(P_CWX+PATM)
          XGA_CWX = (XMGA_CWX*WTMA)/(XMGA_CWX*WTMA + XMGW_CWX*WTMW)
          XGW_CWX = (XMGW_CWX*WTMW)/(XMGA_CWX*WTMA + XMGW_CWX*WTMW)
          CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
          ISRX = 2
          CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
          RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
        ENDIF
!
!---    Load local variable for mass flux, kg/s  ---
!
        IF( M.EQ.1 ) THEN
!
!---      Mass flow rate, kg/s  ---
!
          VAR_CW2X = VAR_CWX(2)
!
!---      Load injection mass flux for use in RSDL_COUP_WELL ---
!
          FX_CW(NCW) = VAR_CW2X
        ENDIF
!
!---    Store top of coupled-well location in previous
!       coupled-well node location  ---
!
        XPX(1) = XTP_CW(1,ID_CW(1,NCW))
        YPX(1) = YTP_CW(1,ID_CW(1,NCW))
        ZPX(1) = ZTP_CW(1,ID_CW(1,NCW))
!
!---    Loop over the nodes in the coupled well ---
!
        DO 200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(NWN)
          I = ID(N)
          INVX = INV_CW(NWN)
          IZN = IZ(N)
          T_CWX = T(2,N)
!
!---      Coupled-well node centroids and projections ---
!
          XLX = PLX_CW(NWN)
          YLX = PLY_CW(NWN)
          ZLX = PLZ_CW(NWN)
          XPX(2) = 5.D-1*(XP_CW(2,NWN)+XP_CW(1,NWN))
          YPX(2) = 5.D-1*(YP_CW(2,NWN)+YP_CW(1,NWN))
          ZPX(2) = 5.D-1*(ZP_CW(2,NWN)+ZP_CW(1,NWN))
!
!---      Cylindrical coordinates with azimuthal symmetry,
!         centrally located wells  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            XPNX = 0.D+0
            YPNX = 0.D+0
            ZPNX = ZP(N)
!
!---      Cylindrical coordinates  ---
!
          ELSEIF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
            XPNX = XP(N)*COS(YP(N))
            YPNX = XP(N)*SIN(YP(N))
            ZPNX = ZP(N)
!
!---      Cartesian or boundary-fitted orthogonal coordinates  ---
!
          ELSE
            XPNX = XP(N)
            YPNX = YP(N)
            ZPNX = ZP(N)
          ENDIF
!
!---      Well pressure using previous coupled-well node density ---
!
          P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOG_CWX
!
!---      Well pressure at the node centroid, used for coupled-well
!         nodal output  ---
!
          IF( M.EQ.1 ) THEN
            NWF = IWP_CW(NWN)
            PF_CW(NWF) = P_CWX - (ZPNX-ZPX(1))*GRAV*RHOG_CWX
          ENDIF
!
!---      Adjust the formation pressure to the coupled-well node
!         centroid  ---
!
          IF( (SG(MF,N)-SGT(MF,N)).GT.EPSL ) THEN
            PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(MF,N)
          ELSE
            PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
          ENDIF
          PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
!
!---      Update coupled-well node density, with fixed solvated
!         water mass fraction  ---
!
          XMGA_CWX = (XGA_CWX/WTMA)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
          XMGW_CWX = (XGW_CWX/WTMW)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
          PVA_CWX = XMGA_CWX*(P_CWX+PATM)
          PVW_CWX = XMGW_CWX*(P_CWX+PATM)
          CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
          ISRX = 2
          CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
          RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---      Gas viscosity  ---
!
          CALL VISC_A( T_CWX,RHOGA_CWX,VISGA_CWX )
          CALL VISC_W( T_CWX,PVW_CWX,RHOGW_CWX,VISGW_CWX )
          CALL VISC_G( VISGA_CWX,VISGW_CWX,XMGA_CWX,XMGW_CWX,
     &      VISG_CWX )
!
!---      Equivalent field node radius components  ---
!
          PERMX = MAX( PERMV(1,N),1.D-20 )
          PERMY = MAX( PERMV(2,N),1.D-20 )
          PERMZ = MAX( PERMV(3,N),1.D-20 )
          RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---      Cylindrical coordinates with azimuthal symmetry,
!         centrally located wells  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            ROZ = RP(I)
            RWX = MIN( RWX,9.999D-1*ROZ )
            PERMX = PERMRF(MF,N)*PERMV(1,N)
            WI_CWX = 2.D+0*GPI*PERMX*ZLX/(LOG(ROZ/RWX)+PAR_CW(1,INVX))
          ELSE
            PERMYZ = SQRT(PERMY/PERMZ)
            PERMZY = SQRT(PERMZ/PERMY)
            DXGFX = DXGF(N)/FF_CW(1,NCW)
            DYGFX = DYGF(N)*RP(I)/FF_CW(2,NCW)
            DZGFX = DZGF(N)/FF_CW(3,NCW)
            ROX = 2.8D-1*SQRT(PERMYZ*(DZGFX**2) + PERMZY*(DYGFX**2))
     &      /(SQRT(PERMYZ)+SQRT(PERMZY))
            PERMZX = SQRT(PERMZ/PERMX)
            PERMXZ = SQRT(PERMX/PERMZ)
            ROY = 2.8D-1*SQRT(PERMZX*(DXGFX**2) + PERMXZ*(DZGFX**2))
     &        /(SQRT(PERMZX)+SQRT(PERMXZ))
            PERMYX = SQRT(PERMY/PERMX)
            PERMXY = SQRT(PERMX/PERMY)
            ROZ = 2.8D-1*SQRT(PERMYX*(DXGFX**2) + PERMXY*(DYGFX**2))
     &        /(SQRT(PERMYX)+SQRT(PERMXY))
!
!---        Well index components  ---
!
            PERMX = PERMRF(MF,N)*PERMV(1,N)
            PERMY = PERMRF(MF,N)*PERMV(2,N)
            PERMZ = PERMRF(MF,N)*PERMV(3,N)
            WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/
     &        (LOG(ROX/RWX)+PAR_CW(1,INVX))
            WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/
     &        (LOG(ROY/RWX)+PAR_CW(1,INVX))
            WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &        (LOG(ROZ/RWX)+PAR_CW(1,INVX))
            WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
          ENDIF
!
!---      Mass fluxes, positive into the node  ---
!
          DPGX = MAX( P_CWX-PGFX,0.D+0 )
          DPLX = MAX( P_CWX-MAX(PLFX,PGFX),0.D+0 )
!
!---      Zero fluxes from well to reservoir  ---
!
          FXW_CW(M,NWN) = 0.D+0
          FXS_CW(M,NWN) = 0.D+0
          FXA_CW(M,NWN) = 0.D+0
!
!---      Flux from well to reservoir  ---
!
          FX_CWX = WI_CWX*RHOG_CWX*DPGX/VISG_CWX
          FXA_CW(M,NWN) = FX_CWX*XGA_CWX
          FXW_CW(M,NWN) = FX_CWX*XGW_CWX
!
!---      Volumetric injection well fluxes  ---
!
          Q_CW(1,NWN) = FX_CWX/RHOG_CWX
          Q_CW(2,NWN) = 0.D+0
          Q_CW(3,NWN) = FX_CWX/RHOG_CWX
!
!---      Store current coupled-well node location in previous
!         coupled-well node location  ---
!
          XPX(1) = XPX(2)
          YPX(1) = YPX(2)
          ZPX(1) = ZPX(2)
  200   CONTINUE
  300 CONTINUE
!
!---  Mass balance residuals for injection type coupled well  ---
!
      NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
      MX = (NWFX*ISVC)+2
      RS_CW(1,NCW) = VAR_CWX(2)
      RS_CW(MX,NCW) = VAR_CWX(2)
      QM_CW(1,NCW) = 0.D+0
      QM_CW(3,NCW) = 0.D+0
      QM_CW(7,NCW) = 0.D+0
      QM_CWX = 0.D+0
!
!---  Loop over coupled-well nodes  ---
!
      DO 330 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
        RS_CW(1,NCW) = RS_CW(1,NCW) - FXA_CW(1,NWN)
     &    - FXW_CW(1,NWN) - FXS_CW(1,NWN)
        RS_CW(MX,NCW) = RS_CW(MX,NCW) - FXA_CW(ISVC+2,NWN)
     &    - FXW_CW(ISVC+2,NWN) - FXS_CW(ISVC+2,NWN)
        QM_CW(1,NCW) = QM_CW(1,NCW) + FXA_CW(1,NWN)
        QM_CW(3,NCW) = QM_CW(3,NCW) + FXW_CW(1,NWN)
        QM_CW(7,NCW) = QM_CW(7,NCW) + FXS_CW(1,NWN)
        QM_CWX = QM_CWX + FXA_CW(ISVC+2,NWN)
     &    + FXW_CW(ISVC+2,NWN) + FXS_CW(ISVC+2,NWN)
  330 CONTINUE
!
!---  Loop over field nodes that contain coupled-well nodes  ---
!
      DO 470 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
        M1 = (NWF-ID_CW(5,NCW))*ISVC + 1
        DO 460 M2 = 1,ISVC
          M3 = M1+M2
          RS_CW(M3,NCW) = VAR_CWX(2)
!
!---      Loop over coupled-well nodes  ---
!
          DO 450 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---        If coupled-well node is within the current field
!           node, use incremented fluxes  ---
!
            IF( IWF_CW(NWF).EQ.IWN_CW(NWN) ) THEN
              RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXA_CW(M2+1,NWN)
     &          - FXW_CW(M2+1,NWN) - FXS_CW(M2+1,NWN)
!
!---        If coupled-well node is outside the current field
!           node, use un-incremented fluxes  ---
!
            ELSE
              RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXA_CW(1,NWN)
     &          - FXW_CW(1,NWN) - FXS_CW(1,NWN)
            ENDIF
  450     CONTINUE
  460   CONTINUE
  470 CONTINUE
!
!---  Exit point  ---
!
  500 CONTINUE
!
!---  Dynamic memory deallocation  ---
!
      IF( ALLOCATED(MCW) ) THEN
      DEALLOCATE( MCW,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: MCW'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      IF( ALLOCATED(MFD) ) THEN
      DEALLOCATE( MFD,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: MFD'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INJ_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INJP_COUP_WELL( NCW )
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
!     STOMP-CO2
!
!     Injection coupled well model
!     
!     Rate controlled or pressure controlled
!
!     Flux of energy, water mass, and CO2 mass
!     from coupled-well nodes to field nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 March 2016.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYST
      USE GRID
      USE FDVP
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
      REAL*8 XPX(2),YPX(2),ZPX(2)
      REAL*8 VAR_CWX(6)
      REAL*8 FXA_CWX(3,LWN_CW),FXW_CWX(3,LWN_CW),FXS_CWX(3,LWN_CW)
      REAL*8 QM_CWX(3),P_CWY(3),VAR_CW2X(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INJP_COUP_WELL'
!
!---  Loop over coupled-well nodes  ---
!
      DO 30 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---    Loop over increment indices  ---
!
        DO 20 M = 1,2
          FXA_CWX(M,NWN) = 0.D+0
          FXS_CWX(M,NWN) = 0.D+0
          FXW_CWX(M,NWN) = 0.D+0
   20   CONTINUE
   30  CONTINUE
!
!---  Injection well time interval ---
!
      TMZ = TM
      IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
      NCT = 0
      DO M = 1,IM_CW(NCW)
        NCT = NCT + IMP_CW(M,NCW)
      ENDDO
      IF( ICC_CW(NCW).EQ.1 ) TMZ = MOD( TM,VAR_CW(1,NCT,NCW) )
!
!---  Coupled well is inactive set well pressure to be in 
!     equilibrium with formation  ---
!
      IF( TMZ.LE.VAR_CW(1,1,NCW) .OR. 
     &  QM_CW(2,NCW).GE.TML_CW(NCW) ) THEN
        CALL EQUIL_COUP_WELL( NCW )
        ID_CW(8,NCW) = 1
        GOTO 500
      ENDIF
      IF( NCT.EQ.1 ) THEN
        DO 80 N = 2,6
          VAR_CWX(N) = VAR_CW(N,1,NCW)
   80   CONTINUE
!
!---    Limit injection rate by total injected mass  ---
!
        VAR_CWX(2) = MIN( VAR_CWX(2),
     &    ((TML_CW(NCW)-QM_CW(2,NCW))*DTI) )
!
!---    Set well state  ---
!
        IT_CWX = IT_CW(NCW)
        IF( IT_CW(NCW).EQ.100 ) IT_CWX = ITS_CW(1,NCW)
      ELSE
        DO 100 M = 2,NCT
          IF( TMZ.LE.VAR_CW(1,M,NCW) ) THEN
            TD_CW = VAR_CW(1,M,NCW)-VAR_CW(1,M-1,NCW)
            DT_CW = MIN( VAR_CW(1,M,NCW)-TMZ,DT )
            TF_CW = (TMZ-VAR_CW(1,M-1,NCW))/TD_CW
            DO 90 N = 2,6
              VAR_CWX(N) = VAR_CW(N,M-1,NCW) + 
     &          TF_CW*(VAR_CW(N,M,NCW)-VAR_CW(N,M-1,NCW))
   90       CONTINUE
!
!---        Limit injection rate by total injected mass  ---
!
            VAR_CWX(2) = MIN( VAR_CWX(2),
     &        ((TML_CW(NCW)-QM_CW(2,NCW))*DTI) )
!
!---        Set well state  ---
!
            IT_CWX = IT_CW(NCW)
            IF( IT_CW(NCW).EQ.100 ) THEN
              NC = 0
              DO N = 1,IM_CW(NCW)
                NC = NC + IMP_CW(N,NCW)
                IF( NC.GE.M ) THEN
                  IT_CWX = ITS_CW(N,NCW)
                  EXIT
                ENDIF
              ENDDO
            ENDIF
            GOTO 122
          ENDIF
  100   CONTINUE
!
!---    Injection well is inactive set well pressure to be in 
!       equilibrium with reservoir  ---
!
        CALL EQUIL_COUP_WELL( NCW )
        ID_CW(8,NCW) = 1
        GOTO 500
      ENDIF
  122 CONTINUE  
!
!---  Load pressure limit ---
!
      PL_CW(NCW) = VAR_CWX(3) - PATM
!
!---  Upper pressure limit ---
!
      PL_CWX = VAR_CWX(3) - PATM
      P_CWY(1) = P_CW(2,NCW)
      P_CWY(2) = PL_CWX
      DP_CWX = 1.D-1
      ICHK_CWX = 0
      ML = 1
      MU = 2
      NC = 0
  124 CONTINUE  
      NC = NC + 1
!
!---  Loop over increment indices ---
!
      DO 300 M = ML,MU
        N = IWN_CW(ID_CW(3,NCW))
        P_CWX = P_CWY(M)
        T_CWX = T(2,N)
!
!---    CO2 injection with zero solvated water  ---
!
        IF( IT_CWX.EQ.1 ) THEN
          PVW_CWX = 0.D+0
          PVA_CWX = MAX( P_CWX+PATM-PVW_CWX,0.D+0 )
          XMGA_CWX = 1.D+0
          XMGW_CWX = 0.D+0
          XGA_CWX = 1.D+0
          XGW_CWX = 0.D+0
          CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
          ISRX = 2
          CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
          RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---    CO2 injection with solvated water concentration 
!       declared as mass fraction  ---
!
        ELSEIF( IT_CWX.EQ.2 ) THEN
          XGW_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
          XGA_CWX = 1.D+0-XGW_CWX
          XMGA_CWX = (XGA_CWX/WTMA)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
          XMGW_CWX = (XGW_CWX/WTMW)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
          PVA_CWX = XMGA_CWX*(P_CWX+PATM)
          PVW_CWX = XMGW_CWX*(P_CWX+PATM)
          CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
          ISRX = 2
          CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
          RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---    CO2 injection with solvated water concentration 
!       declared as relative humidity  ---
!
        ELSEIF( IT_CWX.EQ.3 ) THEN
          CALL SP_W( T_CWX,PSW_CWX )
          PVW_CWX = PSW_CWX*VAR_CWX(4)
          PVA_CWX = MAX( P_CWX+PATM-PVW_CWX,0.D+0 )
          XMGA_CWX = PVA_CWX/(P_CWX+PATM)
          XMGW_CWX = PVW_CWX/(P_CWX+PATM)
          XGA_CWX = (XMGA_CWX*WTMA)/(XMGA_CWX*WTMA + XMGW_CWX*WTMW)
          XGW_CWX = (XMGW_CWX*WTMW)/(XMGA_CWX*WTMA + XMGW_CWX*WTMW)
          CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
          ISRX = 2
          CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
          RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---    Aqueous injection with zero dissolved CO2  ---
!
        ELSEIF( MOD(IT_CWX,10).EQ.4 ) THEN
          PL_CWX = MAX( P_CWX+PATM,0.D+0 )
          XLA_CWX = 0.D+0
!
!---      Aqueous injection with zero dissolved salt  ---
!
          IF( IT_CWX/10.EQ.1 ) THEN
            XLS_CWX = 0.D+0
!
!---      Aqueous injection with dissolved salt declared as 
!         brine mass fraction  ---
!
          ELSEIF( IT_CWX/10.EQ.2 ) THEN
            XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
!
!---      Aqueous injection with dissolved salt declared as 
!         brine relative saturation  ---
!
          ELSEIF( IT_CWX/10.EQ.3 ) THEN
            CALL SOL_LS( T_CWX,XLSMX )
            XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
          ENDIF
          CALL SP_B( T_CWX,XLS_CWX,PSBX )
          PVBX = PSBX
          CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
          PVWX = PVBX
          XLSSX = XLS_CWX
          CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &      XMGWX,XMLAX,XMLSX,XMLWX )
          XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
          XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
          CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
!
!---    Aqueous injection with dissolved CO2 declared as 
!       mass fraction  ---
!
        ELSEIF( MOD(IT_CWX,10).EQ.5 ) THEN
          PL_CWX = MAX( P_CWX+PATM,0.D+0 )
          XLA_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
!
!---      Aqueous injection with zero dissolved salt  ---
!
          IF( IT_CWX/10.EQ.1 ) THEN
            XLS_CWX = 0.D+0
!
!---      Aqueous injection with dissolved salt declared as 
!         brine mass fraction  ---
!
          ELSEIF( IT_CWX/10.EQ.2 ) THEN
            XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
!
!---      Aqueous injection with dissolved salt declared as 
!         brine relative saturation  ---
!
          ELSEIF( IT_CWX/10.EQ.3 ) THEN
            CALL SOL_LS( T_CWX,XLSMX )
            XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
          ENDIF
          CALL SP_B( T_CWX,XLS_CWX,PSBX )
          PVBX = PSBX
          CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
          PVWX = PVBX
          XLSSX = XLS_CWX
          CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &      XMGWX,XMLAX,XMLSX,XMLWX )
          XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
          XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
          CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
!
!---    Aqueous injection with dissolved CO2 declared as 
!       relative saturation  ---
!
        ELSEIF( MOD(IT_CWX,10).EQ.6 ) THEN
          PL_CWX = MAX( P_CWX+PATM,0.D+0 )
          XLA_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
!
!---      Aqueous injection with zero dissolved salt  ---
!
          IF( IT_CWX/10.EQ.1 ) THEN
            XLS_CWX = 0.D+0
!
!---      Aqueous injection with dissolved salt declared as 
!         brine mass fraction  ---
!
          ELSEIF( IT_CWX/10.EQ.2 ) THEN
            XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
 !
!---      Aqueous injection with dissolved salt declared as 
!         brine relative saturation  ---
!
          ELSEIF( IT_CWX/10.EQ.3 ) THEN
            CALL SOL_LS( T_CWX,XLSMX )
            XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
          ENDIF
          CALL SP_B( T_CWX,XLS_CWX,PSBX )
          PVBX = PSBX
          CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
          PVWX = PVBX
          XLSSX = XLS_CWX
          CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &      XMGWX,XMLAX,XMLSX,XMLWX )
          XLA_CWX = XLA_CWX*XLAX
          XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
          XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
          CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
        ENDIF
!
!---    Mass flow rate, kg/s  ---
!
        VAR_CW2X(M) = VAR_CWX(2)
!
!---    Store top of coupled-well location in previous
!       coupled-well node location  ---
!
        XPX(1) = XTP_CW(1,ID_CW(1,NCW))
        YPX(1) = YTP_CW(1,ID_CW(1,NCW))
        ZPX(1) = ZTP_CW(1,ID_CW(1,NCW))
!
!---    Loop over the nodes in the coupled well ---
!
        DO 200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(NWN)
          I = ID(N)
          INVX = INV_CW(NWN)
          IZN = IZ(N)
          T_CWX = T(2,N)
!
!---      Coupled-well node centroids and projections ---
!
          XLX = PLX_CW(NWN)
          YLX = PLY_CW(NWN)
          ZLX = PLZ_CW(NWN)
          XPX(2) = 5.D-1*(XP_CW(2,NWN)+XP_CW(1,NWN))
          YPX(2) = 5.D-1*(YP_CW(2,NWN)+YP_CW(1,NWN))
          ZPX(2) = 5.D-1*(ZP_CW(2,NWN)+ZP_CW(1,NWN))
!
!---      Cylindrical coordinates with azimuthal symmetry,
!         centrally located wells  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            XPNX = 0.D+0
            YPNX = 0.D+0
            ZPNX = ZP(N)
!
!---      Cylindrical coordinates  ---
!
          ELSEIF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
            XPNX = XP(N)*COS(YP(N))
            YPNX = XP(N)*SIN(YP(N))
            ZPNX = ZP(N)
!
!---      Cartesian or boundary-fitted orthogonal coordinates  ---
!
          ELSE
            XPNX = XP(N)
            YPNX = YP(N)
            ZPNX = ZP(N)
          ENDIF
!
!---      Well pressure using previous coupled-well node density ---
!
          IF( IT_CWX.LE.3 ) THEN
            P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOG_CWX
          ELSE
            P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOL_CWX
          ENDIF
!
!---      Adjust the formation pressure to the coupled-well node
!         centroid  ---
!
          IF( (SG(2,N)-SGT(2,N)).GT.EPSL ) THEN
            PGFX = PG(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(2,N)
          ELSE
            PGFX = PG(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
          ENDIF
          PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
!
!---      Update coupled-well node density and viscosity
!         for CO2 injection, by using a fixed solvated
!         water mass fraction  ---
!
          IF( IT_CWX.LE.3 ) THEN
            XMGA_CWX = (XGA_CWX/WTMA)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
            XMGW_CWX = (XGW_CWX/WTMW)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
            PVA_CWX = XMGA_CWX*(P_CWX+PATM)
            PVW_CWX = XMGW_CWX*(P_CWX+PATM)
            CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
            ISRX = 2
            CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
            RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
            CALL VISC_A( T_CWX,RHOGA_CWX,VISGA_CWX )
            CALL VISC_W( T_CWX,PVW_CWX,RHOGW_CWX,VISGW_CWX )
            CALL VISC_G( VISGA_CWX,VISGW_CWX,XMGA_CWX,XMGW_CWX,
     &        VISG_CWX )
!
!---      Update coupled-well node density and viscosity
!         for aqueous injection, by using a fixed CO2 and salt
!         dissolved mass fraction  ---
!
          ELSE
            PL_CWX = MAX( P_CWX+PATM,0.D+0 )
            CALL SP_B( T_CWX,XLS_CWX,PSBX )
            PVAX = MAX( PL_CWX-PSBX,0.D+0 )
            CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
            WTMLX = 1.D+0/(XLA_CWX/WTMA+XLS_CWX/WTMS+XLW_CWX/WTMW)
            XMLA_CWX = XLA_CWX*WTMLX/WTMA
            CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
            ISRX = 1
            CALL DENS_W( T_CWX,PL_CWX,RHOLWX,RHOX,ISRX )
            CALL VISC_W( T_CWX,PL_CWX,RHOLWX,VISLWX )
            CALL VISC_B( T_CWX,XLS_CWX,VISLWX,VISBX )
            CALL DENS_A( T_CWX,PVAX,RHOGAX,I_VX )
            CALL VISC_A( T_CWX,RHOGAX,VISGAX )
            CALL VISC_L( XMLA_CWX,VISBX,VISGAX,VISL_CWX )
          ENDIF
!
!---      Equivalent field node radius components  ---
!
          PERMX = MAX( PERMV(1,N),1.D-20 )
          PERMY = MAX( PERMV(2,N),1.D-20 )
          PERMZ = MAX( PERMV(3,N),1.D-20 )
          RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---      Cylindrical coordinates with azimuthal symmetry,
!         centrally located wells  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            ROZ = RP(I)
            RWX = MIN( RWX,9.999D-1*ROZ )
            PERMX = PERMRF(2,N)*PERMV(1,N)
            WI_CWX = 2.D+0*GPI*PERMX*ZLX/(LOG(ROZ/RWX)+PAR_CW(1,INVX))
          ELSE
            PERMYZ = SQRT(PERMY/PERMZ)
            PERMZY = SQRT(PERMZ/PERMY)
            DXGFX = DXGF(N)/FF_CW(1,NCW)
            DYGFX = DYGF(N)*RP(I)/FF_CW(2,NCW)
            DZGFX = DZGF(N)/FF_CW(3,NCW)
            ROX = 2.8D-1*SQRT(PERMYZ*(DZGFX**2) + PERMZY*(DYGFX**2))
     &      /(SQRT(PERMYZ)+SQRT(PERMZY))
            PERMZX = SQRT(PERMZ/PERMX)
            PERMXZ = SQRT(PERMX/PERMZ)
            ROY = 2.8D-1*SQRT(PERMZX*(DXGFX**2) + PERMXZ*(DZGFX**2))
     &        /(SQRT(PERMZX)+SQRT(PERMXZ))
            PERMYX = SQRT(PERMY/PERMX)
            PERMXY = SQRT(PERMX/PERMY)
            ROZ = 2.8D-1*SQRT(PERMYX*(DXGFX**2) + PERMXY*(DYGFX**2))
     &        /(SQRT(PERMYX)+SQRT(PERMXY))
!
!---        Well index components  ---
!
            PERMX = PERMRF(2,N)*PERMV(1,N)
            PERMY = PERMRF(2,N)*PERMV(2,N)
            PERMZ = PERMRF(2,N)*PERMV(3,N)
            WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/
     &        (LOG(ROX/RWX)+PAR_CW(1,INVX))
            WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/
     &        (LOG(ROY/RWX)+PAR_CW(1,INVX))
            WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &        (LOG(ROZ/RWX)+PAR_CW(1,INVX))
            WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
          ENDIF
!
!---      Mass fluxes, positive into the node  ---
!
          DPGX = MAX( P_CWX-PGFX,0.D+0 )
          DPLX = MAX( P_CWX-MAX(PLFX,PGFX),0.D+0 )
!
!---      Zero fluxes from well to reservoir  ---
!
          FXW_CWX(M,NWN) = 0.D+0
          FXS_CWX(M,NWN) = 0.D+0
          FXA_CWX(M,NWN) = 0.D+0
!
!---      CO2 injection, flux from well to formation  ---
!
          IF( IT_CWX.LE.3 ) THEN
            FX_CWX = WI_CWX*RHOG_CWX*DPGX/VISG_CWX
            FXA_CW(M,NWN) = FX_CWX*XGA_CWX
            FXW_CW(M,NWN) = FX_CWX*XGW_CWX
!
!---      Aqueous injection, flux from well to formation  ---
!
          ELSE
            FX_CWX = WI_CWX*RHOL_CWX*DPLX/VISL_CWX
            FXA_CW(M,NWN) = FX_CWX*XLA_CWX
            FXS_CW(M,NWN) = FX_CWX*XLS_CWX
            FXW_CW(M,NWN) = FX_CWX*XLW_CWX
          ENDIF
!
!---      Store current coupled-well node location in previous
!         coupled-well node location  ---
!
          XPX(1) = XPX(2)
          YPX(1) = YPX(2)
          ZPX(1) = ZPX(2)
  200   CONTINUE
!
!---    Mass balance residuals for injection type coupled well  ---
!
        QM_CWX(M) = 0.D+0
!
!---    Loop over coupled-well nodes  ---
!
        DO 220 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          QM_CWX(M) = QM_CWX(M) + FXW_CWX(M,NWN)
     &     + FXA_CWX(M,NWN) + FXS_CWX(M,NWN)
  220   CONTINUE
  300 CONTINUE
!
!---  Consider current well pressure and well pressure limit  ---
!
      IF( ICHK_CWX.EQ.0 ) THEN
        ICHK_CWX = 1
!
!---    Hold pressure controlled option  ---
!
        IF( ID_CW(8,NCW).EQ.1 .AND. NITER.GE.4 ) THEN
          P_CW(2,NCW) = PL_CW(NCW)
          DNR_CW(NCW) = 1.D-1
          P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
          ID_CW(8,NCW) = 1
          GOTO 500
!
!---    Well-limit pressure insufficient to produce specified
!       rate, well is pressure controlled  ---
!
        ELSEIF( QM_CWX(2).LT.VAR_CW2X(2) ) THEN
          P_CW(2,NCW) = PL_CW(NCW)
          DNR_CW(NCW) = 1.D-1
          P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
          ID_CW(8,NCW) = 1
          GOTO 500
!
!---    Well-limit pressure and current well pressure yield
!       flow above specified rate, well is flow controlled.  Find
!       well pressure that yields positive flow below specified rate ---
!
        ELSEIF( QM_CWX(1).GE.VAR_CW2X(1) ) THEN
          DP_CWX = 1.D+1*DP_CWX
          P_CWY(1) = P_CWY(1) - DP_CWX
          ICHK_CWX = 0
          GOTO 124
!
!---    Well-limit pressure yields flow above specified rate,
!       and current well pressure yields positive flow, well is flow
!       controlled  ---
!
        ELSEIF( QM_CWX(1).GT.EPSL ) THEN
          DNR_CW(NCW) = 1.D-1
          P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
          ID_CW(8,NCW) = 0
          GOTO 500
!
!---    Well limit pressure yields flow above specified rate,
!       and current well pressure yields zero flow, use bisection to
!       determine a well pressure that yields the specified rate  ---
!
        ELSE
          P_CWY(3) = 5.D-1*(P_CWY(1)+P_CWY(2))
          ML = 3
          MU = 3
          GOTO 124
        ENDIF
!
!---    Use bisection to determine a well pressure that yields
!       the specified rate  ---
!
      ELSE
        IF( (ABS(QM_CWX(3)-VAR_CW2X(3)).LT.(1.D-3*VAR_CW2X(3))
     &    .OR. (P_CWY(2)-P_CWY(1)).LT.1.D-1) .OR. NC.GT.32 ) THEN
          P_CW(2,NCW) = P_CWY(3)
          DNR_CW(NCW) = 1.D-1
          P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
          ID_CW(8,NCW) = 0
          GOTO 500
        ELSEIF( ((QM_CWX(3)-VAR_CW2X(3)).LE.0.D+0 .AND.
     &    (QM_CWX(1)-VAR_CW2X(1)).LE.0.D+0) .OR. 
     &    ((QM_CWX(3)-VAR_CW2X(3)).GT.0.D+0 .AND.
     &    (QM_CWX(1)-VAR_CW2X(1)).GT.0.D+0) ) THEN
          P_CWY(1) = P_CWY(3)
          QM_CWX(1) = QM_CWX(3)
          VAR_CW2X(1) = VAR_CW2X(3)
          P_CWY(3) = 5.D-1*(P_CWY(1)+P_CWY(2))
          GOTO 124
        ELSE
          P_CWY(2) = P_CWY(3)
          QM_CWX(2) = QM_CWX(3)
          VAR_CW2X(2) = VAR_CW2X(3)
          P_CWY(3) = 5.D-1*(P_CWY(1)+P_CWY(2))
          GOTO 124
        ENDIF
      ENDIF
!
!---  Exit point  ---
!
  500 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INJP_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCB_COUP_WELL
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
!     STOMP-CO2
!
!     Modify Jacobian matrix for the coupled-well equations
!     and load Jacobian matrix for the coupled-well equations
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 20 April 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FDVP
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCB_COUP_WELL'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Loop over coupled wells ---
!
        DO 500 NCW = 1,N_CW
!
!---      Loop over coupled-well well nodes ---
!
          DO 200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
            N = IWN_CW(NWN)
            NMD = ABS(IXP(N))
!
!---        Water mass balance equation at field node ---
!
            MP = IM(IEQW,NMD)
!
!---        Change in water mass flux into field node with respect
!           to change in field node primary variables  ---
!
            DO 100 M = 1,ISVC
              DNRX = DNR(M,N)
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &          (FXW_CW(1,NWN)-FXW_CW(M+1,NWN))/DNRX
  100       CONTINUE
!
!---        Change in water mass flux into field node with respect
!           to change in coupled-well pressure  ---
!
            MCOL = JM_CW(NCW)
            MROW = MP-MCOL+MDC
            MX = ISVC+2
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &        (FXW_CW(1,NWN)-FXW_CW(MX,NWN))/DNR_CW(NCW)
!
!---        Water mass flux into field node, kg/s  ---
!
            BLU(MP) = BLU(MP) + FXW_CW(1,NWN)
            RSDL(IEQW,N) = BLU(MP)
!
!---        CO2 mass balance equation at field node ---
!
            MP = IM(IEQA,NMD)
!
!---        Change in CO2 mass flux into field node with respect
!           to change in field node primary variables  ---
!
            DO 110 M = 1,ISVC
              DNRX = DNR(M,N)
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &          (FXA_CW(1,NWN)-FXA_CW(M+1,NWN))/DNRX
  110       CONTINUE
!
!---        Change in CO2 mass flux into field node with respect
!           to change in coupled-well pressure  ---
!
            MCOL = JM_CW(NCW)
            MROW = MP-MCOL+MDC
            MX = ISVC+2
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &        (FXA_CW(1,NWN)-FXA_CW(MX,NWN))/DNR_CW(NCW)
!
!---        CO2 mass flux into field node, kg/s  ---
!
            BLU(MP) = BLU(MP) + FXA_CW(1,NWN)
            RSDL(IEQA,N) = BLU(MP)
!
!---        Skip for isobrine simulations  ---
!
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Salt mass balance equation at field node ---
!
              MP = IM(IEQS,NMD)
!
!---          Change in salt mass flux into field node with respect
!             to change in field node primary variables  ---
!
              DO 120 M = 1,ISVC
                DNRX = DNR(M,N)
                MCOL = IM(M,NMD)
                MROW = MP-MCOL+MDC
                ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &            (FXS_CW(1,NWN)-FXS_CW(M+1,NWN))/DNRX
  120         CONTINUE
!
!---          Change in salt mass flux into field node with respect
!             to change in coupled-well pressure  ---
!
              MCOL = JM_CW(NCW)
              MROW = MP-MCOL+MDC
              MX = ISVC+2
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &          (FXS_CW(1,NWN)-FXS_CW(MX,NWN))/DNR_CW(NCW)
!
!---          Salt mass flux into field node, kg/s  ---
!
              BLU(MP) = BLU(MP) + FXS_CW(1,NWN)
              RSDL(IEQS,N) = BLU(MP)
            ENDIF
  200     CONTINUE
!
!---      Coupled-well mass balance  ---
!
          MP = JM_CW(NCW)
          BLU(MP) = -RS_CW(1,NCW)
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) BLU(MP) = 0.D+0
!
!---      Change in coupled-well mass balance with respect to
!         change in coupled-well pressure  ---
!
          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
          MX = (NWFX*ISVC)+2
          MCOL = JM_CW(NCW)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &      (RS_CW(MX,NCW)-RS_CW(1,NCW))/DNR_CW(NCW)
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) ALU(MROW,MCOL) = 1.D+0
!
!---      Loop over field nodes with coupled-well nodes ---
!
          DO 400 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
            N = IWF_CW(NWF)
            NMD = ABS(IXP(N))
            MX = (NWF-ID_CW(5,NCW))*ISVC + 1
!
!---        Change in coupled-well mass balance with respect to
!           change in field node primary variables  ---
!
            DO 300 M = 1,ISVC
              DNRX = DNR(M,N)
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &          (RS_CW(MX+M,NCW)-RS_CW(1,NCW))/DNRX
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) ALU(MROW,MCOL) = 0.D+0
  300       CONTINUE
  400     CONTINUE
  500   CONTINUE
!
!---  SPLib solver  ---
!
      ELSEIF( ILES.EQ.3 ) THEN
!
!---    Loop over coupled wells ---
!
        DO 1500 NCW = 1,N_CW
!
!---      Loop over coupled-well well nodes ---
!
          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
          DO 1200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
            N = IWN_CW(NWN)
            NMD = ABS(IXP(N))
!
!---        Water mass balance equation at field node ---
!
            MP = IM(IEQW,NMD)
!
!---        Change in water mass flux into field node with respect
!           to change in field node primary variables  ---
!
            MA = 0
            DO 1100 M = 1,ISVC
              DNRX = DNR(M,N)
              MCOL = KLU(MP,M+MA)
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXW_CW(1,NWN)-FXW_CW(M+1,NWN))/DNRX
 1100       CONTINUE
!
!---        Change in water fraction of coupled-well flux with respect
!           to change in coupled-well pressure  ---
!
!            MA = NWFX*ISVC + (IWP_CW(NWN)-1)*ISVC + IEQW + 1
            MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + IEQW + 1
            MCOL = KLU_CW(MA,NCW)
            MX = ISVC+2
            DLU(MCOL) = DLU(MCOL) + 
     &        (FXW_CW(1,NWN)-FXW_CW(MX,NWN))/DNR_CW(NCW)
!
!---        Water mass flux into field node, kg/s  ---
!
            BLU(MP) = BLU(MP) + FXW_CW(1,NWN)
            RSDL(IEQW,N) = BLU(MP)
!
!---        CO2 mass balance equation at field node ---
!
            MP = IM(IEQA,NMD)
!
!---        Change in CO2 fraction of coupled-well flux with respect
!           to change in field node primary variables  ---
!
            MA = 0 
            DO 1110 M = 1,ISVC
              MCOL = KLU(MP,M+MA)
              DNRX = DNR(M,N)
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXA_CW(1,NWN)-FXA_CW(M+1,NWN))/DNRX
!            PRINT *,'NWN = ',NWN,' M = ',M,' N = ',N,
!     &        '(FXA_CW(1,NWN)-FXA_CW(M+1,NWN))/DNRX = ',
!     &        (FXA_CW(1,NWN)-FXA_CW(M+1,NWN))/DNRX,
!     &        ' PL(2,N) = ',PL(2,N),' PG(2,N) = ',PG(2,N),
!     &        ' SG(2,N) = ',SG(2,N)

 1110       CONTINUE
!
!---        Change in CO2 fraction of coupled-well flux with respect
!           to change in coupled-well pressure  ---
!
!            MA = NWFX*ISVC + (IWP_CW(NWN)-1)*ISVC + IEQA + 1
            MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + IEQA + 1
            MCOL = KLU_CW(MA,NCW)
            MX = ISVC+2
            DLU(MCOL) = DLU(MCOL) + 
     &        (FXA_CW(1,NWN)-FXA_CW(MX,NWN))/DNR_CW(NCW)
!            PRINT *,'NWN = ',NWN,
!     &        '(FXA_CW(1,NWN)-FXA_CW(MX,NWN))/DNR_CW(NCW) = ',
!     &        (FXA_CW(1,NWN)-FXA_CW(MX,NWN))/DNR_CW(NCW)
!
!---        CO2 mass flux into field node, kg/s  ---
!
            BLU(MP) = BLU(MP) + FXA_CW(1,NWN)
            RSDL(IEQA,N) = BLU(MP)
!            PRINT *,'NWN = ',NWN,'FXA_CW(1,NWN) = ',FXA_CW(1,NWN)
!
!---        Skip for isobrine simulations  ---
!
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Salt mass balance equation at field node ---
!
              MP = IM(IEQS,NMD)
!
!---          Change in salt mass flux into field node with respect
!             to change in field node primary variables  ---
!
              MA = 0
              DO 1120 M = 1,ISVC
                DNRX = DNR(M,N)
                MCOL = KLU(MP,M+MA)
                DLU(MCOL) = DLU(MCOL) + 
     &            (FXS_CW(1,NWN)-FXS_CW(M+1,NWN))/DNRX
 1120         CONTINUE
!
!---          Change in salt fraction of coupled-well flux with respect
!             to change in coupled-well pressure  ---
!
!              MA = NWFX*ISVC + (IWP_CW(NWN)-1)*ISVC + IEQS + 1
              MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC+ IEQS + 1
              MCOL = KLU_CW(MA,NCW)
              MX = ISVC+2
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXS_CW(1,NWN)-FXS_CW(MX,NWN))/DNR_CW(NCW)
!
!---          Water mass flux into field node, kg/s  ---
!
              BLU(MP) = BLU(MP) + FXS_CW(1,NWN)
              RSDL(IEQS,N) = BLU(MP)
            ENDIF
 1200     CONTINUE
!
!---      Coupled-well mass balance  ---
!
          MP = JM_CW(NCW)
          BLU(MP) = -RS_CW(1,NCW)
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) BLU(MP) = 0.D+0
!
!---      Change in coupled-well mass balance with respect to
!         change in coupled-well pressure  ---
!
          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
          MX = (NWFX*ISVC)+2
          MA = 1
          MCOL = KLU_CW(MA,NCW)
          DLU(MCOL) = DLU(MCOL) + 
     &      (RS_CW(MX,NCW)-RS_CW(1,NCW))/DNR_CW(NCW)
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 1.D+0
!
!---      Loop over field nodes with coupled-well nodes ---
!
          DO 1400 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
            N = IWF_CW(NWF)
            NMD = ABS(IXP(N))
            MX = (NWF-ID_CW(5,NCW))*ISVC + 1
            MA = (NWF-ID_CW(5,NCW))*ISVC + 1
!
!---        Change in coupled-well mass balance with respect to
!           change in field node primary variables  ---
!
            DO 1300 M = 1,ISVC
              DNRX = DNR(M,N)
              MCOL = KLU_CW(M+MA,NCW)
              DLU(MCOL) = DLU(MCOL) + 
     &          (RS_CW(MX+M,NCW)-RS_CW(1,NCW))/DNRX
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 0.D+0
 1300       CONTINUE
 1400     CONTINUE
 1500   CONTINUE
#ifdef lis
!
!---  Lis solver  ---
!
      ELSEIF( ILES.EQ.4 ) THEN
!
!---    Loop over coupled wells ---
!
        DO NCW = 1,N_CW
!
!---      Loop over coupled-well well nodes ---
!
          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
          DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
            N = IWN_CW(NWN)
            NMD = ABS(IXP(N))
!
!---        Water mass balance equation at field node ---
!
            MP = IM(IEQW,NMD)
!
!---        Change in water mass flux into field node with respect
!           to change in field node primary variables  ---
!
            MA = 0
            DO M = 1,ISVC
              DNRX = DNR(M,N)
              MCOL = KLU(MP,M+MA)
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXW_CW(1,NWN)-FXW_CW(M+1,NWN))/DNRX
            ENDDO
!
!---        Change in water fraction of coupled-well flux with respect
!           to change in coupled-well pressure  ---
!
!            MA = NWFX*ISVC + (IWP_CW(NWN)-1)*ISVC + IEQW + 1
            MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + IEQW + 1
            MCOL = KLU_CW(MA,NCW)
            MX = ISVC+2
            DLU(MCOL) = DLU(MCOL) + 
     &        (FXW_CW(1,NWN)-FXW_CW(MX,NWN))/DNR_CW(NCW)
!
!---        Water mass flux into field node, kg/s  ---
!
            BLU(MP) = BLU(MP) + FXW_CW(1,NWN)
            RSDL(IEQW,N) = BLU(MP)
!
!---        CO2 mass balance equation at field node ---
!
            MP = IM(IEQA,NMD)
!
!---        Change in CO2 fraction of coupled-well flux with respect
!           to change in field node primary variables  ---
!
            MA = 0 
            DO M = 1,ISVC
              MCOL = KLU(MP,M+MA)
              DNRX = DNR(M,N)
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXA_CW(1,NWN)-FXA_CW(M+1,NWN))/DNRX

            ENDDO
!
!---        Change in CO2 fraction of coupled-well flux with respect
!           to change in coupled-well pressure  ---
!
            MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + IEQA + 1
            MCOL = KLU_CW(MA,NCW)
            MX = ISVC+2
            DLU(MCOL) = DLU(MCOL) + 
     &        (FXA_CW(1,NWN)-FXA_CW(MX,NWN))/DNR_CW(NCW)
!
!---        CO2 mass flux into field node, kg/s  ---
!
            BLU(MP) = BLU(MP) + FXA_CW(1,NWN)
            RSDL(IEQA,N) = BLU(MP)
!            PRINT *,'NWN = ',NWN,'FXA_CW(1,NWN) = ',FXA_CW(1,NWN)
!
!---        Skip for isobrine simulations  ---
!
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Salt mass balance equation at field node ---
!
              MP = IM(IEQS,NMD)
!
!---          Change in salt mass flux into field node with respect
!             to change in field node primary variables  ---
!
              MA = 0
              DO M = 1,ISVC
                DNRX = DNR(M,N)
                MCOL = KLU(MP,M+MA)
                DLU(MCOL) = DLU(MCOL) + 
     &            (FXS_CW(1,NWN)-FXS_CW(M+1,NWN))/DNRX
              ENDDO
!
!---          Change in salt fraction of coupled-well flux with respect
!             to change in coupled-well pressure  ---
!
!              MA = NWFX*ISVC + (IWP_CW(NWN)-1)*ISVC + IEQS + 1
              MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC+ IEQS + 1
              MCOL = KLU_CW(MA,NCW)
              MX = ISVC+2
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXS_CW(1,NWN)-FXS_CW(MX,NWN))/DNR_CW(NCW)
!
!---          Water mass flux into field node, kg/s  ---
!
              BLU(MP) = BLU(MP) + FXS_CW(1,NWN)
              RSDL(IEQS,N) = BLU(MP)
            ENDIF
          ENDDO
!
!---      Coupled-well mass balance  ---
!
          MP = JM_CW(NCW)
          BLU(MP) = -RS_CW(1,NCW)
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) BLU(MP) = 0.D+0
!
!---      Change in coupled-well mass balance with respect to
!         change in coupled-well pressure  ---
!
          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
          MX = (NWFX*ISVC)+2
          MA = 1
          MCOL = KLU_CW(MA,NCW)
          DLU(MCOL) = DLU(MCOL) + 
     &      (RS_CW(MX,NCW)-RS_CW(1,NCW))/DNR_CW(NCW)
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 1.D+0
!
!---      Loop over field nodes with coupled-well nodes ---
!
          DO NWF = ID_CW(5,NCW),ID_CW(6,NCW)
            N = IWF_CW(NWF)
            NMD = ABS(IXP(N))
            MX = (NWF-ID_CW(5,NCW))*ISVC + 1
            MA = (NWF-ID_CW(5,NCW))*ISVC + 1
!
!---        Change in coupled-well mass balance with respect to
!           change in field node primary variables  ---
!
            DO M = 1,ISVC
              DNRX = DNR(M,N)
              MCOL = KLU_CW(M+MA,NCW)
              DLU(MCOL) = DLU(MCOL) + 
     &          (RS_CW(MX+M,NCW)-RS_CW(1,NCW))/DNRX
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 0.D+0
            ENDDO
          ENDDO
        ENDDO
#endif
#ifdef petsc
!
!---  PETSc solver  ---
!
      ELSEIF( ILES.EQ.5 ) THEN
!
!---    Loop over coupled wells ---
!
        DO NCW = 1,N_CW
!
!---      Loop over coupled-well well nodes ---
!
          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
          DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
            N = IWN_CW(NWN)
            NMD = ABS(IXP(N))
!
!---        Water mass balance equation at field node ---
!
            MP = IM(IEQW,NMD)
!
!---        Change in water mass flux into field node with respect
!           to change in field node primary variables  ---
!
            MA = 0
            DO M = 1,ISVC
              DNRX = DNR(M,N)
              MCOL = KLU(MP,M+MA)
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXW_CW(1,NWN)-FXW_CW(M+1,NWN))/DNRX
            ENDDO
!
!---        Change in water fraction of coupled-well flux with respect
!           to change in coupled-well pressure  ---
!
!            MA = NWFX*ISVC + (IWP_CW(NWN)-1)*ISVC + IEQW + 1
            MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + IEQW + 1
            MCOL = KLU_CW(MA,NCW)
            MX = ISVC+2
            DLU(MCOL) = DLU(MCOL) + 
     &        (FXW_CW(1,NWN)-FXW_CW(MX,NWN))/DNR_CW(NCW)
!
!---        Water mass flux into field node, kg/s  ---
!
            BLU(MP) = BLU(MP) + FXW_CW(1,NWN)
            RSDL(IEQW,N) = BLU(MP)
!
!---        CO2 mass balance equation at field node ---
!
            MP = IM(IEQA,NMD)
!
!---        Change in CO2 fraction of coupled-well flux with respect
!           to change in field node primary variables  ---
!
            MA = 0 
            DO M = 1,ISVC
              MCOL = KLU(MP,M+MA)
              DNRX = DNR(M,N)
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXA_CW(1,NWN)-FXA_CW(M+1,NWN))/DNRX

            ENDDO
!
!---        Change in CO2 fraction of coupled-well flux with respect
!           to change in coupled-well pressure  ---
!
            MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + IEQA + 1
            MCOL = KLU_CW(MA,NCW)
            MX = ISVC+2
            DLU(MCOL) = DLU(MCOL) + 
     &        (FXA_CW(1,NWN)-FXA_CW(MX,NWN))/DNR_CW(NCW)
!
!---        CO2 mass flux into field node, kg/s  ---
!
            BLU(MP) = BLU(MP) + FXA_CW(1,NWN)
            RSDL(IEQA,N) = BLU(MP)
!            PRINT *,'NWN = ',NWN,'FXA_CW(1,NWN) = ',FXA_CW(1,NWN)
!
!---        Skip for isobrine simulations  ---
!
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Salt mass balance equation at field node ---
!
              MP = IM(IEQS,NMD)
!
!---          Change in salt mass flux into field node with respect
!             to change in field node primary variables  ---
!
              MA = 0
              DO M = 1,ISVC
                DNRX = DNR(M,N)
                MCOL = KLU(MP,M+MA)
                DLU(MCOL) = DLU(MCOL) + 
     &            (FXS_CW(1,NWN)-FXS_CW(M+1,NWN))/DNRX
              ENDDO
!
!---          Change in salt fraction of coupled-well flux with respect
!             to change in coupled-well pressure  ---
!
!              MA = NWFX*ISVC + (IWP_CW(NWN)-1)*ISVC + IEQS + 1
              MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC+ IEQS + 1
              MCOL = KLU_CW(MA,NCW)
              MX = ISVC+2
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXS_CW(1,NWN)-FXS_CW(MX,NWN))/DNR_CW(NCW)
!
!---          Water mass flux into field node, kg/s  ---
!
              BLU(MP) = BLU(MP) + FXS_CW(1,NWN)
              RSDL(IEQS,N) = BLU(MP)
            ENDIF
          ENDDO
!
!---      Coupled-well mass balance  ---
!
          MP = JM_CW(NCW)
          BLU(MP) = -RS_CW(1,NCW)
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) BLU(MP) = 0.D+0
!
!---      Change in coupled-well mass balance with respect to
!         change in coupled-well pressure  ---
!
          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
          MX = (NWFX*ISVC)+2
          MA = 1
          MCOL = KLU_CW(MA,NCW)
          DLU(MCOL) = DLU(MCOL) + 
     &      (RS_CW(MX,NCW)-RS_CW(1,NCW))/DNR_CW(NCW)
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 1.D+0
!
!---      Loop over field nodes with coupled-well nodes ---
!
          DO NWF = ID_CW(5,NCW),ID_CW(6,NCW)
            N = IWF_CW(NWF)
            NMD = ABS(IXP(N))
            MX = (NWF-ID_CW(5,NCW))*ISVC + 1
            MA = (NWF-ID_CW(5,NCW))*ISVC + 1
!
!---        Change in coupled-well mass balance with respect to
!           change in field node primary variables  ---
!
            DO M = 1,ISVC
              DNRX = DNR(M,N)
              MCOL = KLU_CW(M+MA,NCW)
              DLU(MCOL) = DLU(MCOL) + 
     &          (RS_CW(MX+M,NCW)-RS_CW(1,NCW))/DNRX
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 0.D+0
            ENDDO
          ENDDO
        ENDDO
#endif
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCB_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCB_LEAK_WELL
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
!     Modify Jacobian matrix for flux between the leaky well nodes
!     and field nodes
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 20 October 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE LEAK_WELL
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXS
      USE FLUXP
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
      REAL*8 RTP(LUK),RTA(LUK),FT(LSFV)
      REAL*8 RWP(LUK),RWA(LUK),FW(LSFV)
      REAL*8 RAP(LUK),RAA(LUK),FA(LSFV)
      REAL*8 RSP(LUK),RSA(LUK),FS(LSFV)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCB_LEAK_WELL'
!
!---  Loop over number of leaky wells  ---
!
      DO NLW = 1,N_LW
!
!---    Loop over number of leaky well nodes in the leaky well  ---
!
        DO NWN = ID_LW(3,NLW),ID_LW(4,NLW)
          NW = ND_LW(NWN)
          N = NF_LW(NWN)
!
!---      Leaky well node is not connected to a field node  ---
!
          IF( N.EQ.0 ) CYCLE
          INVX = INV_LW(NWN)
          NPX = NSX(NW)
          NQX = NSX(NW) + 1
          NMD = ABS(IXP(N))
          NMDW = ABS(IXP(NW))
!
!---      Skip for isothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Energy flux between leaky well node and field node
!           expressed in terms of field node energy equation ---
!
            DO M = 1,ISVF
              FT(M) = -AFX(NPX)*UQ(M,NPX)
            ENDDO
!
!---        Energy equation residuals  ---
!
            RTS = FT(1)
            DO M = 1,ISVC
              MM = 2*M
              RTP(M) = FT(MM)
              MM = 2*M + 1
              RTA(M) = FT(MM)
            ENDDO
!
!---        Load Jacobian Matrix  ---
!
            CALL JCBBM_LEAK_WELL( RTS,RTP,RTA,N,NW,NWN,IEQT )
!
!---        Energy flux between leaky well node and field node
!           expressed in terms of leaky well node energy equation ---
!
            DO M = 1,ISVF
              MF = MFLX(M)
              FT(M) = AFX(NQX)*UQ(MF,NQX)
            ENDDO
!
!---        Energy equation residuals  ---
!
            RTS = FT(1)
            DO M = 1,ISVC
              MM = 2*M
              RTP(M) = FT(MM)
              MM = 2*M + 1
              RTA(M) = FT(MM)
            ENDDO
!
!---        Load Jacobian Matrix  ---
!
            CALL JCBBM_LEAK_WELL( RTS,RTP,RTA,NW,N,NWN,IEQT )
          ENDIF
!
!---      Water mass flux between leaky well node and field node
!         expressed in terms of field node water mass equation ---
!
          DO M = 1,ISVF
            MW = MADJ(M)
            MP = MNOD(M)
            FLWW = XLW(MW,NW)*RHOL(MW,NW)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWW,FLWP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            FGWW = XGW(MW,NW)*RHOG(MW,NW)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            FW(M) = -AFX(NPX)*(UL(M,NPX)*FLW + UG(M,NPX)*FGW
     &        + WTMW*UDGW(M,NPX) - WTMW*UDLA(M,NPX)
     &        - WTMW*UDS(M,NPX)/WTMS)
          ENDDO
!
!---      Water mass equation residuals  ---
!
          RWS = FW(1)
          DO M = 1,ISVC
            MM = 2*M
            RWP(M) = FW(MM)
            MM = 2*M + 1
            RWA(M) = FW(MM)
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBBM_LEAK_WELL( RWS,RWP,RWA,N,NW,NWN,IEQW )
!
!---      Water mass flux between leaky well node and field node
!         expressed in terms of leaky well node water mass equation ---
!
          DO M = 1,ISVF
            ME = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLWE = XLW(ME,N)*RHOL(ME,N)
            FLWP = XLW(MP,NW)*RHOL(MP,NW)
            INDX = 2
            FLW = DIFMN( FLWP,FLWE,DXGF(NW),DXGF(N),UL(1,NQX),INDX )
            FGWE = XGW(ME,N)*RHOG(ME,N)
            FGWP = XGW(MP,NW)*RHOG(MP,NW)
            INDX = 3
            FGW = DIFMN( FGWP,FGWE,DXGF(NW),DXGF(N),UG(1,NQX),INDX )
            FW(M) = AFX(NQX)*(UL(MF,NQX)*FLW + UG(MF,NQX)*FGW
     &        + WTMW*UDGW(MF,NQX) - WTMW*UDLA(MF,NQX)
     &        - WTMW*UDS(MF,NQX)/WTMS)
          ENDDO
!
!---      Water mass equation residuals  ---
!
          RWS = FW(1)
          DO M = 1,ISVC
            MM = 2*M
            RWP(M) = FW(MM)
            MM = 2*M + 1
            RWA(M) = FW(MM)
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBBM_LEAK_WELL( RWS,RWP,RWA,NW,N,NWN,IEQW )
!
!---      CO2 mass flux between leaky well node and field node
!         expressed in terms of field node CO2 mass equation ---
!
          DO M = 1,ISVF
            MW = MADJ(M)
            MP = MNOD(M)
            FLAW = XLA(MW,NW)*RHOL(MW,NW)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAW,FLAP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            FGAW = XGA(MW,NW)*RHOG(MW,NW)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAW,FGAP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            FA(M) = -AFX(NPX)*(UL(M,NPX)*FLA + UG(M,NPX)*FGA -
     &        WTMA*UDGW(M,NPX) + WTMA*UDLA(M,NPX))
          ENDDO
!
!---      CO2 mass equation residuals  ---
!
          RAS = FA(1)
          DO M = 1,ISVC
            MM = 2*M
            RAP(M) = FA(MM)
            MM = 2*M + 1
            RAA(M) = FA(MM)
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBBM_LEAK_WELL( RAS,RAP,RAA,N,NW,NWN,IEQA )
!
!---      CO2 mass flux between leaky well node and field node
!         expressed in terms of leaky well node CO2 mass equation ---
!
          DO M = 1,ISVF
            ME = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLAE = XLA(ME,N)*RHOL(ME,N)
            FLAP = XLA(MP,NW)*RHOL(MP,NW)
            INDX = 2
            FLA = DIFMN( FLAP,FLAE,DXGF(NW),DXGF(N),UL(1,NQX),INDX )
            FGAE = XGA(ME,N)*RHOG(ME,N)
            FGAP = XGA(MP,NW)*RHOG(MP,NW)
            INDX = 3
            FGA = DIFMN( FGAP,FGAE,DXGF(NW),DXGF(N),UG(1,NQX),INDX )
            FA(M) = AFX(NQX)*(UL(MF,NQX)*FLA + UG(MF,NQX)*FGA -
     &        WTMA*UDGW(MF,NQX) + WTMA*UDLA(MF,NQX))
          ENDDO
!
!---      CO2 mass equation residuals  ---
!
          RAS = FA(1)
          DO M = 1,ISVC
            MM = 2*M
            RAP(M) = FA(MM)
            MM = 2*M + 1
            RAA(M) = FA(MM)
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBBM_LEAK_WELL( RAS,RAP,RAA,NW,N,NWN,IEQA )
!
!---      Skip for isobrine simulations  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt mass flux between leaky well node and field node
!           expressed in terms of field node salt mass equation ---
!
            DO M = 1,ISVF
              FS(M) = -AFX(NPX)*US(M,NPX)
            ENDDO
!
!---        Salt mass equation residuals  ---
!
            RSS = FS(1)
            DO M = 1,ISVC
              MM = 2*M
              RSP(M) = FS(MM)
              MM = 2*M + 1
              RSA(M) = FS(MM)
            ENDDO
!
!---        Load Jacobian Matrix  ---
!
            CALL JCBBM_LEAK_WELL( RSS,RSP,RSA,N,NW,NWN,IEQS )
!
!---        Salt mass flux between leaky well node and field node
!           expressed in terms of leaky well node salt mass equation ---
!
            DO M = 1,ISVF
              MF = MFLX(M)
              FS(M) = AFX(NQX)*US(MF,NQX)
            ENDDO
!
!---        Salt mass equation residuals  ---
!
            RSS = FS(1)
            DO M = 1,ISVC
              MM = 2*M
              RSP(M) = FS(MM)
              MM = 2*M + 1
              RSA(M) = FS(MM)
            ENDDO
!
!---        Load Jacobian Matrix  ---
!
            CALL JCBBM_LEAK_WELL( RSS,RSP,RSA,NW,N,NWN,IEQS )
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCB_LEAK_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBBM_LEAK_WELL( RSS,RSP,RSA,NP,NA,NWN,MEQ )
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
!     Load the Jacobian matrix (banded matrix solver).
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 October 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE LEAK_WELL
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
      REAL*8 RSP(LUK),RSA(LUK)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBBM_LEAK_WELL'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Node (either field node or leaky well node)  ---
!
        NMD = IXP(NP)
        MP = IM(MEQ,NMD)
        DO M = 1,ISVC
          MCOL = IM(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSP(M)-RSS)/DNR(M,NP)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,NP) = BLU(MP)
!
!---    Adjacent node (either field node or leaky well node) ---
!
        NMD = IXP(NA)
        DO M = 1,ISVC
          MCOL = IM(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M)-RSS)/DNR(M,NA)
        ENDDO
!
!---  SPLib solver  ---
!
      ELSEIF( ILES.EQ.3 ) THEN
!
!---    Node (either field node or leaky well node)  ---
!
        NMD = IXP(NP)
        MP = IM(MEQ,NMD)
        MA = 0
        DO M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSP(M)-RSS)/DNR(M,NP)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,NP) = BLU(MP)
        MA = MA + ISVC
!
!---    Adjacent node (either field node or leaky well node) ---
!
        NMD = IXP(NA)
!
!---    Leaky well node is adjacent node  ---
!
        IF( ND_LW(NWN).EQ.NA ) THEN
          MP = (NWN-1)*ISVC + MEQ
!
!---    Leaky well node is principal node  ---
!
        ELSE
          MP = NWN_LW*ISVC + (NWN-1)*ISVC + MEQ
        ENDIF
        DO M = 1,ISVC
          MCOL = KLU_LW(MP,M)
          DLU(MCOL) = DLU(MCOL) + (RSA(M)-RSS)/DNR(M,NA)
        ENDDO
#ifdef lis
!
!---  Lis solver  ---
!
      ELSEIF( ILES.EQ.4 ) THEN
!
!---    Node (either field node or leaky well node)  ---
!
        NMD = IXP(NP)
        MP = IM(MEQ,NMD)
        MA = 0
        DO M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSP(M)-RSS)/DNR(M,NP)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,NP) = BLU(MP)
        MA = MA + ISVC
!
!---    Adjacent node (either field node or leaky well node) ---
!
        NMD = IXP(NA)
!
!---    Leaky well node is adjacent node  ---
!
        IF( ND_LW(NWN).EQ.NA ) THEN
          MP = (NWN-1)*ISVC + MEQ
!
!---    Leaky well node is principal node  ---
!
        ELSE
          MP = NWN_LW*ISVC + (NWN-1)*ISVC + MEQ
        ENDIF
        DO M = 1,ISVC
          MCOL = KLU_LW(MP,M)
          DLU(MCOL) = DLU(MCOL) + (RSA(M)-RSS)/DNR(M,NA)
        ENDDO
#endif
#ifdef petsc
!
!---  PETSc solver  ---
!
      ELSEIF( ILES.EQ.5 ) THEN
!
!---    Node (either field node or leaky well node)  ---
!
        NMD = IXP(NP)
        MP = IM(MEQ,NMD)
        MA = 0
        DO M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSP(M)-RSS)/DNR(M,NP)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,NP) = BLU(MP)
        MA = MA + ISVC
!
!---    Adjacent node (either field node or leaky well node) ---
!
        NMD = IXP(NA)
!
!---    Leaky well node is adjacent node  ---
!
        IF( ND_LW(NWN).EQ.NA ) THEN
          MP = (NWN-1)*ISVC + MEQ
!
!---    Leaky well node is principal node  ---
!
        ELSE
          MP = NWN_LW*ISVC + (NWN-1)*ISVC + MEQ
        ENDIF
        DO M = 1,ISVC
          MCOL = KLU_LW(MP,M)
          DLU(MCOL) = DLU(MCOL) + (RSA(M)-RSS)/DNR(M,NA)
        ENDDO
#endif
      ELSE
        INDX = 3
        CHMSG = 'Unknown Linear Equation Solver'
        CALL WRMSGS( INDX )
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBBM_LEAK_WELL group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO_COUP_WELL
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
!     STOMP-CO2
!
!     Load old-time-step values for coupled-well arrays.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 21 April 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LDO_COUP_WELL'
!
!---  Loop over coupled wells ---
!
      DO 100 NCW = 1,N_CW
        P_CW(1,NCW) = P_CW(2,NCW)
        PL_CW(NCW) = P_CW(2,NCW)
  100 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROJ_WELL( PLX_WX,PLY_WX,PLZ_WX,XP_WX,YP_WX,ZP_WX,N )
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
      SUB_LOG(ISUB_LOG) = '/PROJ_WELL'
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
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
!---  End of PROJ_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDCOUP_WELL
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
!     STOMP-CO2
!
!     Reads the coupled well card.
!     
!     FF_CW(1,LN_CW) - x-direction well fraction factor for well
!     FF_CW(2,LN_CW) - y-direction well fraction factor for well
!     FF_CW(3,LN_CW) - z-direction well fraction factor for well
!     ICC_CW(LN_CW) - cyclic well time index for well
!     IM_CW(LN_CW) - number of time points/periods for well
!     IMP_CW(LWTP_CW,LN_CW) - number of time points per time period
!     ID_CW(1,LN_CW) - starting well interval index for well
!     ID_CW(2,LN_CW) - ending well interval index for well
!     ID_CW(3,LN_CW) - starting well node index for well
!     ID_CW(4,LN_CW) - ending well node index for well
!     ID_CW(5,LN_CW) - principal well node index for well
!     IT_CW(LN_CW) - type index for well
!     1 - CO2 mass injection no water
!     2 - CO2 mass injection mass fraction water
!     3 - CO2 mass injection relative saturation water
!     4 - Aqueous mass injection no CO2
!     5 - Aqueous mass injection mass fraction CO2
!     6 - Aqueous mass injection relative saturation water
!   +10 - Aqueous mass injection no salt
!   +20 - Aqueous mass injection mass fraction salt
!   +30 - Aqueous mass injection relative saturation salt
!   100 - Alternating state injection
!     ITS_CW(LWTP_CW,LN_CW) - type index for injection well
!     1 - CO2 mass injection no water
!     2 - CO2 mass injection mass fraction water
!     3 - CO2 mass injection relative saturation water
!     4 - Aqueous mass injection no CO2
!     5 - Aqueous mass injection mass fraction CO2
!     6 - Aqueous mass injection relative saturation water
!   +10 - Aqueous mass injection no salt
!   +20 - Aqueous mass injection mass fraction salt
!   +30 - Aqueous mass injection relative saturation salt
!     JM_CW(LN_CW) - location of the well equation in the Jacobian
!     N_CW - number of coupled wells
!     PAR_CW(1,LWI_CW) - skin factor for well interval
!     TML_CW(LN_CW) - total mass injected/withdrawn limit for well, kg
!     VAR_CW(1,LWT_CW,LN_CW) - well time for time point for well, s
!     VAR_CW(2,LWT_CW,LN_CW) - mass rate for time point for well, kg/s
!     VAR_CW(3,LWT_CW,LN_CW) - press. limit for time point for well, Pa
!     VAR_CW(4,LWT_CW,LN_CW) - CO2 water conc. for time point for well
!     VAR_CW(4,LWT_CW,LN_CW) - aqu. CO2 conc. for time point for well
!     VAR_CW(5,LWT_CW,LN_CW) - temperature for time point for well
!     VAR_CW(6,LWT_CW,LN_CW) - aqu. salt conc. for time point for well
!     PAR_CW(2,LWI_CW) - well bore radius for well interval, m
!     XTP_CW(2,LWI_CW) - x-transition points for well interval, m
!     YTP_CW(2,LWI_CW) - y-transition points for well interval, m
!     ZTP_CW(2,LWI_CW) - z-transition points for well interval, m
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 March 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
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
      CHARACTER*6 FORM1
      CHARACTER*64 ADUM,BDUM(3),UNTS
      CHARACTER*512 CHDUM
      INTEGER NCH
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDCOUP_WELL'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Coupled Well Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of wells  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Coupled Wells'
      CALL RDINT(ISTART,ICOMMA,CHDUM,N_CW)
!
!---  Check number of wells parameter  ---
!
      IF( N_CW.GT.LN_CW ) THEN
        INDX = 5
        CHMSG = 'Number of Coupled Wells > LN_CW'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Loop over number of coupled wells  ---
!
      NIT_CW = 0
      DO 400 NCW = 1,N_CW
!
!---    Read well type  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Coupled Well Type'
        CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'co2').NE.0 .AND.
     &    INDEX(ADUM(1:),'injection').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Coupled CO2 Injection Well'
          VARB = 'Water-Vapor Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM(1))
          WNM_CW(NCW) = 'GIW'
          ICX = ICOUNT(NCW)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(WNM_CW(NCW)(4:4+ICX-1),FORM1) NCW
          IF( INDEX(BDUM(1)(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
            IT_CW(NCW) = 3
          ELSEIF( INDEX(BDUM(1)(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1)(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
            IT_CW(NCW) = 2
          ELSE
            WRITE(IWR,'(2X,A)') 'No Water-Vapor'
            IT_CW(NCW) = 1
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'injection').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Coupled Aqueous Injection Well'
          VARB = 'Dissolved CO2 Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM(1))
          WNM_CW(NCW) = 'LIW'
          ICX = ICOUNT(NCW)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(WNM_CW(NCW)(4:4+ICX-1),FORM1) NCW
          IF( INDEX(BDUM(1)(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved CO2 Relative Saturation'
            IT_CW(NCW) = 6
          ELSEIF( INDEX(BDUM(1)(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1)(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved CO2 Mass Fraction'
            IT_CW(NCW) = 5
          ELSE
            WRITE(IWR,'(2X,A)') 'No Dissolved CO2'
            IT_CW(NCW) = 4
          ENDIF
          VARB = 'Dissolved Salt Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM(1))
          IF( INDEX(BDUM(1)(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved Salt Relative Saturation'
            IT_CW(NCW) = IT_CW(NCW) + 10
          ELSEIF( INDEX(BDUM(1)(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1)(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved Salt Mass Fraction'
            IT_CW(NCW) = IT_CW(NCW) + 20
          ELSE
            WRITE(IWR,'(2X,A)') 'No Dissolved Salt'
            IT_CW(NCW) = IT_CW(NCW) + 30
          ENDIF
        ELSEIF( (INDEX(ADUM(1:),'withdrawl').NE.0 .OR.
     &    INDEX(ADUM(1:),'production').NE.0) .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0  ) THEN
          WRITE(IWR,'(2X,A)') 'Coupled Volume Flow Rate Withdrawl Well'
          IT_CW(NCW) = -1
          WNM_CW(NCW) = 'VPW'
          ICX = ICOUNT(NCW)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(WNM_CW(NCW)(4:4+ICX-1),FORM1) NCW
        ELSEIF( (INDEX(ADUM(1:),'withdrawl').NE.0 .OR.
     &    INDEX(ADUM(1:),'production').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0  ) THEN
          WRITE(IWR,'(2X,A)') 'Coupled Mass Flow Rate Withdrawl Well'
          IT_CW(NCW) = -2
          ICX = ICOUNT(NCW)
          WRITE(FORM1(3:3),'(I1)') ICX
          WNM_CW(NCW) = 'MPW'
          WRITE(WNM_CW(NCW)(4:4+ICX-1),FORM1) NCW
        ELSEIF( INDEX(ADUM(1:),'alternat').NE.0 .AND.
     &    INDEX(ADUM(1:),'state').NE.0 .AND.
     &    INDEX(ADUM(1:),'injection').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Coupled Alternating State Injection Well'
          IT_CW(NCW) = 100
          ICX = ICOUNT(NCW)
          WRITE(FORM1(3:3),'(I1)') ICX
          WNM_CW(NCW) = 'AIW'
          WRITE(WNM_CW(NCW)(4:4+ICX-1),FORM1) NCW
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Coupled Well Type: ' // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read x-direction well fraction factor  ---
!
        VARB = 'X-Well Fraction Factor'
        IDFLT = 1
        FF_CW(1,NCW) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,FF_CW(1,NCW))
        WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',FF_CW(1,NCW)
        IF( FF_CW(1,NCW).LT.EPSL ) THEN
          INDX = 16
          CHMSG = 'Zero X-Well Fraction Factor: Well Number'
          IMSG = NCW
          RMSG = FF_CW(1,NCW)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read y-direction well fraction factor  ---
!
        VARB = 'Y-Well Fraction Factor'
        IDFLT = 1
        FF_CW(2,NCW) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,FF_CW(2,NCW))
        WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',FF_CW(2,NCW)
        IF( FF_CW(2,NCW).LT.EPSL ) THEN
          INDX = 16
          CHMSG = 'Zero Y-Well Fraction Factor: Well Number'
          IMSG = NCW
          RMSG = FF_CW(2,NCW)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read z-direction well fraction factor  ---
!
        VARB = 'Z-Well Fraction Factor'
        IDFLT = 1
        FF_CW(3,NCW) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,FF_CW(3,NCW))
        WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',FF_CW(3,NCW)
        IF( FF_CW(3,NCW).LT.EPSL ) THEN
          INDX = 16
          CHMSG = 'Zero Z-Well Fraction Factor: Well Number'
          IMSG = NCW
          RMSG = FF_CW(3,NCW)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Injection/production limits  ---
!
        TML_CW(NCW) = 1.D+20
!
!---    Read total CO2 mass injected limit for alternating state
!       injection wells---
!
        IF( IT_CW(NCW).EQ.100 ) THEN
          VARB = 'Total CO2 Mass Injected Limit'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,TML_CW(NCW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',TML_CW(NCW)
          INDX = 0
          IUNKG = 1
          CALL RDUNIT(UNTS,TML_CW(NCW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',TML_CW(NCW),', kg)'
!
!---    Read total mass injected limit  ---
!
        ELSEIF( IT_CW(NCW).GT.0 ) THEN
          VARB = 'Total Mass Injected Limit'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,TML_CW(NCW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',TML_CW(NCW)
          INDX = 0
          IUNKG = 1
          CALL RDUNIT(UNTS,TML_CW(NCW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',TML_CW(NCW),', kg)'
!
!---    Read total volume withdrawn/produced limit  ---
!
        ELSEIF( IT_CW(NCW).EQ.-1 ) THEN
          VARB = 'Total Volume Withdrawn/Produced Limit'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,TML_CW(NCW))
          TML_CW(NCW) = ABS(TML_CW(NCW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',TML_CW(NCW)
          INDX = 0
          IUNM = 3
          CALL RDUNIT(UNTS,TML_CW(NCW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',TML_CW(NCW),', m^3)'
!
!---    Read total mass withdrawn/produced limit  ---
!
        ELSEIF( IT_CW(NCW).EQ.-2 ) THEN
          VARB = 'Total Mass Withdrawn/Produced Limit'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,TML_CW(NCW))
          TML_CW(NCW) = ABS(TML_CW(NCW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',TML_CW(NCW)
          INDX = 0
          IUNKG = 1
          CALL RDUNIT(UNTS,TML_CW(NCW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',TML_CW(NCW),', kg)'
        ENDIF
!
!---    Check for solutes  ---
!
        IF( IEQC.NE.0 ) THEN
          NC = 0
          DO
            CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
            ISX = ISTART
            ICX = ICOMMA
            IF( INDX.EQ.0 ) EXIT
!
!---        Check for solute names  ---
!
            ADUM(1:) = ' '
            VARB = 'Solute Name'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
            IFIND = 0
            DO NSL = 1,NSOLU
!
!---          Solute name found  ---
!
              IF( SOLUT(NSL).EQ.ADUM ) THEN
                NC = NC + 1
                ISOLU_CW(NC,NCW) = NSL
                IFIND = 1
                EXIT
              ENDIF
            ENDDO
            IF( IFIND.EQ.0 ) EXIT
          ENDDO
          ISTART = ISX
          ICOMMA = ICX
        ENDIF
!
!---    Read well name  ---
!
        CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
        VARB = 'Well Name'
        IF( INDX.EQ.1 ) THEN
          IDFLT = 1
          CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,WNM_CW(NCW) )
        ENDIF
        WRITE(IWR,'(2X,3A)') VARB(1:IVR),': ',WNM_CW(NCW)(1:NCH)
!
!---    Read number of well intervals  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Number of Well Intervals'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NI_CWX)
        NIT_CW = NIT_CW + NI_CWX
!
!---    Check total number of well intervals parameter  ---
!
        IF( NIT_CW.GT.LWI_CW ) THEN
          INDX = 5
          CHMSG = 'Number of Well Intervals > LWI_CW'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Assign well transition pointers  ---
!
        IF( NCW.EQ.1 ) THEN
          ID_CW(1,NCW) = 1
        ELSE
          ID_CW(1,NCW) = ID_CW(2,NCW-1) + 1
        ENDIF
        ID_CW(2,NCW) = ID_CW(1,NCW) + NI_CWX - 1
!
!---    Loop over number of well intervals  ---
!
        DO 100 NICW = ID_CW(1,NCW),ID_CW(2,NCW)
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
!
!---      Read first x-transition point  ---
!
          VARB = 'Interval First X-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,XTP_CW(1,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',XTP_CW(1,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,XTP_CW(1,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',XTP_CW(1,NICW),', m)'
!
!---      Read first y-transition point  ---
!
          VARB = 'Interval First Y-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,YTP_CW(1,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',YTP_CW(1,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,YTP_CW(1,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',YTP_CW(1,NICW),', m)'
!
!---      Cylindrical coordinates with azimuthal symmetry  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            IF( ABS(XTP_CW(1,NICW))/EPSL.GT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Non-Zero Interval First X-Transition Point ' // 
     &          'for Radially Symmetric Domain: XTP_CW ='
              RLMSG = XTP_CW(1,NICW)
              CALL WRMSGS( INDX )
            ENDIF
            IF( ABS(YTP_CW(1,NICW))/EPSL.GT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Non-Zero Interval First Y-Transition Point ' // 
     &          'for Radially Symmetric Domain: YTP_CW ='
              RLMSG = YTP_CW(1,NICW)
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
!
!---      Read first z-transition point  ---
!
          VARB = 'Interval First Z-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ZTP_CW(1,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',ZTP_CW(1,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,ZTP_CW(1,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',ZTP_CW(1,NICW),', m)'
!
!---      Read second x-transition point  ---
!
          VARB = 'Interval Second X-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,XTP_CW(2,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',XTP_CW(2,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,XTP_CW(2,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',XTP_CW(2,NICW),', m)'
!
!---      Read second y-transition point  ---
!
          VARB = 'Interval Second Y-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,YTP_CW(2,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',YTP_CW(2,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,YTP_CW(2,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',YTP_CW(2,NICW),', m)'
!
!---      Cylindrical coordinates with azimuthal symmetry  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            IF( ABS(XTP_CW(2,NICW))/EPSL.GT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Non-Zero Interval Second X-Transition Point ' // 
     &          'for Radially Symmetric Domain: XTP_CW ='
              RLMSG = XTP_CW(2,NICW)
              CALL WRMSGS( INDX )
            ENDIF
            IF( ABS(YTP_CW(2,NICW))/EPSL.GT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Non-Zero Interval Second Y-Transition Point ' // 
     &          'for Radially Symmetric Domain: YTP_CW ='
              RLMSG = YTP_CW(2,NICW)
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
!
!---      Read second z-transition point  ---
!
          VARB = 'Interval Second Z-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ZTP_CW(2,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',ZTP_CW(2,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,ZTP_CW(2,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',ZTP_CW(2,NICW),', m)'
!
!---      Read well-bore radius  ---
!
          VARB = 'Interval Well-Bore Radius'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(2,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',PAR_CW(2,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,PAR_CW(2,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_CW(2,NICW),', m)'
!
!---      Read well skin factor  ---
!
          VARB = 'Interval Skin Factor'
          IDFLT = 1
          PAR_CW(1,NICW) = 0.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(1,NICW))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',PAR_CW(1,NICW)
!
!---      Read well interval screen option  ---
!
          VARB = 'Interval Screen Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,ADUM)
          IF( INDEX(ADUM(1:),'screened').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Screened Well Interval'
            IS_CW(NICW) = 1
          ELSE
            WRITE(IWR,'(2X,A)') 'Unscreened Well Interval'
            IS_CW(NICW) = 0
          ENDIF
  100   CONTINUE
!
!---    Read number of well time points or well time periods  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        IF( IT_CW(NCW).LT.100 ) THEN
          VARB = 'Number of Well Time Points'
        ELSE
          VARB = 'Number of Well Time Periods'
        ENDIF
        CALL RDINT(ISTART,ICOMMA,CHDUM,IM_CW(NCW))
!
!---    Check number of well time points parameter  ---
!
        IF( IM_CW(NCW).GT.LWT_CW ) THEN
          INDX = 5
          IF( IT_CW(NCW).LT.100 ) THEN
            CHMSG = 'Number of Well Time Points > LWT_CW'
          ELSE
            CHMSG = 'Number of Well Time Periods > LWT_CW'
          ENDIF
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for cyclic well times  ---
!
        VARB = 'Cyclic well time option'
        CALL CHKCHR(ISTART,ICOMMA,CHDUM,INDX)
        IF( INDX.EQ.1 ) THEN
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,ADUM)
          IF( INDEX(ADUM(1:),'cyclic').NE.0 ) THEN
            ICC_CW(NCW) = 1
            WRITE(IWR,'(2X,A)') 'Cyclic Well Times'
          ELSE
            ICC_CW(NCW) = 0
            WRITE(IWR,'(2X,A)') 'Noncyclic Well Times'
          ENDIF
        ENDIF
!
!---    Loop over the number of well time periods/points  ---
!
        NC = 0
        DO 300 M = 1,IM_CW(NCW)
!
!---      Alternating state injection well, read number of time
!         points for each time period, and the state of the injection
!         well for the time period  ---
!
          IF( IT_CW(NCW).GE.100 ) THEN
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
            ISTART = 1
            VARB = 'Number of Time Points in Time Period'
            CALL RDINT(ISTART,ICOMMA,CHDUM,IMP_CW(M,NCW))
            VARB = 'Coupled Injection Well State'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM(1))
!
!---        CO2 injection well state  ---
!
            IF( INDEX(BDUM(1)(1:),'co2').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Coupled Well State: CO2 Injection'
              VARB = 'Water-Vapor Option'
              CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM(2))
              IF( INDEX(BDUM(2)(1:),'rel').NE.0 ) THEN
                WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
                ITS_CW(M,NCW) = 3
              ELSEIF( INDEX(BDUM(2)(1:),'mass').NE.0 .AND.
     &          INDEX(BDUM(2)(1:),'frac').NE.0 ) THEN
                WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
                ITS_CW(M,NCW) = 2
              ELSE
                WRITE(IWR,'(2X,A)') 'No Water-Vapor'
                ITS_CW(M,NCW) = 1
              ENDIF
!
!---        Aqueous injection well state  ---
!
            ELSEIF( INDEX(BDUM(1)(1:),'aqueous').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Coupled Well State: Aqu. Injection'
              VARB = 'Dissolved CO2 Option'
              CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM(2))
              ICX = ICOUNT(NCW)
              WRITE(FORM1(3:3),'(I1)') ICX
              WRITE(WNM_CW(NCW)(4:4+ICX-1),FORM1) NCW
              IF( INDEX(BDUM(2)(1:),'rel').NE.0 ) THEN
                WRITE(IWR,'(2X,A)') 'Dissolved CO2 Relative Saturation'
                ITS_CW(M,NCW) = 6
              ELSEIF( INDEX(BDUM(2)(1:),'mass').NE.0 .AND.
     &          INDEX(BDUM(2)(1:),'frac').NE.0 ) THEN
                WRITE(IWR,'(2X,A)') 'Dissolved CO2 Mass Fraction'
                ITS_CW(M,NCW) = 5
              ELSE
                WRITE(IWR,'(2X,A)') 'No Dissolved CO2'
                ITS_CW(M,NCW) = 4
              ENDIF
              VARB = 'Dissolved Salt Option'
              CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM(3))
              IF( INDEX(BDUM(3)(1:),'rel').NE.0 ) THEN
                WRITE(IWR,'(2X,A)') 'Dissolved Salt Relative Saturation'
                ITS_CW(M,NCW) = ITS_CW(M,NCW) + 30
              ELSEIF( INDEX(BDUM(3)(1:),'mass').NE.0 .AND.
     &          INDEX(BDUM(3)(1:),'frac').NE.0 ) THEN
                WRITE(IWR,'(2X,A)') 'Dissolved Salt Mass Fraction'
                ITS_CW(M,NCW) = ITS_CW(M,NCW) + 20
              ELSE
                WRITE(IWR,'(2X,A)') 'No Dissolved Salt'
                ITS_CW(M,NCW) = ITS_CW(M,NCW) + 10
              ENDIF
            ENDIF              
!
!---        Loop over the number of time points in time period  ---
!
            DO 200 MX = 1,IMP_CW(M,NCW)
              NC = NC + 1
!
!---          Check number of well time points parameter  ---
!
              IF( NC.GT.LWT_CW ) THEN
                INDX = 5
                CHMSG = 'Number of Well Time Points > LWT_CW'
                CALL WRMSGS( INDX )
              ENDIF
!
!---          Read new line  ---
!
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              ISTART = 1
!
!---          Read well time  ---
!
              VARB = 'Well Time'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(1,NC,NCW))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_CW(1,NC,NCW)
              INDX = 0
              IUNS = 1
              CALL RDUNIT(UNTS,VAR_CW(1,NC,NCW),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(1,NC,NCW),', s)'
              IF( NC.GT.1 ) THEN
                IF( VAR_CW(1,NC,NCW).LT.VAR_CW(1,NC-1,NCW) ) THEN
                  INDX = 7
                  CHMSG = 'Nonascending Well Times: Well'
                  IMSG = NCW
                  CALL WRMSGS( INDX )
                ENDIF
              ENDIF
!
!---          Read injection mass rate  ---
!
              IF( ITS_CW(M,NCW).GT.0 ) THEN
                VARB = 'Injection Mass Rate'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(2,NC,NCW))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',VAR_CW(2,NC,NCW)
                IUNKG = 1
                IUNS = -1
                CALL RDUNIT(UNTS,VAR_CW(2,NC,NCW),INDX)
                WRITE(IWR,'(A,1PE11.4,A)') 
     &            ' (',VAR_CW(2,NC,NCW),', kg/s)'
!
!---            Read maximum well-top pressure  ---
!
                VARB = 'Maximum Well-Top Pressure'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(3,NC,NCW))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',VAR_CW(3,NC,NCW)
                INDX = 0
                IUNM = -1
                IUNKG = 1
                IUNS = -2
                CALL RDUNIT(UNTS,VAR_CW(3,NC,NCW),INDX)
                WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(3,NC,NCW),', Pa)'
!
!---            Read water concentration dissolved in CO2  ---
!
                IF( ITS_CW(M,NCW).EQ.3 ) THEN
                  VARB = 'Water Relative Humidity'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(4,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(4,NC,NCW)
                ELSEIF( ITS_CW(M,NCW).EQ.2 ) THEN
                  VARB = 'Water Mass Fraction'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(4,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(4,NC,NCW)
                ELSEIF( ITS_CW(M,NCW).EQ.1 ) THEN
                  VAR_CW(4,NC,NCW) = 0.D+0
                ENDIF
!
!---            Read CO2 concentration dissolved in aqueous  ---
!
                IT_CWX = MOD(ITS_CW(M,NCW),10)
                IF( IT_CWX.EQ.4 ) THEN
                  VAR_CW(4,NC,NCW) = 0.D+0
                ELSEIF( IT_CWX.EQ.5 ) THEN
                  VARB = 'Dissolved CO2 Mass Fraction'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(4,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(4,NC,NCW)
                ELSEIF( IT_CWX.EQ.6 ) THEN
                  VARB = 'Dissolved CO2 Relative Saturation'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(4,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(4,NC,NCW)
                ENDIF
!
!---            Read salt concentration dissolved in aqueous  ---
!
                IT_CWX = ITS_CW(M,NCW)-IT_CWX
                IF( IT_CWX.EQ.10 ) THEN
                  VAR_CW(6,NC,NCW) = 0.D+0
                ELSEIF( IT_CWX.EQ.20 ) THEN
                  VARB = 'Dissolved Salt Mass Fraction'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(6,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(6,NC,NCW)
                ELSEIF( IT_CWX.EQ.30 ) THEN
                  VARB = 'Dissolved Salt Relative Saturation'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(6,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(6,NC,NCW)
                ENDIF
              ENDIF
!
!---          Simulation with solute transport  ---
!
              IF( IEQC.NE.0 ) THEN
                DO NSL_CW = 1,LSOLU_CW
!
!---              Solute active in well  ---
!
                  IF( ISOLU_CW(NSL_CW,NCW).EQ.0 ) EXIT
                  NSL = ISOLU_CW(NSL_CW,NCW)
!
!---              Read solute concentration in well fluid  ---
!
                  VARB = 'Solute Concentration in Well Fluid'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VARC_CW(NSL,NC,NCW))
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &              UNTS(1:NCH),': ',VARC_CW(NSL,NC,NCW)
                  INDX = 0
                  IUNM = -3
                  CALL RDUNIT(UNTS,VARC_CW(NSL,NC,NCW),INDX)
                  WRITE(IWR,'(A,1PE11.4,A)') ' (',VARC_CW(NSL,NC,NCW),
     &              ', 1/m^3)'
                ENDDO
              ENDIF
  200       CONTINUE
!
!---      Single state injection well, read one time
!         point for each time period, and the state of the injection
!         well for the time period  ---
!
          ELSE
            IMP_CW(M,NCW) = 1
!
!---        Read new line  ---
!
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
            ISTART = 1
!
!---        Read well time  ---
!
            VARB = 'Well Time'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(1,M,NCW))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR_CW(1,M,NCW)
            INDX = 0
            IUNS = 1
            CALL RDUNIT(UNTS,VAR_CW(1,M,NCW),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(1,M,NCW),', s)'
!
!---        Read injection mass rate  ---
!
            IF( IT_CW(NCW).GT.0 ) THEN
              VARB = 'Injection Mass Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(2,M,NCW))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_CW(2,M,NCW)
              IUNKG = 1
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_CW(2,M,NCW),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(2,M,NCW),', kg/s)'
!
!---          Read maximum well-top pressure  ---
!
              VARB = 'Maximum Well-Top Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(3,M,NCW))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_CW(3,M,NCW)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_CW(3,M,NCW),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(3,M,NCW),', Pa)'
!
!---        Read withdrawl/production volume flow rate  ---
!
            ELSEIF( IT_CW(NCW).EQ.-1 ) THEN
              VARB = 'Withdrawl/Production Volume Flow Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(2,M,NCW))
              VAR_CW(2,M,NCW) = ABS(VAR_CW(2,M,NCW))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_CW(2,M,NCW)
              IUNM = 3
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_CW(2,M,NCW),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(2,M,NCW),', m^3/s)'
!
!---          Read minimum well-bottom pressure  ---
!
              VARB = 'Minimum Well-Bottom Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(3,M,NCW))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_CW(3,M,NCW)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_CW(3,M,NCW),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(3,M,NCW),', Pa)'
!
!---        Read withdrawl/production mass flow rate  ---
!
            ELSEIF( IT_CW(NCW).EQ.-2 ) THEN
              VARB = 'Withdrawl/Production Mass Flow Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(2,M,NCW))
              VAR_CW(2,M,NCW) = ABS(VAR_CW(2,M,NCW))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_CW(2,M,NCW)
              IUNKG = 1
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_CW(2,M,NCW),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(2,M,NCW),', kg/s)'
!
!---          Read minimum well-bottom pressure  ---
!
              VARB = 'Minimum Well-Bottom Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(3,M,NCW))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_CW(3,M,NCW)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_CW(3,M,NCW),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(3,M,NCW),', Pa)'
            ENDIF
!
!---        Read water concentration dissolved in CO2  ---
!
            IF( IT_CW(NCW).EQ.3 ) THEN
              VARB = 'Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(4,M,NCW))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_CW(4,M,NCW)
            ELSEIF( IT_CW(NCW).EQ.2 ) THEN
              VARB = 'Water Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(4,M,NCW))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_CW(4,M,NCW)
            ELSEIF( IT_CW(NCW).EQ.1 ) THEN
              VAR_CW(4,M,NCW) = 0.D+0
            ENDIF
!
!---        Read CO2 concentration dissolved in aqueous  ---
!
            IT_CWX = MOD(IT_CW(NCW),10)
            IF( IT_CWX.EQ.4 ) THEN
              VARB = 'Dissolved CO2 Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(4,M,NCW))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_CW(4,M,NCW)
            ELSEIF( IT_CWX.EQ.5 ) THEN
              VARB = 'Dissolved CO2 Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(4,M,NCW))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_CW(4,M,NCW)
            ELSEIF( IT_CWX.EQ.6 ) THEN
              VAR_CW(4,M,NCW) = 0.D+0
            ENDIF
!
!---        Read salt concentration dissolved in aqueous  ---
!
            IT_CWX = IT_CW(NCW)-IT_CWX
            IF( IT_CWX.EQ.10 ) THEN
              VARB = 'Dissolved Salt Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(6,M,NCW))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_CW(6,M,NCW)
            ELSEIF( IT_CWX.EQ.20 ) THEN
              VARB = 'Dissolved Salt Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(6,M,NCW))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_CW(6,M,NCW)
            ELSEIF( IT_CWX.EQ.30 ) THEN
              VAR_CW(6,M,NCW) = 0.D+0
            ENDIF
!
!---        Simulation with solute transport  ---
!
            IF( IEQC.NE.0 ) THEN
              DO NSL_CW = 1,LSOLU_CW
!
!---            Solute active in well  ---
!
                IF( ISOLU_CW(NSL_CW,NCW).EQ.0 ) EXIT
                NSL = ISOLU_CW(NSL_CW,NCW)
!
!---            Read solute concentration in well fluid  ---
!
                VARB = 'Solute Concentration in Well Fluid'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VARC_CW(NSL,M,NCW))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',VARC_CW(NSL,M,NCW)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VARC_CW(NSL,M,NCW),INDX)
                WRITE(IWR,'(A,1PE11.4,A)') ' (',VARC_CW(NSL,M,NCW),
     &            ', 1/m^3)'
              ENDDO
            ENDIF
          ENDIF
  300   CONTINUE
  400 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDCOUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDLEAK_WELL
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
!     STOMP-CO2
!
!     Reads the leaky well card.
!     
!     FF_LW(1,LN_LW) - x-direction well fraction factor for well
!     FF_LW(2,LN_LW) - y-direction well fraction factor for well
!     FF_LW(3,LN_LW) - z-direction well fraction factor for well
!     ID_LW(1,LN_LW) - starting well interval index for well
!     ID_LW(2,LN_LW) - ending well interval index for well
!     ID_LW(3,LN_LW) - starting well node index for well
!     ID_LW(4,LN_LW) - ending well node index for well
!     ID_LW(5,LN_LW) - principal well node index for well
!     IS_LW(LWI_LW) - type index for well interval
!       0 - Cased
!       1 - Screened
!     IT_LW(LN_LW) - type index for well
!       1 - Grouted borehole
!     IZ_LW(LWI_LW) - rock/soil/grout/fill type for well
!     N_LW - number of leaky wells
!     PAR_LW(1,LWI_LW) - skin factor for well interval
!     PAR_LW(2,LWI_LW) - well bore radius for well interval, m
!     PAR_LW(3,LWI_LW) - nominal well node spacing, m
!     XTP_LW(2,LWI_LW) - x-transition points for well interval, m
!     YTP_LW(2,LWI_LW) - y-transition points for well interval, m
!     ZTP_LW(2,LWI_LW) - z-transition points for well interval, m
!     WN_LW(LN_LW) - name of leaky well
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 October 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE LEAK_WELL
      USE GRID
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
      CHARACTER*6 FORM1
      CHARACTER*64 ADUM,UNTS
      CHARACTER*512 CHDUM
      INTEGER NCH
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDLEAK_WELL'
      NROCKX = NROCK+1
!
!---  Write card information to ouput file  ---
!
      CARD = 'Leaky Well Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of leaky wells  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Leaky Wells'
      CALL RDINT(ISTART,ICOMMA,CHDUM,N_LW)
!
!---  Check number of leaky wells parameter  ---
!
      IF( N_LW.GT.LN_LW ) THEN
        INDX = 5
        CHMSG = 'Number of Leaky Wells > LN_LW'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Loop over number of leaky wells  ---
!
      NIT_LW = 0
      DO 400 NCW = 1,N_LW
!
!---    Read well type  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Leaky Well Type'
        CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'grouted').NE.0 .OR.
     &    INDEX(ADUM(1:),'filled').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Grouted or Filled Leaky Well'
          WNM_LW(NCW) = 'GLW'
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Coupled Well Type: ' // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read x-direction well fraction factor  ---
!
        VARB = 'X-Well Fraction Factor'
        IDFLT = 1
        FF_LW(1,NCW) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,FF_LW(1,NCW))
        WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',FF_LW(1,NCW)
        IF( FF_LW(1,NCW).LT.EPSL ) THEN
          INDX = 16
          CHMSG = 'Zero X-Well Fraction Factor: Well Number'
          IMSG = NCW
          RMSG = FF_LW(1,NCW)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read y-direction well fraction factor  ---
!
        VARB = 'Y-Well Fraction Factor'
        IDFLT = 1
        FF_LW(2,NCW) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,FF_LW(2,NCW))
        WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',FF_LW(2,NCW)
        IF( FF_LW(2,NCW).LT.EPSL ) THEN
          INDX = 16
          CHMSG = 'Zero Y-Well Fraction Factor: Well Number'
          IMSG = NCW
          RMSG = FF_LW(2,NCW)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read z-direction well fraction factor  ---
!
        VARB = 'Z-Well Fraction Factor'
        IDFLT = 1
        FF_LW(3,NCW) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,FF_LW(3,NCW))
        WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',FF_LW(3,NCW)
        IF( FF_LW(3,NCW).LT.EPSL ) THEN
          INDX = 16
          CHMSG = 'Zero Z-Well Fraction Factor: Well Number'
          IMSG = NCW
          RMSG = FF_LW(3,NCW)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read well name  ---
!
        CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
        VARB = 'Leaky Well Name'
        IF( INDX.EQ.1 ) THEN
          IDFLT = 1
          CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,WNM_LW(NCW) )
        ENDIF
        WRITE(IWR,'(2X,3A)') VARB(1:IVR),': ',WNM_LW(NCW)(1:NCH)
!
!---    Read number of well intervals  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Number of Leaky Well Intervals'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NI_LWX)
        NIT_LW = NIT_LW + NI_LWX
!
!---    Check total number of well intervals parameter  ---
!
        IF( NIT_LW.GT.LWI_LW ) THEN
          INDX = 5
          CHMSG = 'Number of Well Intervals > LWI_LW'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Assign well transition pointers  ---
!
        IF( NCW.EQ.1 ) THEN
          ID_LW(1,NCW) = 1
        ELSE
          ID_LW(1,NCW) = ID_LW(2,NCW-1) + 1
        ENDIF
        ID_LW(2,NCW) = ID_LW(1,NCW) + NI_LWX - 1
!
!---    Loop over number of well intervals  ---
!
        DO 300 NICW = ID_LW(1,NCW),ID_LW(2,NCW)
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
!
!---      Read first x-transition point  ---
!
          VARB = 'Interval First X-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,XTP_LW(1,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',XTP_LW(1,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,XTP_LW(1,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',XTP_LW(1,NICW),', m)'
!
!---      Read first y-transition point  ---
!
          VARB = 'Interval First Y-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,YTP_LW(1,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',YTP_LW(1,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,YTP_LW(1,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',YTP_LW(1,NICW),', m)'
!
!---      Cylindrical coordinates with azimuthal symmetry  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            IF( ABS(XTP_LW(1,NICW))/EPSL.GT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Non-Zero Interval First X-Transition Point ' // 
     &          'for Radially Symmetric Domain: XTP_LW ='
              RLMSG = XTP_LW(1,NICW)
              CALL WRMSGS( INDX )
            ENDIF
            IF( ABS(YTP_LW(1,NICW))/EPSL.GT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Non-Zero Interval First Y-Transition Point ' // 
     &          'for Radially Symmetric Domain: YTP_LW ='
              RLMSG = YTP_LW(1,NICW)
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
!
!---      Read first z-transition point  ---
!
          VARB = 'Interval First Z-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ZTP_LW(1,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',ZTP_LW(1,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,ZTP_LW(1,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',ZTP_LW(1,NICW),', m)'
!
!---      Read second x-transition point  ---
!
          VARB = 'Interval Second X-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,XTP_LW(2,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',XTP_LW(2,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,XTP_LW(2,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',XTP_LW(2,NICW),', m)'
!
!---      Read second y-transition point  ---
!
          VARB = 'Interval Second Y-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,YTP_LW(2,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',YTP_LW(2,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,YTP_LW(2,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',YTP_LW(2,NICW),', m)'
!
!---      Cylindrical coordinates with azimuthal symmetry  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            IF( ABS(XTP_LW(2,NICW))/EPSL.GT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Non-Zero Interval Second X-Transition Point ' // 
     &          'for Radially Symmetric Domain: XTP_LW ='
              RLMSG = XTP_LW(2,NICW)
              CALL WRMSGS( INDX )
            ENDIF
            IF( ABS(YTP_LW(2,NICW))/EPSL.GT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Non-Zero Interval Second Y-Transition Point ' // 
     &          'for Radially Symmetric Domain: YTP_LW ='
              RLMSG = YTP_LW(2,NICW)
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
!
!---      Read second z-transition point  ---
!
          VARB = 'Interval Second Z-Transition Point'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ZTP_LW(2,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',ZTP_LW(2,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,ZTP_LW(2,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',ZTP_LW(2,NICW),', m)'
!
!---      Read interval well-bore radius  ---
!
          VARB = 'Interval Well-Bore Radius'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_LW(2,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',PAR_LW(2,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,PAR_LW(2,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_LW(2,NICW),', m)'
!
!---      Read interval well skin factor  ---
!
          VARB = 'Interval Skin Factor'
          IDFLT = 1
          PAR_LW(1,NICW) = 0.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_LW(1,NICW))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',PAR_LW(1,NICW)
!
!---      Read interval nominal well node spacing  ---
!
          VARB = 'Interval Well Node Spacing'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_LW(3,NICW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',PAR_LW(3,NICW)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,PAR_LW(3,NICW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_LW(3,NICW),', m)'
!
!---      Read well interval screen option  ---
!
          VARB = 'Interval Screen Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,ADUM)
          IF( INDEX(ADUM(1:),'screened').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Screened Well Interval'
            IS_LW(NICW) = 1
          ELSE
            WRITE(IWR,'(2X,A)') 'Cased Well Interval'
            IS_LW(NICW) = 0
          ENDIF
!
!---      Read well interval rock/soil/grout/fill  ---
!
          VARB = 'Interval Rock/Soil Name'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,ADUM)
          IF( INDEX(ROCK(1),'indexing').NE.0 ) THEN
            DO 210 M = NROCKX,NROCK
              IF( ROCK(M).EQ.ADUM ) THEN
                IROCK = M
                GOTO 220
              ENDIF
  210       CONTINUE
            NROCK = NROCK+1
            IF( NROCK.GT.LRC ) THEN
              INDX = 5
              CHMSG = 'Number of Rock/Soil Types > Parameter LRC'
              CALL WRMSGS( INDX )
            ENDIF
            ROCK(NROCK) = ADUM
            IROCK = NROCK
  220       CONTINUE
            IZ_LW(NICW) = IROCK
          ELSE
            DO 230 M = 1,NROCK
              IF( ROCK(M).EQ.ADUM ) THEN
                IROCK = M
                GOTO 240
              ENDIF
  230       CONTINUE
            NROCK = NROCK+1
            IF( NROCK.GT.LRC ) THEN
              INDX = 5
              CHMSG = 'Number of Rock/Soil Types > Parameter LRC'
              CALL WRMSGS( INDX )
            ENDIF
            ROCK(NROCK) = ADUM
            IROCK = NROCK
  240       CONTINUE
            IZ_LW(NICW) = IROCK
          ENDIF
  300   CONTINUE
  400 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDLEAK_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL_COUP_WELL
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
!     STOMP-CO2
!
!     Coupled-well equation residuals
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 21 April 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE FILES
      USE COUP_WELL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL_COUP_WELL'
!
!---  Zero maximum residuals  ---
!
      RSD_CW = 0.D+0
      NSD_CW = 0
!
!---  Loop over coupled wells ---
!
      DO 100 NCW = 1,N_CW
        MP = JM_CW(NCW)        
        DP_CWX = BLU(MP)
!
!---    Injection well ---
!
        IF( IT_CW(NCW).GT.0 ) THEN
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) THEN
            RSDX = 0.D+0
!
!---      Flow controlled coupled well  ---
!
          ELSE
!            IF( ABS(QM_CW(7,NCW)).GT.EPSL ) THEN
!              RSDFX = ABS( RS_CW(1,NCW)/QM_CW(7,NCW) )
!            ELSE
              RSDFX = 0.D+0
!            ENDIF
            RSDPX = ABS(DP_CWX)/(ABS(P_CW(2,NCW))+PATM)
            RSDX = MAX( RSDPX,RSDFX )
          ENDIF
!
!---    Withdrawl well ---
!
        ELSEIF( IT_CW(NCW).LT.0 ) THEN
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) THEN
            RSDX = 0.D+0
!
!---      Flow controlled coupled well  ---
!
          ELSE
!            IF( ABS(QM_CW(7,NCW)).GT.EPSL ) THEN
!              RSDFX = ABS( RS_CW(1,NCW)/QM_CW(7,NCW) )
!            ELSE
              RSDFX = 0.D+0
!            ENDIF
            RSDPX = ABS(DP_CWX)/(ABS(P_CW(2,NCW))+PATM)
            RSDX = MAX( RSDPX,RSDFX )
          ENDIF
        ENDIF
        IF( RSDX.GT.RSD_CW ) THEN
          RSD_CW = RSDX
          NSD_CW = NCW
        ENDIF
  100 CONTINUE
      IF( RSD_CW.GT.RSDMX ) ICNV = 2
!
!---  Unconverged solution Newton-Raphson iteration limit exceeded  ---
!
      IF( ICNV.EQ.2 .AND. NITER.GE.NRIMX ) THEN
        IF( RSD_CW.GE.1.D+2 ) THEN
          WRITE(ISC,'(10X,A)') '---  Excessive Residual  ---'
          WRITE(IWR,'(10X,A)') '---  Excessive Residual  ---'
        ELSE
          WRITE(ISC,'(10X,A)') '---  Convergence Failure  ---'
          WRITE(IWR,'(10X,A)') '---  Convergence Failure  ---'
        ENDIF
        WRITE(ISC,'(4X,A,1PE11.4,A,I6)')
     &    'Coupled Well Maximum Residual = ',RSD_CW,
     &    ': Coupled Well Number = ',NSD_CW
        WRITE(IWR,'(4X,A,1PE11.4,A,I6)')
     &    'Coupled Well Maximum Residual = ',RSD_CW,
     &    ': Coupled Well Number = ',NSD_CW
!!
!!---    Reduce time step  ---
!!
!        IF( NTSR.LT.4 .OR. (DTCF*DT).GT.DTMN ) THEN
!          NTSR = NTSR + 1
!          DTX = DT
!          TM = TM - (1.D+0-DTCF)*DT
!          DT = DTCF*DT
!          DTO = DT
!          DTI = 1.D+0/DT
!          VAR = DT
!          VARX = DTX
!          IF( UNTM.NE.'null' ) THEN
!            INDX = 1
!            IUNS = 1
!            CALL RDUNIT(UNTM,VAR,INDX)
!            IUNS = 1
!            CALL RDUNIT(UNTM,VARX,INDX)
!            NCH = INDEX( UNTM,'  ')-1
!          ENDIF
!          WRITE(ISC,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
!     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
!     &      VAR,UNTM(1:NCH)
!          WRITE(IWR,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
!     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
!     &      VAR,UNTM(1:NCH)
!          DO 400 N = 1,NFLD
!            PL(2,N) = PL(1,N)
!            PG(2,N) = PG(1,N)
!            PVW(2,N) = PVW(1,N)
!            XLA(2,N) = XLA(1,N)
!            SG(2,N) = SG(1,N)
!            SGT(2,N) = SGT(1,N)
!            ASLMIN(2,N) = ASLMIN(1,N)
!            YLS(2,N) = YLS(1,N)
!            TMS(2,N) = TMS(1,N)
!            NPHAZ(2,N) = NPHAZ(1,N)
!  400     CONTINUE
!!
!!---      Coupled-well pressure  ---
!!
!          DO 402 NCW = 1,N_CW
!            P_CW(2,NCW) = P_CW(1,NCW)
!  402     CONTINUE
!!
!!---  Number of time step reductions failure: stop simulation  ---
!!
!        ELSE
!          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit Exceeded
!     & ---'
!          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit Exceeded
!     & ---'
!          ICNV = 4
!        ENDIF
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDL_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT_COUP_WELL( NSL )
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
!     STOMP-EOR
!
!     Compute solute source transport terms for coupled wells.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 12 January 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FDVP
      USE COUP_WELL
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
      SUB_LOG(ISUB_LOG) = '/SORT_COUP_WELL'
!
!---  Loop over coupled wells  ---
!
      L1: DO NCW = 1,N_CW
!
!---    Skip for conservation or kinetic components  ---
!
        IF( NSL.GT.NSOLU ) CYCLE L1
!
!---    Check for solute in well  ---
!
        IFIND = 0
        L2: DO M = 1,NSOLU
          IF( ISOLU_CW(M,NCW).EQ.0 ) EXIT L2
          IF( ISOLU_CW(M,NCW).EQ.NSL ) THEN
            IFIND = 1
            NC = M
            EXIT L2
          ENDIF
        ENDDO L2
!
!---    Solute not found in well, cycle to next well ---
!
        IF( IFIND.EQ.0 ) CYCLE L1
!
!---    Coupled well time interval ---
!
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---    Cyclic time periods  ---
!
        IF( ICC_CW(NCW).EQ.1 ) THEN
!
!---      Loop over the coupled well time periods, 
!         to find the final well time  ---
!
          NTX = 0
          L3: DO NTP = 1,IM_CW(NCW)
            NTX = NTX + IMP_CW(NTP,NCW)
          ENDDO L3
!
!---      Determine time with the cyclic time period  ---
!
          TMZ = MOD( TM,VAR_CW(1,NTX,NCW) )
          IF( TM.GT.VAR_CW(1,NTX,NCW) ) THEN
            IF( TMZ.LT.EPSL ) TMZ = VAR_CW(1,NTX,NCW)
          ENDIF
        ENDIF
!
!---    Coupled well is inactive  ---
!
        IF( TMZ.LE.VAR_CW(1,1,NCW) ) CYCLE L1
!
!---    Loop over the coupled well time periods  ---
!
        NS = 1
        IFIND = 0
        L4: DO NTP = 1,IM_CW(NCW)
!
!---      Coupled well time period only has one time (start time)  ---
!
          IF( IMP_CW(NTP,NCW).EQ.1 ) THEN
!
!---        Time prior to start time, coupled well is inactive,
!           cycle to next well  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) ) CYCLE L1
!
!---        Time after start time, coupled well is active  ---
!
            VARC_CWX = VARC_CW(NSL,1,NCW)
            IFIND = 1
            EXIT L4
!
!---      Coupled well time period has multiple times  ---
!
          ELSE
            NE = NS + IMP_CW(NTP,NCW) - 1
!
!---        Time outside of coupled well time period, go to next 
!           coupled well time period  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) .OR. 
     &        TMZ.GT.VAR_CW(1,NE,NCW) ) THEN
              NS = NS + IMP_CW(NTP,NCW)
              CYCLE L4
            ENDIF
!
!---        Coupled well time period has multiple time points, use  
!           linear interpolation of well parameters between 
!           time points  ---
!
            L5: DO M = 2,IMP_CW(NTP,NCW)
              MX = NS + M - 1
              IF( TMZ.LE.VAR_CW(1,MX,NCW) ) THEN
                TD_CW = VAR_CW(1,MX,NCW)-VAR_CW(1,MX-1,NCW)
                DT_CW = MIN( VAR_CW(1,MX,NCW)-TMZ,DT )
                TF_CW = (TMZ-VAR_CW(1,MX-1,NCW))/TD_CW
                VARC_CWX = VARC_CW(NSL,MX-1,NCW) + 
     &            TF_CW*(VARC_CW(NSL,MX,NCW)-VARC_CW(NSL,MX-1,NCW))
                IFIND = 1
                EXIT L4
              ENDIF
            ENDDO L5
          ENDIF
          NS = NS + IMP_CW(NTP,NCW)
        ENDDO L4
!
!---    Coupled well is inactive, cycle to next well  ---
!
        IF( IFIND.EQ.0 ) CYCLE L1
!
!---    Loop over coupled-well nodes  ---
!
        L6: DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(NWN)
          IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
          MP = IXP(N)
          IF( ILES.EQ.1 ) THEN
            MCOL = MP
            MROW = MDT
          ELSEIF( ILES.EQ.3 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#ifdef lis
          ELSEIF( ILES.EQ.4 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#endif
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#endif
          ENDIF
          SORTX = 0.D+0
!
!---      Injection well (volumetric fluxes are positive from well) ---
!
          IF( IT_CW(NCW).GT.0 ) THEN
            BLU(MP) = BLU(MP) + Q_CW(1,NWN)*VARC_CWX
!
!---      Production well (volumetric fluxes are positive into well) ---
!
          ELSEIF( IT_CW(NCW).LT.0 ) THEN
!
!---        Solute produced via aqueous phase production  ---
!
            CLX = PORD(2,N)*SL(2,N)
            IF( CLX.GT.EPSL ) SORTX = SORTX + Q_CW(2,NWN)/CLX
!
!---        Solute produced via gas phase production  ---
!
            CGX = PORD(2,N)*SG(2,N)
            IF( CGX.GT.EPSL ) SORTX = SORTX + Q_CW(3,NWN)/CGX
          ENDIF
!
!---      Load Jacobian  ---
!
          IF( ILES.EQ.1 ) THEN
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + SORTX
          ELSEIF( ILES.EQ.3 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#ifdef lis
          ELSEIF( ILES.EQ.4 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#endif
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#endif
          ENDIF
        ENDDO L6
      ENDDO L1
#ifdef ecke
      NEQ = NSL - NSOLU
!
!---  Loop over coupled wells ---
!
      L11: DO NCW = 1,N_CW
!
!---    Skip for passive solutes  ---
!
        IF( NSL.LE.NSOLU ) CYCLE L11
!
!---    Skip for non-conservation solutes  ---
!
        IF( NEQ.LT.1 .OR. NEQ.GT.NEQC ) CYCLE L11
!
!---    Check for conservation component solute in well  ---
!
        IF( NEQ.GT.0 .AND. NEQ.LE.NEQC ) THEN
          IF( ISOLC_CW(NEQ,NCW).EQ.0 ) CYCLE L11
        ENDIF
!
!---    Coupled well time interval ---
!
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---    Cyclic time periods  ---
!
        IF( ICC_CW(NCW).EQ.1 ) THEN
!
!---      Loop over the coupled well time periods, 
!         to find the final well time  ---
!
          NTX = 0
          L13: DO NTP = 1,IM_CW(NCW)
            NTX = NTX + IMP_CW(NTP,NCW)
          ENDDO L13
!
!---      Determine time with the cyclic time period  ---
!
          TMZ = MOD( TM,VAR_CW(1,NTX,NCW) )
          IF( TM.GT.VAR_CW(1,NTX,NCW) ) THEN
            IF( TMZ.LT.EPSL ) TMZ = VAR_CW(1,NTX,NCW)
          ENDIF
        ENDIF
!
!---    Coupled well is inactive  ---
!
        IF( TMZ.LE.VAR_CW(1,1,NCW) ) CYCLE L11
!
!---    Loop over the coupled well time periods  ---
!
        NS = 1
        IFIND = 0
        L14: DO NTP = 1,IM_CW(NCW)
!
!---      Coupled well time period only has one time (start time)  ---
!
          IF( IMP_CW(NTP,NCW).EQ.1 ) THEN
            ITS_CWX = MOD( ITS_CW(NTP,NCW),100 )
!
!---        Time prior to start time, coupled well is inactive,
!           cycle to next well  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) ) CYCLE L11
!
!---        Time after start time, coupled well is active  ---
!
            VARC_CWX = 0.D+0
!
!---        Aqueous species only
!
            IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
              DO NSPCW = 1,NSP_CW(NCW)
                NSP = ISPC_CW(NSPCW,NCW)
                IF( NSP.EQ.0 ) EXIT
                IF( NSP.LE.NSPL ) THEN
                  DO IX = 1,IEQ_C(1,NEQ)
                    NSP_C = IEQ_C(IX+1,NEQ)
                    IF( NSP_C.EQ.NSP ) THEN
                      VARC_CWX = VARC_CWX + 
     &                  EQ_C(IX,NEQ)*VARSP_CW(NSPCW,1,NCW)
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF   
              ENDDO             
!
!---        Gas species only
!      
            ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN
              DO NSPCW = 1,NSP_CW(NCW)
                NSP = ISPC_CW(NSPCW,NCW)
                IF( NSP.EQ.0 ) EXIT
                IF( NSP.GT.NSPL+NSPS .AND. 
     &            NSP.LE.NSPL+NSPS+NSPG ) THEN
                  DO IX = 1,IEQ_C(1,NEQ)
                    NSP_C = IEQ_C(IX+1,NEQ)
                    IF( NSP_C.EQ.NSP ) THEN
                      VARC_CWX = VARC_CWX + 
     &                  EQ_C(IX,NEQ)*VARSP_CW(NSPCW,1,NCW)
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF   
              ENDDO             
            ENDIF
            IFIND = 1
            EXIT L14
!
!---      Coupled well time period has multiple times  ---
!
          ELSE
            NE = NS + IMP_CW(NTP,NCW) - 1
            ITS_CWX = MOD( ITS_CW(NTP,NCW),100 )
!
!---        Time outside of coupled well time period, go to next 
!           coupled well time period  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) .OR. 
     &        TMZ.GT.VAR_CW(1,NE,NCW) ) THEN
              NS = NS + IMP_CW(NTP,NCW)
              CYCLE L14
            ENDIF
!
!---        Coupled well time period has multiple time points, use  
!           linear interpolation of well parameters between 
!           time points  ---
!
            L15: DO M = 2,IMP_CW(NTP,NCW)
              MX = NS + M - 1
              IF( TMZ.LE.VAR_CW(1,MX,NCW) ) THEN
                TD_CW = VAR_CW(1,MX,NCW)-VAR_CW(1,MX-1,NCW)
                DT_CW = MIN( VAR_CW(1,MX,NCW)-TMZ,DT )
                TF_CW = (TMZ-VAR_CW(1,MX-1,NCW))/TD_CW
!
!---            Time after start time, coupled well is active  ---
!
                VARC_CWX = 0.D+0
!
!---            Aqueous species only
!
                IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
                  DO NSPCW = 1,NSP_CW(NCW)
                    NSP = ISPC_CW(NSPCW,NCW)
                    IF( NSP.EQ.0 ) EXIT
                    IF( NSP.LE.NSPL ) THEN
                      DO IX = 1,IEQ_C(1,NEQ)
                        NSP_C = IEQ_C(IX+1,NEQ)
                        IF( NSP_C.EQ.NSP ) THEN
                          VARC_CWX = VARC_CWX + EQ_C(IX,NEQ)*
     &                      (VARSP_CW(NSPCW,MX-1,NCW) + TF_CW*
     &                      (VARSP_CW(NSPCW,MX,NCW) - 
     &                      VARSP_CW(NSPCW,MX-1,NCW)))
                          EXIT
                        ENDIF
                      ENDDO
                    ENDIF   
                  ENDDO             
!
!---            Gas species only
!
                ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN
                  DO NSPCW = 1,NSP_CW(NCW)
                    NSP = ISPC_CW(NSPCW,NCW)
                    IF( NSP.EQ.0 ) EXIT
                    IF( NSP.GT.NSPL+NSPS .AND. 
     &                NSP.LE.NSPL+NSPS+NSPG ) THEN
                      DO IX = 1,IEQ_C(1,NEQ)
                        NSP_C = IEQ_C(IX+1,NEQ)
                        IF( NSP_C.EQ.NSP ) THEN
                          VARC_CWX = VARC_CWX + EQ_C(IX,NEQ)*
     &                    (VARSP_CW(NSPCW,MX-1,NCW) + TF_CW*
     &                    (VARSP_CW(NSPCW,MX,NCW) - 
     &                    VARSP_CW(NSPCW,MX-1,NCW)))
                          EXIT
                        ENDIF
                      ENDDO
                    ENDIF   
                  ENDDO             
                ENDIF
                IFIND = 1
                EXIT L14
              ENDIF
            ENDDO L15
          ENDIF
          NS = NS + IMP_CW(NTP,NCW)
        ENDDO L14
!
!---    Coupled well is inactive, cycle to next well  ---
!
        IF( IFIND.EQ.0 ) CYCLE L11
!
!---    Loop over coupled-well nodes  ---
!
        L16: DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(NWN)
          IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
          MP = IXP(N)
          IF( ILES.EQ.1 ) THEN
            MCOL = MP
            MROW = MDT
          ELSEIF( ILES.EQ.3 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#ifdef lis
          ELSEIF( ILES.EQ.4 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#endif
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#endif
          ENDIF
          SORTX = 0.D+0
!
!---      Injection well (volumetric fluxes are positive from well) ---
!
!         Q_CW(1,NWN) - total volumetric flux, m^3/s
!         Q_CW(2,NWN) - aqueous volumetric flux, m^3/s
!         Q_CW(3,NWN) - gas volumetric flux, m^3/s
!
          IF( IT_CW(NCW).GT.0 ) THEN
!
!---        Aqueous species only
!
            IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
              BLU(MP) = BLU(MP) + Q_CW(2,NWN)*VARC_CWX
!
!---        Gas species only
!
            ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN
              BLU(MP) = BLU(MP) + Q_CW(3,NWN)*VARC_CWX
            ENDIF
!
!---      Production well (volumetric fluxes are positive into well) ---
!
          ELSEIF( IT_CW(NCW).LT.0 ) THEN
!
!---        Solute produced via aqueous phase production  ---
!
            CLX = PORD(2,N)*SL(2,N)
            IF( CLX.GT.EPSL ) SORTX = SORTX + Q_CW(2,NWN)/CLX
!
!---        Solute produced via gas phase production  ---
!
            CGX = PORD(2,N)*SG(2,N)
            IF( CGX.GT.EPSL ) SORTX = SORTX + Q_CW(3,NWN)/CGX
          ENDIF
!
!---      Load Jacobian  ---
!
          IF( ILES.EQ.1 ) THEN
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + SORTX
          ELSEIF( ILES.EQ.3 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#ifdef lis
          ELSEIF( ILES.EQ.4 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#endif
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#endif
          ENDIF
        ENDDO L16
      ENDDO L11
      NEQ = NEQ - NEQC
!
!---  Loop over coupled wells ---
!
      L21: DO NCW = 1,N_CW
!
!---    Skip for passive solutes  ---
!
        IF( NSL.LE.NSOLU ) CYCLE L21
!
!---    Skip for non-kinetic solutes  ---
!
        IF( NEQ.LT.1 .OR. NEQ.GT.NEQK ) CYCLE L21
!
!---    Coupled well time interval ---
!
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---    Cyclic time periods  ---
!
        IF( ICC_CW(NCW).EQ.1 ) THEN
!
!---      Loop over the coupled well time periods, 
!         to find the final well time  ---
!
          NTX = 0
          L23: DO NTP = 1,IM_CW(NCW)
            NTX = NTX + IMP_CW(NTP,NCW)
          ENDDO L23
!
!---      Determine time with the cyclic time period  ---
!
          TMZ = MOD( TM,VAR_CW(1,NTX,NCW) )
          IF( TM.GT.VAR_CW(1,NTX,NCW) ) THEN
            IF( TMZ.LT.EPSL ) TMZ = VAR_CW(1,NTX,NCW)
          ENDIF
        ENDIF
!
!---    Coupled well is inactive  ---
!
        IF( TMZ.LE.VAR_CW(1,1,NCW) ) CYCLE L21
!
!---    Loop over the coupled well time periods  ---
!
        NS = 1
        IFIND = 0
        L24: DO NTP = 1,IM_CW(NCW)
!
!---      Coupled well time period only has one time (start time)  ---
!
          IF( IMP_CW(NTP,NCW).EQ.1 ) THEN
            ITS_CWX = MOD( ITS_CW(NTP,NCW),100 )
!
!---        Time prior to start time, coupled well is inactive,
!           cycle to next well  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) ) CYCLE L21
!
!---        Time after start time, coupled well is active  ---
!
            VARC_CWX = 0.D+0
!
!---        Aqueous species only
!
            IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
              DO NSPCW = 1,NSP_CW(NCW)
                NSP = ISPC_CW(NSPCW,NCW)
                IF( NSP.EQ.0 ) EXIT
                IF( NSP.LE.NSPL ) THEN
                  DO IX = 1,IEQ_K(1,NEQ)
                    NSP_K = IEQ_K(IX+1,NEQ)
                    IF( NSP_K.EQ.NSP ) THEN
                      VARC_CWX = VARC_CWX + 
     &                  EQ_K(IX,NEQ)*VARSP_CW(NSPCW,1,NCW)
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF   
              ENDDO             
!
!---        Gas species only
!
            ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN
              DO NSPCW = 1,NSP_CW(NCW)
                NSP = ISPC_CW(NSPCW,NCW)
                IF( NSP.EQ.0 ) EXIT
                IF( NSP.GT.NSPL+NSPS .AND. 
     &            NSP.LE.NSPL+NSPS+NSPG ) THEN
                  DO IX = 1,IEQ_K(1,NEQ)
                    NSP_K = IEQ_K(IX+1,NEQ)
                    IF( NSP_K.EQ.NSP ) THEN
                      VARC_CWX = VARC_CWX + 
     &                  EQ_K(IX,NEQ)*VARSP_CW(NSPCW,1,NCW)
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF   
              ENDDO             
            ENDIF
            IFIND = 1
            EXIT L24
!
!---      Coupled well time period has multiple times  ---
!
          ELSE
            NE = NS + IMP_CW(NTP,NCW) - 1
            ITS_CWX = MOD( ITS_CW(NTP,NCW),100 )
!
!---        Time outside of coupled well time period, go to next 
!           coupled well time period  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) .OR. 
     &        TMZ.GT.VAR_CW(1,NE,NCW) ) THEN
              NS = NS + IMP_CW(NTP,NCW)
              CYCLE L24
            ENDIF
!
!---        Coupled well time period has multiple time points, use  
!           linear interpolation of well parameters between 
!           time points  ---
!
            L25: DO M = 2,IMP_CW(NTP,NCW)
              MX = NS + M - 1
              IF( TMZ.LE.VAR_CW(1,MX,NCW) ) THEN
                TD_CW = VAR_CW(1,MX,NCW)-VAR_CW(1,MX-1,NCW)
                DT_CW = MIN( VAR_CW(1,MX,NCW)-TMZ,DT )
                TF_CW = (TMZ-VAR_CW(1,MX-1,NCW))/TD_CW
!
!---            Time after start time, coupled well is active  ---
!
                VARC_CWX = 0.D+0
!
!---            Aqueous species only
!
                IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
                  DO NSPCW = 1,NSP_CW(NCW)
                    NSP = ISPC_CW(NSPCW,NCW)
                    IF( NSP.EQ.0 ) EXIT
                    IF( NSP.LE.NSPL ) THEN
                      DO IX = 1,IEQ_K(1,NEQ)
                        NSP_K = IEQ_K(IX+1,NEQ)
                        IF( NSP_K.EQ.NSP ) THEN
                          VARC_CWX = VARC_CWX + EQ_K(IX,NEQ)*
     &                      (VARSP_CW(NSPCW,MX-1,NCW) + TF_CW*
     &                      (VARSP_CW(NSPCW,MX,NCW) - 
     &                      VARSP_CW(NSPCW,MX-1,NCW)))
                          EXIT
                        ENDIF
                      ENDDO
                    ENDIF   
                  ENDDO             
!
!---            Gas species only
!
                ELSEIF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
                  DO NSPCW = 1,NSP_CW(NCW)
                    NSP = ISPC_CW(NSPCW,NCW)
                    IF( NSP.EQ.0 ) EXIT
                    IF( NSP.GT.NSPL+NSPS .AND. 
     &                NSP.LE.NSPL+NSPS+NSPG ) THEN
                      DO IX = 1,IEQ_K(1,NEQ)
                        NSP_K = IEQ_K(IX+1,NEQ)
                        IF( NSP_K.EQ.NSP ) THEN
                          VARC_CWX = VARC_CWX + EQ_K(IX,NEQ)*
     &                      (VARSP_CW(NSPCW,MX-1,NCW) + TF_CW*
     &                      (VARSP_CW(NSPCW,MX,NCW) - 
     &                      VARSP_CW(NSPCW,MX-1,NCW)))
                          EXIT
                        ENDIF
                      ENDDO
                    ENDIF   
                  ENDDO             
                ENDIF
                IFIND = 1
                EXIT L24
              ENDIF
            ENDDO L25
          ENDIF
          NS = NS + IMP_CW(NTP,NCW)
        ENDDO L24
!
!---    Coupled well is inactive, cycle to next well  ---
!
        IF( IFIND.EQ.0 ) CYCLE L21
!
!---    Loop over coupled-well nodes  ---
!
        L26: DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(NWN)
          IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
          MP = IXP(N)
          IF( ILES.EQ.1 ) THEN
            MCOL = MP
            MROW = MDT
          ELSEIF( ILES.EQ.3 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#ifdef lis
          ELSEIF( ILES.EQ.4 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#endif
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#endif
          ENDIF
          SORTX = 0.D+0
!
!---      Injection well (volumetric fluxes are positive from well) ---
!
!         Q_CW(1,NWN) - total volumetric flux, m^3/s
!         Q_CW(2,NWN) - aqueous volumetric flux, m^3/s
!         Q_CW(3,NWN) - gas volumetric flux, m^3/s
!
          IF( IT_CW(NCW).GT.0 ) THEN
!
!---        Aqueous species only
!
            IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
              BLU(MP) = BLU(MP) + Q_CW(2,NWN)*VARC_CWX
!
!---        Gas species only
!
            ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN
              BLU(MP) = BLU(MP) + Q_CW(3,NWN)*VARC_CWX
            ENDIF
!
!---      Production well (volumetric fluxes are positive into well) ---
!
          ELSEIF( IT_CW(NCW).LT.0 ) THEN
!
!---        Solute produced via aqueous phase production  ---
!
            CLX = PORD(2,N)*SL(2,N)
            IF( CLX.GT.EPSL ) SORTX = SORTX + Q_CW(2,NWN)/CLX
!
!---        Solute produced via gas phase production  ---
!
            CGX = PORD(2,N)*SG(2,N)
            IF( CGX.GT.EPSL ) SORTX = SORTX + Q_CW(3,NWN)/CGX
          ENDIF
!
!---      Load Jacobian  ---
!
          IF( ILES.EQ.1 ) THEN
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + SORTX
          ELSEIF( ILES.EQ.3 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#ifdef lis
          ELSEIF( ILES.EQ.4 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#endif
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#endif
          ENDIF
        ENDDO L26
      ENDDO L21
#endif
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORT_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_COUP_WELL
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
!     STOMP-CO2
!
!     Update coupled-well pressure.  Injection wells are limited
!     by a high-pressure limit, and withdrawl wells are limited by a 
!     low-pressure limit.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 20 April 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_COUP_WELL'
!
!---  Loop over coupled wells ---
!
      DO 100 NCW = 1,N_CW
        MP = JM_CW(NCW)
        DP_CWX = BLU(MP)
        DPX = 5.D+5
        DP_CWX = SIGN( MIN(ABS(DPX),ABS(DP_CWX)),DP_CWX )
        P_CW(2,NCW) = P_CW(2,NCW) + DP_CWX
!
!---    Limit coupled-well pressure to upper limit for injection
!       wells or lower limit for withdrawl wells  ---
!
        IF( IT_CW(NCW).GT.0 ) THEN
          P_CW(2,NCW) = MIN( PL_CW(NCW),P_CW(2,NCW) )
        ELSEIF( IT_CW(NCW).LT.0 ) THEN
          P_CW(2,NCW) = MAX( PL_CW(NCW),P_CW(2,NCW) )
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WR_WELL
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
!     STOMP-CO2E
!
!     Write well.dat file.
!     
!     ID_CW(1,LN_CW) - starting well interval index for well
!     ID_CW(2,LN_CW) - ending well interval index for well
!     ID_CW(3,LN_CW) - starting well node index for well
!     ID_CW(4,LN_CW) - ending well node index for well
!     ID_CW(5,LN_CW) - principal well node index for well
!
!     IT_CW(LN_CW) - type index for well
!     Well types
!     1 - Injection Mass Rate, kg/s
!     2 - Injection Volume Rate, m^3/s
!     -1 - Fluid Production Constant Bottomhole Pressure, Pa
!     -2 - Liquid Production Constant Bottomhole Pressure, Pa
!
!     N_CW - number of coupled wells
!     XTP_CW(2,LWI_CW) - x-transition points for well interval, m
!     YTP_CW(2,LWI_CW) - y-transition points for well interval, m
!     ZTP_CW(2,LWI_CW) - z-transition points for well interval, m
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
      USE OUTPU
      USE LEAK_WELL
      USE FDVP
      USE COUP_WELL
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
      SUB_LOG(ISUB_LOG) = '/WR_WELL'
!
!---  Open 'well.dat' file for writing  ---
!
      OPEN(UNIT=26, FILE='well.dat', STATUS='UNKNOWN', FORM='FORMATTED')
      CLOSE(UNIT=26,STATUS='DELETE')
      OPEN(UNIT=26, FILE='well.dat', STATUS='NEW', FORM='FORMATTED')
!
!---  Length conversion for output ---
!
      VARX = 1.D+0
      INDX = 4
      IUNM = 1
      CALL RDUNIT(UNLN,VARX,INDX)
!
!---  Loop over number of coupled wells  ---
!
      DO 400 NCW = 1,N_CW
!
!---    Write a purple well name text tag and geometry data to  
!       "well.dat" for an injection well  ---
!
        IF( IT_CW(NCW).GT.0 ) THEN
          NCH = INDEX( WNM_CW(NCW),'  ')-1
          NICW = ID_CW(1,NCW)
          WRITE(26,'(2A,3(1PE12.5,A),2A)') 
     &      'TEXT C=PURPLE HU=POINT H=10 CS=GRID3D',
     &      ', X=',VARX*XTP_CW(1,NICW),
     &      ', Y=',VARX*YTP_CW(1,NICW),
     &      ', Z=',VARX*ZTP_CW(1,NICW),
     &      ', T="',WNM_CW(NCW)(1:NCH),'"'
          WRITE(26,'(A)') 'GEOMETRY T=LINE3D, C=PURPLE LT=0.2'
          WRITE(26,'(A)') '1'
          NC = 0
          DO 110 NICW = ID_CW(1,NCW),ID_CW(2,NCW)
            IF( NICW.GT.ID_CW(1,NCW) ) THEN
              IF( ABS(XTP_CW(1,NICW)-XTP_CW(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(YTP_CW(1,NICW)-YTP_CW(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(ZTP_CW(1,NICW)-ZTP_CW(2,NICW-1)).LT.1.D-3 ) THEN
                NC = NC + 1
              ELSE
                NC = NC + 2
              ENDIF
            ELSE
              NC = NC + 1
            ENDIF
            IF( NICW.EQ.ID_CW(2,NCW) ) THEN
              NC = NC + 1
            ENDIF
  110     CONTINUE
          WRITE(FORM1(3:3),'(I1)') ICOUNT(NC)
          WRITE(26,FORM1) NC
          DO 112 NICW = ID_CW(1,NCW),ID_CW(2,NCW)
            IF( NICW.EQ.ID_CW(1,NCW) ) THEN
              WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_CW(1,NICW),
     &          VARX*YTP_CW(1,NICW),VARX*ZTP_CW(1,NICW)
            ENDIF
            IF( NICW.GT.ID_CW(1,NCW) ) THEN
              IF( ABS(XTP_CW(1,NICW)-XTP_CW(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(YTP_CW(1,NICW)-YTP_CW(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(ZTP_CW(1,NICW)-ZTP_CW(2,NICW-1)).LT.1.D-3 ) THEN
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_CW(1,NICW),
     &            VARX*YTP_CW(1,NICW),VARX*ZTP_CW(1,NICW)
              ELSE
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)') 
     &            VARX*XTP_CW(2,NICW-1),
     &            VARX*YTP_CW(2,NICW-1),VARX*ZTP_CW(2,NICW-1)
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_CW(1,NICW),
     &            VARX*YTP_CW(1,NICW),VARX*ZTP_CW(1,NICW)
              ENDIF
            ENDIF
            IF( NICW.EQ.ID_CW(2,NCW) ) THEN
              WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_CW(2,NICW),
     &          VARX*YTP_CW(2,NICW),VARX*ZTP_CW(2,NICW)
            ENDIF
  112     CONTINUE
!
!---    Write a cyan well name text tag and geometry data to  
!       "well.dat" for a production well  ---
!
        ELSEIF( IT_CW(NCW).LT.0 ) THEN
          NCH = INDEX( WNM_CW(NCW),'  ')-1
          NICW = ID_CW(1,NCW)
          WRITE(26,'(2A,3(1PE12.5,A),2A)') 
     &      'TEXT C=CYAN HU=POINT H=10 CS=GRID3D',
     &      ', X=',VARX*XTP_CW(1,NICW),
     &      ', Y=',VARX*YTP_CW(1,NICW),
     &      ', Z=',VARX*ZTP_CW(1,NICW),
     &      ', T="',WNM_CW(NCW)(1:NCH),'"'
          WRITE(26,'(A)')  'GEOMETRY T=LINE3D, C=CYAN LT=0.2'
          WRITE(26,'(A)') '1'
          NC = 0
          DO 120 NICW = ID_CW(1,NCW),ID_CW(2,NCW)
            IF( NICW.GT.ID_CW(1,NCW) ) THEN
              IF( ABS(XTP_CW(1,NICW)-XTP_CW(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(YTP_CW(1,NICW)-YTP_CW(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(ZTP_CW(1,NICW)-ZTP_CW(2,NICW-1)).LT.1.D-3 ) THEN
                NC = NC + 1
              ELSE
                NC = NC + 2
              ENDIF
            ELSE
              NC = NC + 1
            ENDIF
            IF( NICW.EQ.ID_CW(2,NCW) ) THEN
              NC = NC + 1
            ENDIF
  120     CONTINUE
          WRITE(FORM1(3:3),'(I1)') ICOUNT(NC)
          WRITE(26,FORM1) NC
          DO 122 NICW = ID_CW(1,NCW),ID_CW(2,NCW)
            IF( NICW.EQ.ID_CW(1,NCW) ) THEN
              WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_CW(1,NICW),
     &          VARX*YTP_CW(1,NICW),VARX*ZTP_CW(1,NICW)
            ENDIF
            IF( NICW.GT.ID_CW(1,NCW) ) THEN
              IF( ABS(XTP_CW(1,NICW)-XTP_CW(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(YTP_CW(1,NICW)-YTP_CW(2,NICW-1)).LT.1.D-3 .AND.
     &          ABS(ZTP_CW(1,NICW)-ZTP_CW(2,NICW-1)).LT.1.D-3 ) THEN
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_CW(1,NICW),
     &            VARX*YTP_CW(1,NICW),VARX*ZTP_CW(1,NICW)
              ELSE
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)')
     &            VARX*XTP_CW(2,NICW-1),
     &            VARX*YTP_CW(2,NICW-1),VARX*ZTP_CW(2,NICW-1)
                WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_CW(1,NICW),
     &            VARX*YTP_CW(1,NICW),VARX*ZTP_CW(1,NICW)
              ENDIF
            ENDIF
            IF( NICW.EQ.ID_CW(2,NCW) ) THEN
              WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_CW(2,NICW),
     &          VARX*YTP_CW(2,NICW),VARX*ZTP_CW(2,NICW)
            ENDIF
  122     CONTINUE
        ENDIF
  400 CONTINUE
!
!---  Loop over number of leaky wells  ---
!
      DO NLW = 1,N_LW
!
!---  Write a white well name text tag and geometry data to  
!     "well.dat" for an leaky well  ---
!
        NCH = INDEX( WNM_LW(NLW),'  ')-1
        NILW = ID_LW(1,NLW)
        WRITE(26,'(2A,3(1PE12.5,A),2A)') 
     &    'TEXT C=CYAN HU=POINT H=10 CS=GRID3D',
     &    ', X=',VARX*XTP_LW(1,NILW),
     &    ', Y=',VARX*YTP_LW(1,NILW),
     &    ', Z=',VARX*ZTP_LW(1,NILW),
     &    ', T="',WNM_LW(NLW)(1:NCH),'"'
        WRITE(26,'(A)') 'GEOMETRY T=LINE3D, C=CYAN LT=0.2'
        WRITE(26,'(A)') '1'
        NC = 0
        DO NILW = ID_LW(1,NLW),ID_LW(2,NLW)
          IF( NILW.GT.ID_LW(1,NLW) ) THEN
            IF( ABS(XTP_LW(1,NILW)-XTP_LW(2,NILW-1)).LT.1.D-3 .AND.
     &        ABS(YTP_LW(1,NILW)-YTP_LW(2,NILW-1)).LT.1.D-3 .AND.
     &        ABS(ZTP_LW(1,NILW)-ZTP_LW(2,NILW-1)).LT.1.D-3 ) THEN
              NC = NC + 1
            ELSE
              NC = NC + 2
            ENDIF
          ELSE
            NC = NC + 1
          ENDIF
          IF( NILW.EQ.ID_LW(2,NLW) ) THEN
            NC = NC + 1
          ENDIF
        ENDDO
        WRITE(FORM1(3:3),'(I1)') ICOUNT(NC)
        WRITE(26,FORM1) NC
        DO NILW = ID_LW(1,NLW),ID_LW(2,NLW)
          IF( NILW.EQ.ID_LW(1,NLW) ) THEN
            WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_LW(1,NILW),
     &        VARX*YTP_LW(1,NILW),VARX*ZTP_LW(1,NILW)
          ENDIF
          IF( NILW.GT.ID_LW(1,NLW) ) THEN
            IF( ABS(XTP_LW(1,NILW)-XTP_LW(2,NILW-1)).LT.1.D-3 .AND.
     &        ABS(YTP_LW(1,NILW)-YTP_LW(2,NILW-1)).LT.1.D-3 .AND.
     &        ABS(ZTP_LW(1,NILW)-ZTP_LW(2,NILW-1)).LT.1.D-3 ) THEN
              WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_LW(1,NILW),
     &          VARX*YTP_LW(1,NILW),VARX*ZTP_LW(1,NILW)
            ELSE
              WRITE(26,'(2(1PE12.5,1X),1PE12.5)') 
     &          VARX*XTP_LW(2,NILW-1),
     &          VARX*YTP_LW(2,NILW-1),VARX*ZTP_LW(2,NILW-1)
              WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_LW(1,NILW),
     &          VARX*YTP_LW(1,NILW),VARX*ZTP_LW(1,NILW)
            ENDIF
          ENDIF
          IF( NILW.EQ.ID_LW(2,NLW) ) THEN
            WRITE(26,'(2(1PE12.5,1X),1PE12.5)') VARX*XTP_LW(2,NILW),
     &        VARX*YTP_LW(2,NILW),VARX*ZTP_LW(2,NILW)
          ENDIF
        ENDDO
      ENDDO
!
!---  Close 'well.dat' file  ---
!
      CLOSE(UNIT=26)
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WR_WELL group  ---
!
      RETURN
      END



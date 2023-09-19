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
!     STOMP-EOR
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
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FILES
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
      REAL*8 XPX(5),YPX(5),ZPX(5)
      REAL*8 XIX(2),YIX(2),ZIX(2)
      REAL*8 AJ(3,3),BJ(3)
      INTEGER MSX(4,6),IJ(3)
      INTEGER N1X(4),N2X(4)
      CHARACTER*9 FORM1
      CHARACTER*8 FORM2
      CHARACTER*25 FORM5
      CHARACTER*40 FORM22
!
!----------------------Data Statements---------------------------------!
!
      DATA MSX / 1,2,4,3,1,5,6,2,1,3,7,5,2,6,8,4,3,4,8,7,5,7,8,6 /
      DATA N1X / 2,3,4,1 /
      DATA N2X / 1,2,3,4 /
      DATA FORM1 /'(2X,I1,$)'/
      DATA FORM2 /'(2X,I1)'/
      DATA FORM5 /'(A,I3,A,I3,A,I3,A,I8,A)'/
      DATA FORM22 /'(A,I3,A,I3,A,I3,A,I3,A,I3,A,I3,A,I8,A)'/
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
            IRX = 2**IBR(1,N)
            JRX = 2**IBR(2,N)
            KRX = 2**IBR(3,N)
            DO 470 KX = 1,KRX
            DO 460 JX = 1,JRX
            DO 450 IX = 1,IRX
              NX = IBR(4,N) + NDBR(IX,IRX,JX,JRX,KX)
              NC = 0
!
!---          Determine whether the transition points are within
!             the hexahedron node volume or on the hexahedron
!             surfaces  ---
!
              DO 190 NPT = 1,2
!
!---            Cylindrical coordinates with azimuthal symmetry,
!               centrally located wells  ---
!
                IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1
     &            .AND. I.EQ.1 ) THEN
                  XPX(1) = 0.D+0
                  XPX(2) = 0.D+0
                  YPX(1) = 0.D+0
                  YPX(2) = 0.D+0
                  ZPX(1) = ZE(1,NX)
                  ZPX(2) = ZE(5,NX)
!
!---              Node height greater than EPSLX  ---
!
                  IF( ABS(ZPX(1)-ZPX(2)).GT.EPSLX ) THEN
                    DZPX1 = ZTP_CW(NPT,NICW)-ZPX(1)
                    DZPX2 = ZPX(2)-ZTP_CW(NPT,NICW)
                    IF( ABS(DZPX1).LT.EPSLX ) DZPX1 = 0.D+0
                    IF( ABS(DZPX2).LT.EPSLX ) DZPX2 = 0.D+0
!
!---                Transition point within vertical limits of node  ---
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
!---            Check for point with hexahedron  ---
!
                CALL WITHIN( XTP_CW(NPT,NICW),YTP_CW(NPT,NICW),
     &            ZTP_CW(NPT,NICW),ICWX,NX )
!
!---            Opposing rotations found, point outside hexahedron  ---
!
                IF( ICWX.EQ.0 ) GOTO 190
!
!---            No opposing rotations found, point inside 
!               hexahedron  ---
!
                NC = NC+1
                XIX(NC) = XTP_CW(NPT,NICW)
                YIX(NC) = YTP_CW(NPT,NICW)
                ZIX(NC) = ZTP_CW(NPT,NICW)
  190         CONTINUE
!
!---          Both transition points inside hexahedron, skip
!             search for well path crossing hexahedron surfaces  ---
!
              IF( NC.EQ.2 ) GOTO 232
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1
     &          .AND. I.EQ.1 ) THEN
!
!---            Interval crosses lower node surface  ---
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
!---            Interval crosses upper node surface  ---
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
!---          Loop over node surfaces,
!             (bottom,south,west,east,north,top)  ---
!
              DO 230 NS = 1,6
!
!---            Define the five surface points, four corners
!               and one centroid---
!
                DO 200 NP = 1,4
                  MX = MSX(NP,NS)
!
!---              Cylindrical coordinates---
!
                  IF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
                    XPX(NP) = XE(MX,NX)*COS(YE(MX,NX))
                    YPX(NP) = XE(MX,NX)*SIN(YE(MX,NX))
                    ZPX(NP) = ZE(MX,NX)
                  ELSE
                    XPX(NP) = XE(MX,NX)
                    YPX(NP) = YE(MX,NX)
                    ZPX(NP) = ZE(MX,NX)
                  ENDIF
  200           CONTINUE
                NP = 4
                CALL PGCNTRD( NP,XPX(1),YPX(1),ZPX(1),
     &            XPX(5),YPX(5),ZPX(5) )
!
!
!---            Loop over the four triangular planes on the 
!               surface face  ---
!
                DO 220 NT = 1,4
!
!---              Load plane-line intersection matrix and 
!                 problem vector  ---
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
!---              Check for no intersection  ---
!
                  DO 210 IP = 1,3
                    AAMAX = 0.D+0
                    DO 202 JP = 1,3
                      IF( ABS(AJ(IP,JP)).GT.AAMAX ) AAMAX = 
     &                  ABS(AJ(IP,JP))
  202               CONTINUE
!
!---                No intersection go to next triangle on the 
!                   surface  ---
!
                    IF( ABS(AAMAX)/EPSL.LT.EPSL ) GOTO 220
  210             CONTINUE
!
!---              Find plane-line intersection matrix inverse  ---
!
                  JP = 3
                  KP = 3
                  CALL LUDCMP( AJ,JP,KP,IJ,DJ )
                  CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---              Find plane-line intersection point  ---
!
                  TX = BJ(1)
                  UX = BJ(2)
                  VX = BJ(3)
                  IF( ABS(TX).LT.EPSL ) TX = 0.D+0
                  IF( ABS(UX).LT.EPSL ) UX = 0.D+0
                  IF( ABS(VX).LT.EPSL ) VX = 0.D+0
!
!---              Line crosses surface, within the triangle  ---
!
                  IF( TX.GE.0.D+0 .AND. TX.LE.1.D+0 .AND.
     &              UX.GE.0.D+0 .AND. UX.LE.1.D+0 .AND. 
     &              VX.GE.0.D+0 .AND. VX.LE.1.D+0 .AND.
     &              (UX+VX).LE.1.D+0 ) THEN
                    XTX = XTP_CW(1,NICW) 
     &                + (XTP_CW(2,NICW)-XTP_CW(1,NICW))*TX
                    YTX = YTP_CW(1,NICW)
     &                + (YTP_CW(2,NICW)-YTP_CW(1,NICW))*TX
                    ZTX = ZTP_CW(1,NICW)
     &                + (ZTP_CW(2,NICW)-ZTP_CW(1,NICW))*TX
!
!---                Check for non-distinct points  ---
!
                    IF( NC.GE.1 ) THEN
                      DO 212 NDP = 1,NC
                        DPX = SQRT( ((XTX-XIX(NDP))**2) 
     &                    + ((YTX-YIX(NDP))**2) + ((ZTX-ZIX(NDP))**2) )
!
!---                    Duplicate point found  ---
!
                        IF( DPX.LT.EPSLX ) GOTO 214
  212                 CONTINUE
                    ENDIF
                    IF( NC.EQ.2 ) THEN
                      INDX = 7
                      CHMSG = 'Three Distinct Coupled Well Points'
     &                  //  ' at Node'
                      IMSG = NX
                      CALL WRMSGS( INDX )
                    ENDIF
                    NC = NC + 1
                    XIX(NC) = XTX
                    YIX(NC) = YTX
                    ZIX(NC) = ZTX
                  ENDIF
  214             CONTINUE
  220           CONTINUE
  230         CONTINUE
  232         CONTINUE
!
!---          Check that two well points are distinct  ---
!
              IF( NC.EQ.2 ) THEN
                DPX = SQRT( ((XIX(2)-XIX(1))**2) + ((YIX(2)-YIX(1))**2)
     &            + ((ZIX(2)-ZIX(1))**2) )
              ENDIF
!
!---          Two distinct well points within node define or update
!             well node trajectory points  ---
!
              IF( NC.EQ.2 .AND. DPX.GT.EPSLX ) THEN
!
!---            Cylindrical coordinates with azimuthal symmetry  ---
!
                IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) GOTO 362
!
!---            Check if line between well points is contained in
!               a node surface plane, loop over node surfaces,
!               (bottom,south,west,east,north,top)  ---
!
                DO 360 NS = 1,6
!
!---              Define the five surface points, four corners
!                 and one centroid---
!
                  DO 352 NP = 1,4
                    MX = MSX(NP,NS)
!
!---                Cylindrical coordinates---
!
                    IF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
                      XPX(NP) = XE(MX,NX)*COS(YE(MX,NX))
                      YPX(NP) = XE(MX,NX)*SIN(YE(MX,NX))
                      ZPX(NP) = ZE(MX,NX)
                    ELSE
                      XPX(NP) = XE(MX,NX)
                      YPX(NP) = YE(MX,NX)
                      ZPX(NP) = ZE(MX,NX)
                    ENDIF
  352             CONTINUE
                  NP = 4
                  CALL PGCNTRD( NP,XPX(1),YPX(1),ZPX(1),
     &              XPX(5),YPX(5),ZPX(5) )
!
!---              Loop over the four triangular planes on the 
!                 surface face  ---
!
                  DO 358 NT = 1,4
!
!---                Loop over trajectory points  ---
!
                    DO 356 NPT = 1,2
!
!---                  Load point-plane matrix  ---
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
!---                  Check for singular matrix  ---
!
                      DO 310 IP = 1,3
                        AAMAX = 0.D+0
                        DO 302 JP = 1,3
                          IF( ABS(AJ(IP,JP)).GT.AAMAX ) 
     &                    AAMAX = ABS(AJ(IP,JP))
  302                   CONTINUE
!
!---                    Singular matrix, trajectory point is within
!                       the surface plane  ---
!
                        IF( ABS(AAMAX)/EPSL.LT.EPSL ) GOTO 356
  310                 CONTINUE
!
!---                  Find matrix determinant  ---
!
                      JP = 3
                      KP = 3
                      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
                      DO 354 M = 1,JP
                        DJ = DJ*AJ(M,M)
  354                 CONTINUE
!
!---                  If determinant equals zero trajectory point
!                     is within the surface plane  ---
!
                      IF( ABS(DJ).GT.EPSL ) GOTO 358
  356               CONTINUE
                    INDX = 7
                    CHMSG = 'Both Well Trajectory Points within '
     &                //  'Surface Plane at Node'
                    IMSG = NX
                    CALL WRMSGS( INDX )
  358             CONTINUE
  360           CONTINUE
  362           CONTINUE
!
!---            Loop over well nodes in current well  ---
!
                DO 370 NWN = ID_CW(3,NCW),NWN_CW
!!
!!---            Loop over all well nodes  ---
!!
!                DO 370 NWN = 1,NWN_CW
!
!---              Well node previously counted  ---
!
                  IF( IWN_CW(NWN).EQ.NX ) THEN
                    GOTO 372
                  ENDIF
  370           CONTINUE
!
!---            Increment the number of field nodes with
!               coupled-well nodes  ---
!
                NWF_CW = NWF_CW + 1
                IF( NWF_CW.GT.LWF_CW ) THEN
                 INDX = 5
                 CHMSG = 'Number of Field Nodes with Coupled Wells ' //
     &             '> Parameter LWF_CW'
                 CALL WRMSGS( INDX )
                ENDIF
  372           CONTINUE
!
!---            Define a new well node and set the trajectory
!               points  ---
!
                NWN_CW = NWN_CW + 1
                IF( NWN_CW.GT.LWN_CW ) THEN
                 INDX = 5
                 CHMSG = 'Number of Coupled-Well Nodes ' //
     &             '> Parameter LWN_CW'
                 CALL WRMSGS( INDX )
                ENDIF
                XP_CW(1,NWN_CW) = XIX(1)
                YP_CW(1,NWN_CW) = YIX(1)
                ZP_CW(1,NWN_CW) = ZIX(1)
                XP_CW(2,NWN_CW) = XIX(2)
                YP_CW(2,NWN_CW) = YIX(2)
                ZP_CW(2,NWN_CW) = ZIX(2)
                IWN_CW(NWN_CW) = NX
                INV_CW(NWN_CW) = NICW
!
!---            Well projections on local grid  ---
!
                NWN = NWN_CW
                CALL PROJ_COUP_WELL( NX,NWN )
              ENDIF
  450       CONTINUE
  460       CONTINUE
  470       CONTINUE
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
          IF( NX.GT.NFLD ) THEN
            IRX = IBR(1,NX)
            JRX = IBR(2,NX)
            KRX = IBR(3,NX)
            IX = ID(INP(NX-NFLD))
            JX = JD(INP(NX-NFLD))
            KX = KD(INP(NX-NFLD))
            WRITE(FORM22(5:5),'(I1)') ICOUNT(IX)
            WRITE(FORM22(10:10),'(I1)') ICOUNT(JX)
            WRITE(FORM22(15:15),'(I1)') ICOUNT(KX)
            WRITE(FORM22(20:20),'(I1)') ICOUNT(IRX)
            WRITE(FORM22(25:25),'(I1)') ICOUNT(JRX)
            WRITE(FORM22(30:30),'(I1)') ICOUNT(KRX)
            WRITE(FORM22(35:35),'(I1)') ICOUNT(NX)
            WRITE(IWR,FORM22) ' (',IX,',',JX,',',KX,':',
     &        IRX,',',JRX,',',KRX,':',NX,')'
          ELSE
            IX = ID(NX)
            JX = JD(NX)
            KX = KD(NX)
            WRITE(FORM5(5:5),'(I1)') ICOUNT(IX)
            WRITE(FORM5(10:10),'(I1)') ICOUNT(JX)
            WRITE(FORM5(15:15),'(I1)') ICOUNT(KX)
            WRITE(FORM5(20:20),'(I1)') ICOUNT(NX)
            WRITE(IWR,FORM5) ' (',IX,',',JX,',',KX,':',NX,')'
          ENDIF
  710   CONTINUE
  720 CONTINUE
!
!---  Record coupled-well nodes to output file ---
!
      WRITE(IWR,'(//,A,/)') ' --- Coupled-Well Transmissibility ' // 
     &  'Record  ---'
      DO 740 NCW = 1,N_CW
        NCH = INDEX( WNM_CW(NCW),'  ') - 1
        WRITE(IWR,'(2X,A,I4,2A)') 'Nodes Containing Coupled-Well ' //
     &    'Number: ',NCW,' Coupled-Well Name: ',WNM_CW(NCW)(1:NCH)
        WRITE(IWR,'(A,$)') '    I     '
        WRITE(IWR,'(A,$)') '    J     '
        WRITE(IWR,'(A,$)') '    K     '
        WRITE(IWR,'(A,$)') '   Xs (m)    '
        WRITE(IWR,'(A,$)') '   Xe (m)    '
        WRITE(IWR,'(A,$)') '   Ys (m)    '
        WRITE(IWR,'(A,$)') '   Ye (m)    '
        WRITE(IWR,'(A,$)') '   Zs (m)    '
        WRITE(IWR,'(A,$)') '   Ze (m)    '
        WRITE(IWR,'(A,$)') '  Tx (m^3)   '
        WRITE(IWR,'(A,$)') '  Ty (m^3)   '
        WRITE(IWR,'(A,$)') '  Tz (m^3)   '
        WRITE(IWR,'(A,$)') ' Teff (m^3)  '
        WRITE(IWR,'(A,$)') '  khx (m^3)  '
        WRITE(IWR,'(A,$)') '  khy (m^3)  '
        WRITE(IWR,'(A,$)') '  khz (m^3)  '
        WRITE(IWR,'(A,$)') '  rox (m)    '
        WRITE(IWR,'(A,$)') '  roy (m)    '
        WRITE(IWR,'(A,$)') '  roz (m)    '
        WRITE(IWR,'(A,$)') '   rw (m)    '
        WRITE(IWR,'(A)')   '     s       '
!
!---    Loop over the nodes in the coupled well ---
!
        DO 730 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(NWN)
          INVX = INV_CW(NWN)
          IZN = IZ(N)
          I = ID(N)
          J = JD(N)
          K = KD(N)
          XLX = PLX_CW(NWN)
          YLX = PLY_CW(NWN)
          ZLX = PLZ_CW(NWN)
          PERMX = MAX( PERM(1,IZN),1.D-20 )
          PERMY = MAX( PERM(2,IZN),1.D-20 )
          PERMZ = MAX( PERM(3,IZN),1.D-20 )
          RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---      Cylindrical coordinates with azimuthal symmetry,
!         centrally located wells  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            ROZ = RP(I)          
            
!
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROZ = PAR_CW(3,INVX)
            ENDIF
            ROX = ROZ
            ROY = ROZ
            RWX = MIN( RWX,9.999D-1*ROZ )
            PERMX = PERMRF(2,N)*PERM(1,IZN)
            WI_CWX = 2.D+0*GPI*PERMX*ZLX/(LOG(ROZ/RWX)+PAR_CW(1,INVX))
            WIX = WI_CWX
            WIY = WI_CWX
            WIZ = WI_CWX
            EKHX = PERMX*ZLX
            EKHY = PERMX*ZLX
            EKHZ = PERMX*ZLX
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
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROX = PAR_CW(3,INVX)
              ROY = PAR_CW(3,INVX)
              ROZ = PAR_CW(3,INVX)
            ENDIF
!
!---        Well index components  ---
!
            PERMX = PERMRF(2,N)*PERM(1,IZN)
            PERMY = PERMRF(2,N)*PERM(2,IZN)
            PERMZ = PERMRF(2,N)*PERM(3,IZN)
            EKHX = SQRT(PERMY*PERMZ)*XLX
            EKHY = SQRT(PERMX*PERMZ)*YLX
            EKHZ = SQRT(PERMX*PERMY)*ZLX
            WIX = 2.D+0*GPI*EKHX/(LOG(ROX/RWX)+PAR_CW(1,INVX))
            WIY = 2.D+0*GPI*EKHY/(LOG(ROY/RWX)+PAR_CW(1,INVX))
            WIZ = 2.D+0*GPI*EKHZ/(LOG(ROZ/RWX)+PAR_CW(1,INVX))
            WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
          ENDIF
          IF( N.GT.NFLD ) THEN
            IX = ID(INP(N-NFLD))
            JX = JD(INP(N-NFLD))
            KX = KD(INP(N-NFLD))
          ELSE
            IX = ID(N)
            JX = JD(N)
            KX = KD(N)
          ENDIF
          WRITE(IWR,'(1X,I8,1X,$)') IX
          WRITE(IWR,'(1X,I8,1X,$)') JX
          WRITE(IWR,'(1X,I8,1X,$)') KX
          WRITE(IWR,'(1X,1PE12.5,1X,$)') XP_CW(1,NWN)
          WRITE(IWR,'(1X,1PE12.5,1X,$)') XP_CW(2,NWN)
          WRITE(IWR,'(1X,1PE12.5,1X,$)') YP_CW(1,NWN)
          WRITE(IWR,'(1X,1PE12.5,1X,$)') YP_CW(2,NWN)
          WRITE(IWR,'(1X,1PE12.5,1X,$)') ZP_CW(1,NWN)
          WRITE(IWR,'(1X,1PE12.5,1X,$)') ZP_CW(2,NWN)
          WRITE(IWR,'(1X,1PE12.5,1X,$)') WIX
          WRITE(IWR,'(1X,1PE12.5,1X,$)') WIY
          WRITE(IWR,'(1X,1PE12.5,1X,$)') WIZ
          WRITE(IWR,'(1X,1PE12.5,1X,$)') WI_CWX
          WRITE(IWR,'(1X,1PE12.5,1X,$)') EKHX
          WRITE(IWR,'(1X,1PE12.5,1X,$)') EKHY
          WRITE(IWR,'(1X,1PE12.5,1X,$)') EKHZ
          WRITE(IWR,'(1X,1PE12.5,1X,$)') ROX
          WRITE(IWR,'(1X,1PE12.5,1X,$)') ROY
          WRITE(IWR,'(1X,1PE12.5,1X,$)') ROZ
          WRITE(IWR,'(1X,1PE12.5,1X,$)') RWX
          WRITE(IWR,'(1X,1PE12.5,1X)') PAR_CW(1,INVX)
  730   CONTINUE
  740 CONTINUE
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
!     STOMP-EOR
!
!     Equilibrate coupled-well pressure with reservoir.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 May 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
!---  Saturation weighted fluid density  ---
!
      IF( SG(2,N).GT.0.D+0 ) THEN
        RHOX = RHOG(2,N)
      ELSE
        RHOX = SN(2,N)*RHON(2,N) + SL(2,N)*RHOL(2,N)
      ENDIF
      P_CW(2,NCW) = PSO(2,N) - (ZCWX-ZP(N))*GRAV*RHOX
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
!     STOMP-EOR
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
        IF( IT_CW(NCW).GT.0 ) THEN
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
      SUBROUTINE INJ_BO_WELL( NCW )
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
!     Injection black-oil well model
!     
!     Rate controlled or pressure controlled
!
!     Flux of water mass, gas mass, and oil mass from coupled-well 
!     nodes to field nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June 2014.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE JACOB
      USE HYST
      USE GRID
      USE GLB_PAR
      USE FDVP
      USE COUP_WELL
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XPX(2),YPX(2),ZPX(2)
      REAL*8 VAR_CWX(6+LNGC)
      REAL*8 XGC_CWX(LNGC),XNC_CWX(LNGC)
      INTEGER, SAVE :: IALLOC
      DATA IALLOC /0/
      INTEGER, SAVE, DIMENSION(:), ALLOCATABLE :: MCW,MFD
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INJ_BO_WELL'
      IF( IALLOC.EQ.0 ) THEN
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
   10   CONTINUE        
        IALLOC = 1
      ENDIF
!
!---  Zero injection well fluxes ---
!
      QM_CW(1,NCW) = 0.D+0
      QM_CW(3,NCW) = 0.D+0
      QM_CW(5,NCW) = 0.D+0
      DO 11 IGC = 1,NGC+2
        QM_CW((7+(IGC-1)*2),NCW) = 0.D+0
   11 CONTINUE
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
!       Q_CW(4,NWN) - nonaqueous-liquid volumetric flux, m^3/s  ---
!
        DO M = 1,4
          Q_CW(M,NWN) = 0.D+0
        ENDDO
!
!---    Loop over increment indices  ---
!
        DO 20 M = 1,ISVC+2
          DO 12 IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = 0.D+0
   12     CONTINUE
          FXW_CW(M,NWN) = 0.D+0
   20   CONTINUE
   30  CONTINUE
!
!---  Injection well time interval ---
!
      DQ_CWX = 1.D-6
      TMZ = TM
      IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---  Cyclic time periods  ---
!
      IF( ICC_CW(NCW).EQ.1 ) THEN
!
!---    Loop over the injection well time periods, to find the 
!       final well time  ---
!
        NTX = 0
        DO 40 NTP = 1,IM_CW(NCW)
          NTX = NTX + IMP_CW(NTP,NCW)
   40   CONTINUE
!
!---    Determine time with the cyclic time period  ---
!
        TMZ = MOD( TM,VAR_CW(1,NTX,NCW) )
        IF( TM.GT.VAR_CW(1,NTX,NCW) ) THEN
          IF( TMZ.LT.EPSL ) TMZ = VAR_CW(1,NTX,NCW)
        ENDIF
      ENDIF
!
!---  Injection well is inactive set well pressure to be in 
!     equilibrium with reservoir  ---
!
      IF( TMZ.LE.VAR_CW(1,1,NCW) ) THEN
        CALL EQUIL_COUP_WELL( NCW )
        ID_CW(8,NCW) = 1
        GOTO 500
      ENDIF
!
!---  Loop over the injection well time periods  ---
!
      NS = 1
      DO 120 NTP = 1,IM_CW(NCW)
!
!---    Injection well time period only has one time (start time)  ---
!
        IF( IMP_CW(NTP,NCW).EQ.1 ) THEN
!
!---      Time prior to start time, injection well is inactive, set 
!         well pressure to be in equilibrium with reservoir  ---
!
          IF( TMZ.LE.VAR_CW(1,NS,NCW) ) THEN
            CALL EQUIL_COUP_WELL( NCW )
            ID_CW(8,NCW) = 1
            GOTO 500
!
!---      Time after start time, injection well is active  ---
!
          ELSE
            DO 80 N = 2,7
              VAR_CWX(N) = VAR_CW(N,1,NCW)
   80       CONTINUE
            GOTO 122
          ENDIF
!
!---    Injection well time period only has multiple times  ---
!
        ELSE
          NE = NS + IMP_CW(NTP,NCW) - 1
!
!---      Time outside of injection well time period, go to next 
!         injection well time period  ---
!
          IF( TMZ.LE.VAR_CW(1,NS,NCW) .OR. 
     &      TMZ.GT.VAR_CW(1,NE,NCW) ) GOTO 110 
!
!---      Injection well time period has multiple time points, use  
!         linear interpolation of well parameters between time points  ---
!
          DO 100 M = 2,IMP_CW(NTP,NCW)
            MX = NS + M - 1
            IF( TMZ.LE.VAR_CW(1,MX,NCW) ) THEN
              TD_CW = VAR_CW(1,MX,NCW)-VAR_CW(1,MX-1,NCW)
              DT_CW = MIN( VAR_CW(1,MX,NCW)-TMZ,DT )
              TF_CW = (TMZ-VAR_CW(1,MX-1,NCW))/TD_CW
              DO 90 N = 2,7
                VAR_CWX(N) = VAR_CW(N,MX-1,NCW) + 
     &            TF_CW*(VAR_CW(N,MX,NCW)-VAR_CW(N,MX-1,NCW))
   90         CONTINUE
              GOTO 122
            ENDIF
  100     CONTINUE
        ENDIF
  110   CONTINUE
        NS = NS + IMP_CW(NTP,NCW)
  120 CONTINUE
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
        NPHAZX = 0
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          T_CWX = VAR_CWX(4)
        ELSE
          T_CWX = T(2,N)
        ENDIF
!
!---    Well state #1  ---
!
!       SL = 1.0
!       SG+SN = 0.0
!
!       Declared variables:
!
!       VAR_CWX(6) - aqueous salt relative saturation, or
!                    aqueous salt mass fraction
!
        IF( ITS_CW(NTP,NCW)/100.EQ.1 ) THEN
          PVAX = VAR_CWX(5)
          TMSX = VAR_CWX(6)
          ISR5X = MOD(ITS_CW(NTP,NCW)/10,10)
          CALL BO_IC1( P_CWX,RHOL_CWX,TMSX,T_CWX,
     &      XLS_CWX,XLW_CWX,YLSX,ISR5X )
          SG_CWX = 0.D+0
          SN_CWX = 0.D+0
          SL_CWX = 1.D+0
          RHOF_CWX = RHOL_CWX
!
!---    Well state #2  ---
!
!       SL < 1.0
!       SG+SN > 0.0
!
!       Declared variables:
!
!       VAR_CWX(5) - aqueous saturation
!       VAR_CWX(6) - aqueous salt relative saturation, or
!                    aqueous salt mass fraction
!       VAR_CWX(7) - gas mass fraction of oil + gas
!
        ELSEIF( ITS_CW(NTP,NCW)/100.EQ.2 ) THEN
          SL_CWX = VAR_CWX(5)
          TMSX = VAR_CWX(6)
          XGX = VAR_CWX(7)
          RHOGSX = 1.225D+0*GCPP(37,1)
          RHONSX = 1.D+3*GCPP(37,2)
          PX = P_CWX + PATM
          IF( 1.D+0-XGX.GT.EPSL ) THEN
            POSM_CWX = XGX*RHONSX/((1.D+0-XGX)*RHOGSX)
            CALL BO_BP( PBPX,POSM_CWX,T_CWX )
            IF( PX.GT.PBPX ) THEN
              SN_CWX = 1.D+0
              SG_CWX = 0.D+0
            ELSE
              PVWX = 0.D+0
              CALL BO_EQ_GAS( PX,PVWX,RHOGX,T_CWX,XGC_CWX,XGWX )
              CALL BO_EQ_OIL( RHONX,PX,POSM_CWX,T_CWX,XNC_CWX )
              SN_CWX = RHOGX*(XGX*(1.D+0-SL_CWX)-XGC_CWX(1)*
     &          (1.D+0-SL_CWX))/(RHONX*(XNC_CWX(1)-XGX) - 
     &          RHOGX*(XGC_CWX(1)-XGX))
              SG_CWX = MAX( 1.D+0-SL_CWX-SN_CWX,0.D+0 )
            ENDIF
          ELSE
            SG_CWX = 1.D+0
            SN_CWX = 0.D+0
          ENDIF
          ISR5X = MOD(ISRT(NS)/10,10)
          IZNX = 0
          CALL BO_IC2( POSM_CWX,P_CWX,RHOG_CWX,RHOL_CWX,RHON_CWX,SL_CWX,
     &      TMSX,T_CWX,XGC_CWX,XGW_CWX,XLS_CWX,XLW_CWX,XNC_CWX,
     &      YLSX,ISR5X,IZNX )
          RHOF_CWX = SG_CWX*RHOG_CWX + SL_CWX*RHOL_CWX + 
     &      SN_CWX*RHON_CWX
!
!---    Well state #3  ---
!
!       SL = 0.0
!       SG+SN = 1.0
!
!       Declared variables:
!
!       VAR_CWX(5) - water vapor relative saturation
!       VAR_CWX(7) - gas mass fraction of oil + gas
!      
        ELSEIF( ITS_CW(NTP,NCW)/100.EQ.3 ) THEN
          SL_CWX = 0.D+0
          PVWX = VAR_CWX(5)
          YLSX = 0.D+0
          XGX = VAR_CWX(7)
          RHOGSX = 1.225D+0*GCPP(37,1)
          RHONSX = 1.D+3*GCPP(37,2)
          PX = P_CWX + PATM
          IF( 1.D+0-XGX.GT.EPSL ) THEN
            POSM_CWX = XGX*RHONSX/((1.D+0-XGX)*RHOGSX)
            CALL BO_BP( PBPX,POSM_CWX,T_CWX )
            IF( PX.GT.PBPX ) THEN
              SN_CWX = 1.D+0
              SG_CWX = 0.D+0
            ELSE
              PVWX = 0.D+0
              CALL BO_EQ_GAS( PX,PVWX,RHOGX,T_CWX,XGC_CWX,XGWX )
              CALL BO_EQ_OIL( RHONX,PX,POSM_CWX,T_CWX,XNC_CWX )
              SN_CWX = RHOGX*(XGX*(1.D+0-SL_CWX)-XGC_CWX(1)*
     &          (1.D+0-SL_CWX))/(RHONX*(XNC_CWX(1)-XGX) - 
     &          RHOGX*(XGC_CWX(1)-XGX))
              SG_CWX = MAX( 1.D+0-SL_CWX-SN_CWX,0.D+0 )
            ENDIF
          ELSE
            SG_CWX = 1.D+0
            SN_CWX = 0.D+0
          ENDIF
          IZNX = 0
          CALL BO_IC3( POSM_CWX,PVWX,P_CWX,RHOG_CWX,RHON_CWX,
     &      T_CWX,XGC_CWX,XGW_CWX,XNC_CWX,YLSX,IZNX )
          RHOF_CWX = SG_CWX*RHOG_CWX + SN_CWX*RHON_CWX
        ENDIF
!
!---    Load local variable for mass flux, kg/s  ---
!
        IF( M.EQ.1 ) THEN
!
!---      Convert volumetric flow rate, m^3/s to mass flow rate, kg/s,
!         using fluid density at top of screened interval  ---
!
          IF( IT_CW(NCW).EQ.2 ) THEN
            VAR_CW2X = VAR_CWX(2)*RHOF_CWX
          ELSE
            VAR_CW2X = VAR_CWX(2)
          ENDIF
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
          NPHAZX = 0
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            T_CWX = VAR_CWX(4)
          ELSE
            T_CWX = T(2,N)
          ENDIF
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
          P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOF_CWX
!
!---      Well pressure at the node centroid, used for coupled-well
!         nodal output  ---
!
          IF( M.EQ.1 ) THEN
            NWF = IWP_CW(NWN)
            PF_CW(NWF) = P_CWX - (ZPNX-ZPX(1))*GRAV*RHOF_CWX
          ENDIF
!
!---      Adjust the reservoir pressure to the coupled-well node
!         centroid  ---
!
!
!---      Phase condition #1  ---
!
!         SL = 1.0, SG+SN = 0.0
!
          IF( MOD(NPHAZ(2,N),10).EQ.1 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #2  ---
!
!         SL < 1.0, ST > 0.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.2 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #3  ---
!
!         SL < 1.0, ST > 0.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.3 ) THEN
            PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(MF,N)
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
!
!---      Phase condition #4  ---
!
!         SL = 0.0, ST = 1.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.4 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #5  ---
!
!         SL = 0.0, ST = 1.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.5 ) THEN
            PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(MF,N)
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
          ENDIF
!
!---      Well state #1  ---
!
!         SL = 1.0
!         SG+SN = 0.0
!
!         Declared variables:
!
!         VAR_CWX(6) - aqueous salt relative saturation, or
!                           aqueous salt mass fraction
!
          IF( ITS_CW(NTP,NCW)/100.EQ.1 ) THEN
            PVAX = VAR_CWX(5)
            TMSX = VAR_CWX(6)
            ISR5X = MOD(ITS_CW(NTP,NCW)/10,10)
            CALL BO_IC1( P_CWX,RHOL_CWX,TMSX,T_CWX,
     &        XLS_CWX,XLW_CWX,YLSX,ISR5X )
            SG_CWX = 0.D+0
            SN_CWX = 0.D+0
            SL_CWX = 1.D+0
            RHOF_CWX = RHOL_CWX
!
!---      Well state #2  ---
!
!         SL < 1.0
!         SG+SN > 0.0
!
!         Declared variables:
!
!         VAR_CWX(5) - aqueous saturation
!         VAR_CWX(6) - aqueous salt relative saturation, or
!                      aqueous salt mass fraction
!         VAR_CWX(7) - gas mass fraction of oil + gas
!
          ELSEIF( ITS_CW(NTP,NCW)/100.EQ.2 ) THEN
            SL_CWX = VAR_CWX(5)
            TMSX = VAR_CWX(6)
            XGX = VAR_CWX(7)
            RHOGSX = 1.225D+0*GCPP(37,1)
            RHONSX = 1.D+3*GCPP(37,2)
            PX = P_CWX + PATM
            IF( 1.D+0-XGX.GT.EPSL ) THEN
              POSM_CWX = XGX*RHONSX/((1.D+0-XGX)*RHOGSX)
              CALL BO_BP( PBPX,POSM_CWX,T_CWX )
              IF( PX.GT.PBPX ) THEN
                SN_CWX = 1.D+0
                SG_CWX = 0.D+0
              ELSE
                PVWX = 0.D+0
                CALL BO_EQ_GAS( PX,PVWX,RHOGX,T_CWX,XGC_CWX,XGWX )
                CALL BO_EQ_OIL( RHONX,PX,POSM_CWX,T_CWX,XNC_CWX )
                SN_CWX = RHOGX*(XGX*(1.D+0-SL_CWX)-XGC_CWX(1)*
     &            (1.D+0-SL_CWX))/(RHONX*(XNC_CWX(1)-XGX) - 
     &            RHOGX*(XGC_CWX(1)-XGX))
                SG_CWX = MAX( 1.D+0-SL_CWX-SN_CWX,0.D+0 )
              ENDIF
            ELSE
              SG_CWX = 1.D+0
              SN_CWX = 0.D+0
            ENDIF
            ISR5X = MOD(ISRT(NS)/10,10)
            IZNX = 0
            CALL BO_IC2( POSM_CWX,P_CWX,RHOG_CWX,RHOL_CWX,RHON_CWX,
     &        SL_CWX,TMSX,T_CWX,XGC_CWX,XGW_CWX,XLS_CWX,XLW_CWX,XNC_CWX,
     &        YLSX,ISR5X,IZNX )
            RHOF_CWX = SG_CWX*RHOG_CWX + SL_CWX*RHOL_CWX + 
     &        SN_CWX*RHON_CWX
!
!---      Well state #3  ---
!
!         SL = 0.0
!         SG+SN = 1.0
!
!         Declared variables:
!
!         VAR_CWX(5) - water vapor relative saturation
!         VAR_CWX(7) - gas mass fraction of oil + gas
!        
          ELSEIF( ITS_CW(NTP,NCW)/100.EQ.3 ) THEN
            SL_CWX = 0.D+0
            PVWX = VAR_CWX(5)
            YLSX = 0.D+0
            XLS_CWX = 0.D+0
            XGX = VAR_CWX(7)
            RHOGSX = 1.225D+0*GCPP(37,1)
            RHONSX = 1.D+3*GCPP(37,2)
            PX = P_CWX + PATM
            IF( 1.D+0-XGX.GT.EPSL ) THEN
              POSM_CWX = XGX*RHONSX/((1.D+0-XGX)*RHOGSX)
              CALL BO_BP( PBPX,POSM_CWX,T_CWX )
              IF( PX.GT.PBPX ) THEN
                SN_CWX = 1.D+0
                SG_CWX = 0.D+0
              ELSE
                PVWX = 0.D+0
                CALL BO_EQ_GAS( PX,PVWX,RHOGX,T_CWX,XGC_CWX,XGWX )
                CALL BO_EQ_OIL( RHONX,PX,POSM_CWX,T_CWX,XNC_CWX )
                SN_CWX = RHOGX*(XGX*(1.D+0-SL_CWX)-XGC_CWX(1)*
     &            (1.D+0-SL_CWX))/(RHONX*(XNC_CWX(1)-XGX) - 
     &            RHOGX*(XGC_CWX(1)-XGX))
                SG_CWX = MAX( 1.D+0-SL_CWX-SN_CWX,0.D+0 )
              ENDIF
            ELSE
              SG_CWX = 1.D+0
              SN_CWX = 0.D+0
            ENDIF
            IZNX = 0
            CALL BO_IC3( POSM_CWX,PVWX,P_CWX,RHOG_CWX,RHON_CWX,
     &        T_CWX,XGC_CWX,XGW_CWX,XNC_CWX,YLSX,IZNX )
              RHOF_CWX = SG_CWX*RHOG_CWX + SN_CWX*RHON_CWX
          ENDIF
!
!---      Well gas viscosity, Pa s  ---
!
          IF( SG_CWX.GT.0.D+0 ) THEN
            CALL BO_VISG( RHOG_CWX,VISG_CWX )
          ENDIF
!
!---      Well aqueous viscosity, Pa s
!         concentration of CO2 in mol/kg,
!         concentration of salt in mol/kg  ---
!
          IF( SL_CWX.GT.0.D+0 ) THEN
            CMLAX = 0.D+0
            CMLSX = 1.D+3*XLS_CWX*RHOL_CWX/WTMS
            CALL VISC_L( CMLAX,CMLSX,P_CWX,T_CWX,VISL_CWX )
          ENDIF
!
!---      Well nonaqueous-liquid viscosity, Pa s  ---
!
          IF( SN_CWX.GT.0.D+0 ) THEN
            CALL BO_VISN( PX,POSM_CWX,T_CWX,VISN_CWX )
          ENDIF
!
!---      Equivalent field node radius components  ---
!
          PERMX = MAX( PERM(1,IZN),1.D-20 )
          PERMY = MAX( PERM(2,IZN),1.D-20 )
          PERMZ = MAX( PERM(3,IZN),1.D-20 )
          RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---      Cylindrical coordinates with azimuthal symmetry,
!         centrally located wells  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            ROZ = RP(I)
!
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROZ = PAR_CW(3,INVX)
            ENDIF
            RWX = MIN( RWX,9.999D-1*ROZ )
            PERMX = PERMRF(MF,N)*PERM(1,IZN)
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
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROX = PAR_CW(3,INVX)
              ROY = PAR_CW(3,INVX)
              ROZ = PAR_CW(3,INVX)
            ENDIF
!
!---        Well index components  ---
!
            PERMX = PERMRF(MF,N)*PERM(1,IZN)
            PERMY = PERMRF(MF,N)*PERM(2,IZN)
            PERMZ = PERMRF(MF,N)*PERM(3,IZN)
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
          DPNX = MAX( P_CWX-PGFX,0.D+0 )
          DPLX = MAX( P_CWX-PGFX,0.D+0 )
!
!---      Zero fluxes from well to reservoir  ---
!
          FXW_CW(M,NWN) = 0.D+0
          DO 150 IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = 0.D+0
  150     CONTINUE
!
!---      Gas flux from well to reservoir  ---
!
          IF( SG_CWX.GT.0.D+0 ) THEN
            FGX_CWX = SG_CWX*WI_CWX*RHOG_CWX*DPGX/VISG_CWX
!
!---        Water component gas flux from well to reservoir  ---
!
            FXW_CW(M,NWN) = FXW_CW(M,NWN) + FGX_CWX*XGW_CWX
!
!---        Gas component gas flux from well to reservoir  ---
!
            FXC_CW(1,M,NWN) = FXC_CW(1,M,NWN) +  FGX_CWX*XGC_CWX(1)
          ENDIF
!
!---      Aqueous flux from well to reservoir  ---
!
          IF( SL_CWX.GT.0.D+0 ) THEN
            FLX_CWX = SL_CWX*WI_CWX*RHOL_CWX*DPLX/VISL_CWX
!
!---        Water component aqueous flux from well to reservoir  ---
!
            FXW_CW(M,NWN) = FXW_CW(M,NWN) + FLX_CWX*XLW_CWX
          ENDIF
!
!---      Nonaqueous-liquid flux from well to reservoir  ---
!
          IF( SN_CWX.GT.0.D+0 ) THEN
            FNX_CWX = SN_CWX*WI_CWX*RHON_CWX*DPNX/VISN_CWX
!
!---        Gas component Nonaqueous-liquid flux from well to 
!           reservoir  ---
!
            FXC_CW(1,M,NWN) = FXC_CW(1,M,NWN) +  FNX_CWX*XNC_CWX(1)
!
!---        Oil component Nonaqueous-liquid flux from well to 
!           reservoir  ---
!
            FXC_CW(2,M,NWN) = FXC_CW(2,M,NWN) +  FNX_CWX*XNC_CWX(2)
          ENDIF
!
!---      Store current coupled-well node location in previous
!         coupled-well node location  ---
!
          XPX(1) = XPX(2)
          YPX(1) = YPX(2)
          ZPX(1) = ZPX(2)
!
!---      Volumetric injection well fluxes  ---
!
          IF( M.EQ.1 ) THEN
            IF( SL_CWX.GT.0.D+0 ) Q_CW(2,NWN) = FLX_CWX/RHOL_CWX
            IF( SG_CWX.GT.0.D+0 ) Q_CW(3,NWN) = FGX_CWX/RHOG_CWX
            IF( SN_CWX.GT.0.D+0 ) Q_CW(4,NWN) = FNX_CWX/RHON_CWX
            Q_CW(1,NWN) = Q_CW(2,NWN) + Q_CW(3,NWN) + Q_CW(4,NWN)
          ENDIF
  200   CONTINUE
  300 CONTINUE
!
!---  Mass balance residuals for injection type coupled well  ---
!
      NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
      MX = (NWFX*ISVC)+2
      RS_CW(1,NCW) = VAR_CW2X
      RS_CW(MX,NCW) = VAR_CW2X
      QM_CW(1,NCW) = 0.D+0
      QM_CW(3,NCW) = 0.D+0
      QM_CW(5,NCW) = 0.D+0
      DO 310 IGC = 1,NGC+2
        QM_CW((7+(IGC-1)*2),NCW) = 0.D+0
  310 CONTINUE
      QM_CWX = 0.D+0
      QTM_CWX = 0.D+0
!
!---  Loop over coupled-well nodes  ---
!
      DO 330 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
        RS_CW(1,NCW) = RS_CW(1,NCW) - FXW_CW(1,NWN)
        RS_CW(MX,NCW) = RS_CW(MX,NCW) - FXW_CW(ISVC+2,NWN)
        QTM_CWX = QTM_CWX + FXW_CW(1,NWN)
        QM_CW(3,NCW) = QM_CW(3,NCW) + FXW_CW(1,NWN)
        QM_CWX = QM_CWX + FXW_CW(ISVC+2,NWN)
        DO 320 IGC = 1,NGC+2
          RS_CW(1,NCW) = RS_CW(1,NCW) - FXC_CW(IGC,1,NWN)
          RS_CW(MX,NCW) = RS_CW(MX,NCW) - FXC_CW(IGC,ISVC+2,NWN)
          QTM_CWX = QTM_CWX + FXC_CW(IGC,1,NWN)
          QM_CW(1,NCW) = QM_CW(1,NCW) + FXC_CW(IGC,1,NWN)
          QM_CW((7+(IGC-1)*2),NCW) = QM_CW((7+(IGC-1)*2),NCW) + 
     &      FXC_CW(IGC,1,NWN)
          QM_CWX = QM_CWX + FXC_CW(IGC,ISVC+2,NWN)
  320   CONTINUE
  330 CONTINUE
!
!---  Loop over field nodes that contain coupled-well nodes  ---
!
      DO 470 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
        M1 = (NWF-ID_CW(5,NCW))*ISVC + 1
        DO 460 M2 = 1,ISVC
          M3 = M1+M2
          RS_CW(M3,NCW) = VAR_CW2X
!
!---      Loop over coupled-well nodes  ---
!
          DO 450 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---        If coupled-well node is within the current field
!           node, use incremented fluxes  ---
!
            IF( IWF_CW(NWF).EQ.IWN_CW(NWN) ) THEN
              RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXW_CW(M2+1,NWN) 
              DO 440 IGC = 1,NGC+2
                RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXC_CW(IGC,M2+1,NWN)
  440         CONTINUE
!
!---        If coupled-well node is outside the current field
!           node, use un-incremented fluxes  ---
!
            ELSE
              RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXW_CW(1,NWN)
              DO 442 IGC = 1,NGC+2
                RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXC_CW(IGC,1,NWN)
  442         CONTINUE
            ENDIF
  450     CONTINUE
  460   CONTINUE
  470 CONTINUE
  500 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INJ_BO_WELL group  ---
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
!     STOMP-EOR
!
!     Injection coupled well model
!     
!     Rate controlled or pressure controlled
!
!     Flux of energy, water mass, CO2 mass, CH4 mass, and petroleum
!     component masses from coupled-well nodes to field nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 April 2014.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE JACOB
      USE HYST
      USE GRID
      USE GLB_PAR
      USE FDVP
      USE FDVGC
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
      REAL*8 VAR_CWX(6+LNGC)
      REAL*8 ZMCX(LNGC)
      REAL*8 XGC_CWX(LNGC),XNC_CWX(LNGC),XMGCX(LNGC),XMNCX(LNGC)
      INTEGER, SAVE :: IALLOC
      DATA IALLOC /0/
      INTEGER, SAVE, DIMENSION(:), ALLOCATABLE :: MCW,MFD
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INJ_COUP_WELL'
!
!---  Black-oil option  ---
!
      IF( ISLC(68).EQ.0 ) THEN
        CALL INJ_BO_WELL( NCW )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compositional option  ---
!
      IF( IALLOC.EQ.0 ) THEN
!
!---    Dynamic memory allocation  ---
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
   10   CONTINUE        
        IALLOC = 1
      ENDIF
!
!---  K-factors computed from fugacity coefficients  ---
!
      IFK = 0
!
!---  Flow controlled well ---
!
      ID_CW(8,NCW) = 0
!
!---  Zero injection well fluxes ---
!
      QM_CW(1,NCW) = 0.D+0
      QM_CW(3,NCW) = 0.D+0
      QM_CW(5,NCW) = 0.D+0
      DO 11 IGC = 1,NGC+2
        QM_CW((7+(IGC-1)*2),NCW) = 0.D+0
   11 CONTINUE
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
!       Q_CW(4,NWN) - nonaqueous-liquid volumetric flux, m^3/s  ---
!
        DO M = 1,4
          Q_CW(M,NWN) = 0.D+0
        ENDDO
!
!---    Loop over increment indices  ---
!
        DO 20 M = 1,ISVC+2
          DO 12 IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = 0.D+0
   12     CONTINUE
          FXE_CW(M,NWN) = 0.D+0
          FXW_CW(M,NWN) = 0.D+0
   20   CONTINUE
   30  CONTINUE
!
!---  Injection well time interval ---
!
      DQ_CWX = 1.D-6
      TMZ = TM
      IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---  Cyclic time periods  ---
!
      IF( ICC_CW(NCW).EQ.1 ) THEN
!
!---    Loop over the injection well time periods, to find the 
!       final well time  ---
!
        NTX = 0
        DO 40 NTP = 1,IM_CW(NCW)
          NTX = NTX + IMP_CW(NTP,NCW)
   40   CONTINUE
!
!---    Determine time with the cyclic time period  ---
!
        TMZ = MOD( TM,VAR_CW(1,NTX,NCW) )
        IF( TM.GT.VAR_CW(1,NTX,NCW) ) THEN
          IF( TMZ.LT.EPSL ) TMZ = VAR_CW(1,NTX,NCW)
        ENDIF
      ENDIF
!
!---  Injection well is inactive set well pressure to be in 
!     equilibrium with reservoir  ---
!
      IF( TMZ.LE.VAR_CW(1,1,NCW) ) THEN
        CALL EQUIL_COUP_WELL( NCW )
        ID_CW(8,NCW) = 1
        GOTO 500
      ENDIF
!
!---  Loop over the injection well time periods  ---
!
      NS = 1
      DO 120 NTP = 1,IM_CW(NCW)
!
!---    Injection well time period only has one time (start time)  ---
!
        IF( IMP_CW(NTP,NCW).EQ.1 ) THEN
!
!---      Time prior to start time, injection well is inactive, set 
!         well pressure to be in equilibrium with reservoir  ---
!
          IF( TMZ.LE.VAR_CW(1,NS,NCW) ) THEN
            CALL EQUIL_COUP_WELL( NCW )
            ID_CW(8,NCW) = 1
            GOTO 500
!
!---      Time after start time, injection well is active  ---
!
          ELSE
            DO 80 N = 2,NGC+8
              VAR_CWX(N) = VAR_CW(N,1,NCW)
   80       CONTINUE
            QMO_CWX = 0.D+0
            GOTO 122
          ENDIF
!
!---    Injection well time period has multiple times  ---
!
        ELSE
          NE = NS + IMP_CW(NTP,NCW) - 1
!
!---      Time outside of injection well time period, go to next 
!         injection well time period  ---
!
          IF( TMZ.LE.VAR_CW(1,NS,NCW) .OR. 
     &      TMZ.GT.VAR_CW(1,NE,NCW) ) GOTO 110 
!
!---      Injection well time period has multiple time points, use  
!         linear interpolation of well parameters between time points  ---
!
          DO 100 M = 2,IMP_CW(NTP,NCW)
            MX = NS + M - 1
            IF( TMZ.LE.VAR_CW(1,MX,NCW) ) THEN
              TD_CW = VAR_CW(1,MX,NCW)-VAR_CW(1,MX-1,NCW)
              DT_CW = MIN( VAR_CW(1,MX,NCW)-TMZ,DT )
              TF_CW = (TMZ-VAR_CW(1,MX-1,NCW))/TD_CW
              DO 90 N = 2,NGC+8
                VAR_CWX(N) = VAR_CW(N,MX-1,NCW) + 
     &            TF_CW*(VAR_CW(N,MX,NCW)-VAR_CW(N,MX-1,NCW))
   90         CONTINUE
              QMO_CWX = 0.D+0
              GOTO 122
            ENDIF
  100     CONTINUE
        ENDIF
  110   CONTINUE
        NS = NS + IMP_CW(NTP,NCW)
  120 CONTINUE
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
!---  Pressure controlled well ---
!
      IF( PL_CW(NCW)-P_CW(2,NCW).LT.EPSL ) THEN
        ID_CW(8,NCW) = 1
      ENDIF
!
!---  Excessive flow rate, pressure controlled well ---
!
      IF( VAR_CWX(2).GT.1.D+5 ) THEN
        ID_CW(8,NCW) = 1
        P_CW(2,NCW) = PL_CW(NCW)
      ENDIF
!
!---  Loop over increment indices ---
!
      DO 300 M = 1,ISVC+2
        MW = MCW(M)
        MF = MFD(M)
        N = IWN_CW(ID_CW(3,NCW))
        P_CWX = P_CW(MW,NCW)
        NPHAZX = 0
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          T_CWX = VAR_CWX(4)
        ELSE
          T_CWX = T(2,N)
        ENDIF
!
!---    Well state #1  ---
!
!       SL = 1.0
!       SG+SN = 0.0
!
!       Declared variables:
!
!       VAR_CWX(5) - aqueous CO2 relative saturation, or
!                    aqueous CO2 mass fraction
!       VAR_CWX(6) - aqueous salt relative saturation, or
!                         aqueous salt mass fraction
!
        IF( ITS_CW(NTP,NCW)/100.EQ.1 ) THEN
          PVAX = VAR_CWX(5)
          TMSX = VAR_CWX(6)
          ISR4X = MOD(ITS_CW(NTP,NCW),10)
          ISR5X = MOD(ITS_CW(NTP,NCW)/10,10)
          CALL FLH_IC1( P_CWX,PVAX,RHOL_CWX,RHOMLX,TMSX,T_CWX,XLA_CWX,
     &      XLS_CWX,XLW_CWX,XMLAX,XMLSX,XMLWX,YLSX,ZMCX,ISR4X,ISR5X )
          SG_CWX = 0.D+0
          SN_CWX = 0.D+0
          SL_CWX = 1.D+0
          RHOF_CWX = RHOL_CWX
!
!---    Well state #2  ---
!
!       SL < 1.0
!       SG+SN > 0.0
!
!       Declared variables:
!
!       VAR_CWX(5) - aqueous saturation
!       VAR_CWX(6) - aqueous salt relative saturation, or
!                    aqueous salt mass fraction
!       VAR_CWX(6+NGC) - total mole fraction of petroleum comp.
!
        ELSEIF( ITS_CW(NTP,NCW)/100.EQ.2 ) THEN
          SL_CWX = VAR_CWX(5)
          TMSX = VAR_CWX(6)
          ISR5X = MOD(ISRT(NS)/10,10)
          DO 130 IGC = 1,NGC+2
            JGC = 6 + IGC
            ZMCX(IGC) = VAR_CWX(JGC)
  130     CONTINUE
          IZNX = 0
          CALL FLH_IC2( P_CWX,RHOG_CWX,RHOMGX,RHOL_CWX,RHOMLX,
     &      RHON_CWX,RHOMNX,SG_CWX,SL_CWX,SN_CWX,TMSX,T_CWX,
     &      XGC_CWX,
     &      XMGCX,XGW_CWX,XMGWX,XLA_CWX,XLS_CWX,XLW_CWX,XMLAX,XMLSX,
     &      XMLWX,XNC_CWX,XMNCX,YLSX,ZGX,ZMCX,ZNX,ISR5X,IZNX,NPHAZX )
          RHOF_CWX = SG_CWX*RHOG_CWX + SL_CWX*RHOL_CWX + 
     &      SN_CWX*RHON_CWX
!
!---    Well state #3  ---
!
!       SL = 0.0
!       SG+SN = 1.0
!
!       Declared variables:
!
!       VAR_CWX(5) - water vapor relative saturation
!       VAR_CWX(6+) - total mole fraction of petroleum comp.
!      
        ELSEIF( ITS_CW(NTP,NCW)/100.EQ.3 ) THEN
          DO 132 IGC = 1,NGC+2
            JGC = 6 + IGC
            ZMCX(IGC) = VAR_CWX(JGC)
  132     CONTINUE
          SL_CWX = 0.D+0
          PVWX = VAR_CWX(5)
          YLSX = 0.D+0
          IZNX = 0
          CALL FLH_IC3( PVWX,P_CWX,RHOG_CWX,RHOMGX,RHON_CWX,RHOMNX,
     &      SG_CWX,SN_CWX,T_CWX,XGC_CWX,XGW_CWX,XMGCX,
     &      XMGWX,XNC_CWX,XMNCX,YLSX,ZGX,ZMCX,ZNX,IZNX,NPHAZX )
          RHOF_CWX = SG_CWX*RHOG_CWX + SN_CWX*RHON_CWX
        ENDIF
!
!---    Load local variable for mass flux, kg/s  ---
!
        IF( M.EQ.1 ) THEN
!
!---      Convert volumetric flow rate, m^3/s to mass flow rate, kg/s,
!         using fluid density at top of screened interval  ---
!
          IF( IT_CW(NCW).EQ.2 ) THEN
            VAR_CW2X = VAR_CWX(2)*RHOF_CWX
          ELSE
            VAR_CW2X = VAR_CWX(2)
          ENDIF
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
          NPHAZX = 0
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            T_CWX = VAR_CWX(4)
          ELSE
            T_CWX = T(2,N)
          ENDIF
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
          P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOF_CWX
!
!---      Well pressure at the node centroid, used for coupled-well
!         nodal output  ---
!
          IF( M.EQ.1 ) THEN
            NWF = IWP_CW(NWN)
            PF_CW(NWF) = P_CWX - (ZPNX-ZPX(1))*GRAV*RHOF_CWX
          ENDIF
!
!---      Adjust the reservoir pressure to the coupled-well node
!         centroid  ---
!
!---      Phase condition #1  ---
!
!         SL = 1.0, SG+SN = 0.0
!
          IF( MOD(NPHAZ(2,N),10).EQ.1 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #2  ---
!
!         SL < 1.0, ST > 0.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.2 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #3  ---
!
!         SL < 1.0, ST > 0.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.3 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #4  ---
!
!         SL = 0.0, ST = 1.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.4 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #5  ---
!
!         SL = 0.0, ST = 1.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.5 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
          ENDIF
!
!---      Well state #1  ---
!
!         SL = 1.0
!         SG+SN = 0.0
!
!         Declared variables:
!
!         VAR_CWX(5) - aqueous CO2 relative saturation, or
!                      aqueous CO2 mass fraction
!         VAR_CWX(6) - aqueous salt relative saturation, or
!                           aqueous salt mass fraction
!
          IF( ITS_CW(NTP,NCW)/100.EQ.1 ) THEN
            PVAX = VAR_CWX(5)
            TMSX = VAR_CWX(6)
            ISR4X = MOD(ITS_CW(NTP,NCW),10)
            ISR5X = MOD(ITS_CW(NTP,NCW)/10,10)
            CALL FLH_IC1( P_CWX,PVAX,RHOL_CWX,RHOMLX,TMSX,T_CWX,XLA_CWX,
     &      XLS_CWX,XLW_CWX,XMLAX,XMLSX,XMLWX,YLSX,ZMCX,ISR4X,ISR5X )
            SG_CWX = 0.D+0
            SN_CWX = 0.D+0
            SL_CWX = 1.D+0
            RHOF_CWX = RHOL_CWX
!
!---      Well state #2  ---
!
!         SL < 1.0
!         SG+SN > 0.0
!
!         Declared variables:
!
!         VAR_CWX(5) - aqueous saturation
!         VAR_CWX(6) - aqueous salt relative saturation, or
!                      aqueous salt mass fraction
!         VAR_CWX(6+NGC) - total mole fraction of petroleum comp.
!
          ELSEIF( ITS_CW(NTP,NCW)/100.EQ.2 ) THEN
            SL_CWX = VAR_CWX(5)
            TMSX = VAR_CWX(6)
            ISR5X = MOD(ISRT(NS)/10,10)
            DO 140 IGC = 1,NGC+2
              JGC = 6 + IGC
              ZMCX(IGC) = VAR_CWX(JGC)
  140       CONTINUE
            IZNX = 0
            CALL FLH_IC2( P_CWX,RHOG_CWX,RHOMGX,RHOL_CWX,RHOMLX,
     &        RHON_CWX,RHOMNX,SG_CWX,SL_CWX,SN_CWX,TMSX,T_CWX,XGC_CWX,
     &        XMGCX,XGW_CWX,XMGWX,XLA_CWX,XLS_CWX,XLW_CWX,XMLAX,XMLSX,
     &        XMLWX,XNC_CWX,XMNCX,YLSX,ZGX,ZMCX,ZNX,ISR5X,IZNX,NPHAZX )
            RHOF_CWX = SG_CWX*RHOG_CWX + SL_CWX*RHOL_CWX + 
     &        SN_CWX*RHON_CWX
!
!---      Well state #3  ---
!
!         SL = 0.0
!         SG+SN = 1.0
!
!         Declared variables:
!
!         VAR_CWX(5) - water vapor relative saturation
!         VAR_CWX(6+) - total mole fraction of petroleum comp.
!        
          ELSEIF( ITS_CW(NTP,NCW)/100.EQ.3 ) THEN
            DO 142 IGC = 1,NGC+2
              JGC = 6 + IGC
              ZMCX(IGC) = VAR_CWX(JGC)
  142       CONTINUE
            SL_CWX = 0.D+0
            PVWX = VAR_CWX(5)
            YLSX = 0.D+0
            IZNX = 0
            CALL FLH_IC3( PVWX,P_CWX,RHOG_CWX,RHOMGX,RHON_CWX,RHOMNX,
     &        SG_CWX,SN_CWX,T_CWX,XGC_CWX,XGW_CWX,XMGCX,XMGWX,
     &        XNC_CWX,
     &        XMNCX,YLSX,ZGX,ZMCX,ZNX,IZNX,NPHAZX )
            RHOF_CWX = SG_CWX*RHOG_CWX + SN_CWX*RHON_CWX
          ENDIF
!
!---      Well gas viscosity, Pa s  ---
!
          IF( SG_CWX.GT.0.D+0 ) THEN
            CALL VISC_G( RHOMGX,T_CWX,VISG_CWX,XMGCX,XMGWX )
          ENDIF
!
!---      Well aqueous viscosity, Pa s
!         concentration of CO2 in mol/kg,
!         concentration of salt in mol/kg  ---
!
          IF( SL_CWX.GT.0.D+0 ) THEN
            CMLAX = 0.D+0
            CMLSX = 1.D+3*XLS_CWX*RHOL_CWX/WTMS
            CALL VISC_L( CMLAX,CMLSX,P_CWX,T_CWX,VISL_CWX )
          ENDIF
!
!---      Well nonaqueous-liquid viscosity, Pa s  ---
!
          IF( SN_CWX.GT.0.D+0 ) THEN
            CALL VISC_N( RHOMNX,T_CWX,VISN_CWX,XMNCX )
          ENDIF
!
!---      Fluid enthalpy, unless nonisothermal simulation  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Well gas enthalpy and internal energy, W/kg  ---
!
            IF( SG_CWX.GT.0.D+0 ) THEN
              CALL ENTH_G( HG_CWX,HGA_CWX,HGW_CWX,P_CWX,T_CWX,
     &          RHOG_CWX,UEGX,XGC_CWX,XGW_CWX,XMGCX,ZGX )
            ENDIF
!
!---        Well nonaqueous-liquid enthalpy and 
!           internal energy, W/kg  ---
!
            IF( SN_CWX.GT.0.D+0 ) THEN
              CALL ENTH_N( HN_CWX,HNAX,P_CWX,T_CWX,RHON_CWX,
     &           UENX,XNC_CWX,XMNCX,ZNX )
            ENDIF
!
!---        Well aqueous enthalpy and internal energy, W/kg  ---
!
            IF( SL_CWX.GT.0.D+0 ) THEN
!
!---          Pure water liquid enthalpy  ---
!
              CALL SP_W( T_CWX,PSWX )
              PEX = MAX( PSWX,P_CWX )
              ISRX = 1
              CALL ENTH_W( T_CWX,PEX,HLWX,HGWX,ISRX )
!
!---          Brine enthalpy  ---
!
              CALL ENTH_B( T_CWX,XLS_CWX,HLWX,HBX )
!
!---          Dissolved CO2 contribution  ---
!
              IF( SG_CWX.GT.0.D+0 ) THEN
                FGAX = SG_CWX*RHOG_CWX*XGC_CWX(1)
              ELSE
                FGAX = 0.D+0
              ENDIF
              IF( SN_CWX.GT.0.D+0 ) THEN
                FNAX = SN_CWX*RHON_CWX*XNC_CWX(1)
              ELSE
                FNAX = 0.D+0
              ENDIF
              IF( (FGAX+FNAX)/EPSL.LT.EPSL ) THEN
                HAX = 0.D+0
              ELSE
                HAX = (HGA_CWX*FGAX + HNAX*FNAX)/(FGAX + FNAX)
              ENDIF
!
!---          Aqueous enthalpy and internal energy  ---
!
              CALL ENTH_L( HAX,HBX,HL_CWX,T_CWX,XLA_CWX,XLS_CWX )
            ENDIF
          ENDIF
!
!---      Equivalent field node radius components  ---
!
          PERMX = MAX( PERM(1,IZN),1.D-20 )
          PERMY = MAX( PERM(2,IZN),1.D-20 )
          PERMZ = MAX( PERM(3,IZN),1.D-20 )
          RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---      Cylindrical coordinates with azimuthal symmetry,
!         centrally located wells  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            ROZ = RP(I)
!
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROZ = PAR_CW(3,INVX)
            ENDIF
            RWX = MIN( RWX,9.999D-1*ROZ )
            PERMX = PERMRF(MF,N)*PERM(1,IZN)
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
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROX = PAR_CW(3,INVX)
              ROY = PAR_CW(3,INVX)
              ROZ = PAR_CW(3,INVX)
            ENDIF
!
!---        Well index components  ---
!
            PERMX = PERMRF(MF,N)*PERM(1,IZN)
            PERMY = PERMRF(MF,N)*PERM(2,IZN)
            PERMZ = PERMRF(MF,N)*PERM(3,IZN)
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
          DPNX = MAX( P_CWX-PGFX,0.D+0 )
          DPLX = MAX( P_CWX-PGFX,0.D+0 )
!
!---      Zero fluxes from well to reservoir  ---
!
          FXE_CW(M,NWN) = 0.D+0
          FXW_CW(M,NWN) = 0.D+0
          DO IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = 0.D+0
          ENDDO
!
!---      Gas flux from well to reservoir  ---
!
          IF( SG_CWX.GT.0.D+0 ) THEN
            FGX_CWX = SG_CWX*WI_CWX*RHOG_CWX*DPGX/VISG_CWX
            IF( ISLC(30).EQ.0 )
     &        FXE_CW(M,NWN) = FXE_CW(M,NWN) + FGX_CWX*HG_CWX
            FXW_CW(M,NWN) = FXW_CW(M,NWN) + FGX_CWX*XGW_CWX
            DO IGC = 1,NGC+2
              IF( ISLC(45).NE.0 .AND. IGC.EQ.1 ) CYCLE
              IF( ISLC(64).NE.0 .AND. IGC.EQ.2 ) CYCLE
              FXC_CW(IGC,M,NWN) = FXC_CW(IGC,M,NWN) + 
     &          FGX_CWX*XGC_CWX(IGC)
            ENDDO
          ELSE
            FGX_CWX = 0.D+0
          ENDIF
!
!---      Aqueous flux from well to reservoir  ---
!
          IF( SL_CWX.GT.0.D+0 ) THEN
            FLX_CWX = SL_CWX*WI_CWX*RHOL_CWX*DPLX/VISL_CWX
            IF( ISLC(30).EQ.0 )
     &        FXE_CW(M,NWN) = FXE_CW(M,NWN) + FLX_CWX*HL_CWX
            FXW_CW(M,NWN) = FXW_CW(M,NWN) + FLX_CWX*XLW_CWX
            IF( ISLC(45).EQ.0 ) 
     &        FXC_CW(1,M,NWN) = FXC_CW(1,M,NWN) + FLX_CWX*XLA_CWX
          ELSE
            FLX_CWX = 0.D+0
          ENDIF
!
!---      Nonaqueous-liquid flux from well to reservoir  ---
!
          IF( SN_CWX.GT.0.D+0 ) THEN
            FNX_CWX = SN_CWX*WI_CWX*RHON_CWX*DPNX/VISN_CWX
            IF( ISLC(30).EQ.0 )
     &        FXE_CW(M,NWN) = FXE_CW(M,NWN) + FNX_CWX*HN_CWX
            DO IGC = 1,NGC+2
              IF( ISLC(45).NE.0 .AND. IGC.EQ.1 ) CYCLE
              IF( ISLC(64).NE.0 .AND. IGC.EQ.2 ) CYCLE
              FXC_CW(IGC,M,NWN) = FXC_CW(IGC,M,NWN) + 
     &          FNX_CWX*XNC_CWX(IGC)
            ENDDO
          ELSE
            FNX_CWX = 0.D+0
          ENDIF
!
!---      Store current coupled-well node location in previous
!         coupled-well node location  ---
!
          XPX(1) = XPX(2)
          YPX(1) = YPX(2)
          ZPX(1) = ZPX(2)
!
!---      Volumetric injection well fluxes  ---
!
          IF( M.EQ.1 ) THEN
            IF( SL_CWX.GT.0.D+0 ) Q_CW(2,NWN) = FLX_CWX/RHOL_CWX
            IF( SG_CWX.GT.0.D+0 ) Q_CW(3,NWN) = FGX_CWX/RHOG_CWX
            IF( SN_CWX.GT.0.D+0 ) Q_CW(4,NWN) = FNX_CWX/RHON_CWX
            Q_CW(1,NWN) = Q_CW(2,NWN) + Q_CW(3,NWN) + Q_CW(4,NWN)
          ENDIF
  200   CONTINUE
  300 CONTINUE
!
!---  Mass balance residuals for injection type coupled well  ---
!
      NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
      MX = (NWFX*ISVC)+2
      RS_CW(1,NCW) = VAR_CW2X
      RS_CW(MX,NCW) = VAR_CW2X
      QM_CW(1,NCW) = 0.D+0
      QM_CW(3,NCW) = 0.D+0
      QM_CW(5,NCW) = 0.D+0
      DO IGC = 1,NGC+2
        QM_CW((7+(IGC-1)*2),NCW) = 0.D+0
      ENDDO
      QM_CWX = 0.D+0
      QTM_CWX = 0.D+0
!
!---  Loop over coupled-well nodes  ---
!
      DO 330 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
        RS_CW(1,NCW) = RS_CW(1,NCW) - FXW_CW(1,NWN)
        RS_CW(MX,NCW) = RS_CW(MX,NCW) - FXW_CW(ISVC+2,NWN)
        QTM_CWX = QTM_CWX + FXW_CW(1,NWN)
        QM_CW(3,NCW) = QM_CW(3,NCW) + FXW_CW(1,NWN)
        IF( ISLC(30).EQ.0 )
     &    QM_CW(5,NCW) = QM_CW(5,NCW) + FXE_CW(1,NWN)
        QM_CWX = QM_CWX + FXW_CW(ISVC+2,NWN)
        DO 320 IGC = 1,NGC+2
          RS_CW(1,NCW) = RS_CW(1,NCW) - FXC_CW(IGC,1,NWN)
          RS_CW(MX,NCW) = RS_CW(MX,NCW) - FXC_CW(IGC,ISVC+2,NWN)
          QM_CW(1,NCW) = QM_CW(1,NCW) + FXC_CW(IGC,1,NWN)
          QTM_CWX = QTM_CWX + FXC_CW(IGC,1,NWN)
          QM_CWX = QM_CWX + FXC_CW(IGC,ISVC+2,NWN)
          QM_CW((7+(IGC-1)*2),NCW) = QM_CW((7+(IGC-1)*2),NCW) + 
     &      FXC_CW(IGC,1,NWN)
  320   CONTINUE
  330 CONTINUE
!
!---  Hold pressure controlled option  ---
!
      IF( ID_CW(8,NCW).EQ.1 .AND. NITER.GE.3 ) THEN
        P_CW(2,NCW) = PL_CW(NCW)
        DNR_CW(NCW) = 1.D-1
        P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
        ID_CW(8,NCW) = 1
        GOTO 500
!
!---  Insufficient increment in coupled-well pressure to create
!     flow from well, increase increment and well pressure  ---
!
      ELSEIF( ABS(QM_CWX/DNR_CW(NCW)).LT.1.D-7 ) THEN
        DNR_CW(NCW) = 1.25D+0*DNR_CW(NCW)
        P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
        ID_CW(8,NCW) = 0
!
!---    Incremented pressure exceeds upper pressure limit,
!       transition to pressure-controlled well  ---
!
        IF( P_CW(3,NCW).GE.PL_CW(NCW) ) THEN
          P_CW(2,NCW) = PL_CW(NCW)
          DNR_CW(NCW) = 1.D-1
          P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
          ID_CW(8,NCW) = 1
          GOTO 500
        ENDIF
!
!---    Zero coupled-well fluxes ---
!
        QM_CW(1,NCW) = 0.D+0
        QM_CW(3,NCW) = 0.D+0
        QM_CW(5,NCW) = 0.D+0
        DO IGC = 1,NGC+2
          QM_CW((7+(IGC-1)*2),NCW) = 0.D+0
        ENDDO
!
!---    Loop over coupled-well nodes  ---
!
        DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---      Zero fluxes from well to reservoir  ---
!
          FXE_CW(M,NWN) = 0.D+0
          FXW_CW(M,NWN) = 0.D+0
          DO IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = 0.D+0
          ENDDO
        ENDDO
!
!---    Recompute well fluxes at higher well pressure increment ---
!
        QMO_CWX = QM_CWX
        GOTO 122
!
!---  Excessive well flow for pressure controlled well,
!     reduce well pressure  ---
!
      ELSEIF( ABS(ID_CW(8,NCW)).EQ.1 .AND. QM_CWX.GT.VAR_CW2X ) THEN
        DQ_CWX = 1.D+1*DQ_CWX
        DP_CWX = -DQ_CWX*PATM
        P_CW(2,NCW) = P_CW(2,NCW) + DP_CWX
        P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
!
!---    Zero coupled-well fluxes ---
!
        QM_CW(1,NCW) = 0.D+0
        QM_CW(3,NCW) = 0.D+0
        QM_CW(5,NCW) = 0.D+0
        DO IGC = 1,NGC+2
          QM_CW((7+(IGC-1)*2),NCW) = 0.D+0
        ENDDO
!
!---    Pressure correction for pressure controlled well ---
!
        ID_CW(8,NCW) = -1
!
!---    Loop over coupled-well nodes  ---
!
        DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---      Zero fluxes from well to reservoir  ---
!
          FXE_CW(M,NWN) = 0.D+0
          FXW_CW(M,NWN) = 0.D+0
          DO IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = 0.D+0
          ENDDO
        ENDDO
!
!---    Recompute well fluxes lower well pressure ---
!
        QMO_CWX = QM_CWX
        GOTO 122
!
!---  Acceptable pressure reduction for pressure controlled well, 
!     switch to flow controlled well  ---
!
      ELSEIF( ID_CW(8,NCW).EQ.-1 .AND. QM_CWX.LE.VAR_CW2X ) THEN
        ID_CW(8,NCW) = 0
      ENDIF
!
!---  Loop over field nodes that contain coupled-well nodes  ---
!
      DO 470 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
        M1 = (NWF-ID_CW(5,NCW))*ISVC + 1
        DO 460 M2 = 1,ISVC
          M3 = M1+M2
          RS_CW(M3,NCW) = VAR_CW2X
!
!---      Loop over coupled-well nodes  ---
!
          DO 450 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---        If coupled-well node is within the current field
!           node, use incremented fluxes  ---
!
            IF( IWF_CW(NWF).EQ.IWN_CW(NWN) ) THEN
              RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXW_CW(M2+1,NWN)                
              DO 440 IGC = 1,NGC+2
                RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXC_CW(IGC,M2+1,NWN)
  440         CONTINUE
!
!---        If coupled-well node is outside the current field
!           node, use un-incremented fluxes  ---
!
            ELSE
              RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXW_CW(1,NWN)
              DO 442 IGC = 1,NGC+2
                RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXC_CW(IGC,1,NWN)
  442         CONTINUE
            ENDIF
  450     CONTINUE
  460   CONTINUE
  470 CONTINUE
  500 CONTINUE
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
      SUBROUTINE INJP_BO_WELL( NCW )
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
!     Injection coupled well model
!     
!     Rate controlled or pressure controlled
!
!     Flux of energy, water mass, CO2 mass, CH4 mass, and petroleum
!     component masses from coupled-well nodes to field nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 April 2014.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE HYST
      USE GRID
      USE GLB_PAR
      USE FDVP
      USE COUP_WELL
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XPX(2),YPX(2),ZPX(2)
      REAL*8 VAR_CWX(6+LNGC)
      REAL*8 XGC_CWX(LNGC),XNC_CWX(LNGC)
      REAL*8 FXC_CWX(LNGC,3,LWN_CW),FXW_CWX(3,LWN_CW)
      REAL*8 QM_CWX(3),P_CWY(3),VAR_CW2X(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INJP_BO_WELL'
!
!---  Loop over coupled-well nodes  ---
!
      DO 30 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---    Loop over increment indices  ---
!
        DO 20 M = 1,2
!
!---      Loop over petroleum components  ---
!
          DO 12 IGC = 1,NGC+2
            FXC_CWX(IGC,M,NWN) = 0.D+0
   12     CONTINUE
          FXW_CWX(M,NWN) = 0.D+0
   20   CONTINUE
   30  CONTINUE
!
!---  Injection well time interval ---
!
      TMZ = TM
      IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---  Cyclic time periods  ---
!
      IF( ICC_CW(NCW).EQ.1 ) THEN
!
!---    Loop over the injection well time periods, to find the 
!       final well time  ---
!
        NTX = 0
        DO 40 NTP = 1,IM_CW(NCW)
          NTX = NTX + IMP_CW(NTP,NCW)
   40   CONTINUE
!
!---    Determine time with the cyclic time period  ---
!
        TMZ = MOD( TM,VAR_CW(1,NTX,NCW) )
        IF( TM.GT.VAR_CW(1,NTX,NCW) ) THEN
          IF( TMZ.LT.EPSL ) TMZ = VAR_CW(1,NTX,NCW)
        ENDIF
      ENDIF
!
!---  Injection well is inactive set well pressure to be in 
!     equilibrium with reservoir  ---
!
      IF( TMZ.LE.VAR_CW(1,1,NCW) ) THEN
        CALL EQUIL_COUP_WELL( NCW )
        ID_CW(8,NCW) = 1
        GOTO 500
      ENDIF
!
!---  Loop over the injection well time periods  ---
!
      NS = 1
      DO 120 NTP = 1,IM_CW(NCW)
!
!---    Injection well time period only has one time (start time)  ---
!
        IF( IMP_CW(NTP,NCW).EQ.1 ) THEN
!
!---      Time prior to start time, injection well is inactive, set 
!         well pressure to be in equilibrium with reservoir  ---
!
          IF( TMZ.LE.VAR_CW(1,NS,NCW) ) THEN
            CALL EQUIL_COUP_WELL( NCW )
            ID_CW(8,NCW) = 1
            GOTO 500
!
!---      Time after start time, injection well is active  ---
!
          ELSE
            DO 80 N = 2,NGC+8
              VAR_CWX(N) = VAR_CW(N,1,NCW)
   80       CONTINUE
            GOTO 122
          ENDIF
!
!---    Injection well time period has multiple times  ---
!
        ELSE
          NE = NS + IMP_CW(NTP,NCW) - 1
!
!---      Time outside of injection well time period, go to next 
!         injection well time period  ---
!
          IF( TMZ.LE.VAR_CW(1,NS,NCW) .OR. 
     &      TMZ.GT.VAR_CW(1,NE,NCW) ) GOTO 110 
!
!---      Injection well time period has multiple time points, use  
!         linear interpolation of well parameters between time points  ---
!
          DO 100 M = 2,IMP_CW(NTP,NCW)
            MX = NS + M - 1
            IF( TMZ.LE.VAR_CW(1,MX,NCW) ) THEN
              TD_CW = VAR_CW(1,MX,NCW)-VAR_CW(1,MX-1,NCW)
              DT_CW = MIN( VAR_CW(1,MX,NCW)-TMZ,DT )
              TF_CW = (TMZ-VAR_CW(1,MX-1,NCW))/TD_CW
              DO 90 N = 2,NGC+8
                VAR_CWX(N) = VAR_CW(N,MX-1,NCW) + 
     &            TF_CW*(VAR_CW(N,MX,NCW)-VAR_CW(N,MX-1,NCW))
   90         CONTINUE
              GOTO 122
            ENDIF
  100     CONTINUE
        ENDIF
  110   CONTINUE
        NS = NS + IMP_CW(NTP,NCW)
  120 CONTINUE
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
        NPHAZX = 0
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          T_CWX = VAR_CWX(4)
        ELSE
          T_CWX = T(2,N)
        ENDIF
!
!---    Well state #1  ---
!
!       SL = 1.0
!       SG+SN = 0.0
!
!       Declared variables:
!
!       VAR_CWX(6) - aqueous salt relative saturation, or
!                    aqueous salt mass fraction
!
        IF( ITS_CW(NTP,NCW)/100.EQ.1 ) THEN
          PVAX = VAR_CWX(5)
          TMSX = VAR_CWX(6)
          ISR5X = MOD(ITS_CW(NTP,NCW)/10,10)
          CALL BO_IC1( P_CWX,RHOL_CWX,TMSX,T_CWX,
     &      XLS_CWX,XLW_CWX,YLSX,ISR5X )
          SG_CWX = 0.D+0
          SN_CWX = 0.D+0
          SL_CWX = 1.D+0
          RHOF_CWX = RHOL_CWX
!
!---    Well state #2  ---
!
!       SL < 1.0
!       SG+SN > 0.0
!
!       Declared variables:
!
!       VAR_CWX(5) - aqueous saturation
!       VAR_CWX(6) - aqueous salt relative saturation, or
!                    aqueous salt mass fraction
!       VAR_CWX(7) - gas mass fraction of oil + gas
!
        ELSEIF( ITS_CW(NTP,NCW)/100.EQ.2 ) THEN
          SL_CWX = VAR_CWX(5)
          TMSX = VAR_CWX(6)
          XGX = VAR_CWX(7)
          RHOGSX = 1.225D+0*GCPP(37,1)
          RHONSX = 1.D+3*GCPP(37,2)
          PX = P_CWX + PATM
          IF( 1.D+0-XGX.GT.EPSL ) THEN
            POSM_CWX = XGX*RHONSX/((1.D+0-XGX)*RHOGSX)
            CALL BO_BP( PBPX,POSM_CWX,T_CWX )
            IF( PX.GT.PBPX ) THEN
              SN_CWX = 1.D+0
              SG_CWX = 0.D+0
            ELSE
              PVWX = 0.D+0
              CALL BO_EQ_GAS( PX,PVWX,RHOGX,T_CWX,XGC_CWX,XGWX )
              CALL BO_EQ_OIL( RHONX,PX,POSM_CWX,T_CWX,XNC_CWX )
              SN_CWX = RHOGX*(XGX*(1.D+0-SL_CWX)-XGC_CWX(1)*
     &          (1.D+0-SL_CWX))/(RHONX*(XNC_CWX(1)-XGX) - 
     &          RHOGX*(XGC_CWX(1)-XGX))
              SG_CWX = MAX( 1.D+0-SL_CWX-SN_CWX,0.D+0 )
            ENDIF
          ELSE
            SG_CWX = 1.D+0
            SN_CWX = 0.D+0
          ENDIF
          ISR5X = MOD(ISRT(NS)/10,10)
          IZNX = 0
          CALL BO_IC2( POSM_CWX,P_CWX,RHOG_CWX,RHOL_CWX,RHON_CWX,SL_CWX,
     &      TMSX,T_CWX,XGC_CWX,XGW_CWX,XLS_CWX,XLW_CWX,XNC_CWX,
     &      YLSX,ISR5X,IZNX )
          RHOF_CWX = SG_CWX*RHOG_CWX + SL_CWX*RHOL_CWX + 
     &      SN_CWX*RHON_CWX
!
!---    Well state #3  ---
!
!       SL = 0.0
!       SG+SN = 1.0
!
!       Declared variables:
!
!       VAR_CWX(5) - water vapor relative saturation
!       VAR_CWX(7) - gas mass fraction of oil + gas
!      
        ELSEIF( ITS_CW(NTP,NCW)/100.EQ.3 ) THEN
          SL_CWX = 0.D+0
          PVWX = VAR_CWX(5)
          YLSX = 0.D+0
          XGX = VAR_CWX(7)
          RHOGSX = 1.225D+0*GCPP(37,1)
          RHONSX = 1.D+3*GCPP(37,2)
          PX = P_CWX + PATM
          IF( 1.D+0-XGX.GT.EPSL ) THEN
            POSM_CWX = XGX*RHONSX/((1.D+0-XGX)*RHOGSX)
            CALL BO_BP( PBPX,POSM_CWX,T_CWX )
            IF( PX.GT.PBPX ) THEN
              SN_CWX = 1.D+0
              SG_CWX = 0.D+0
            ELSE
              PVWX = 0.D+0
              CALL BO_EQ_GAS( PX,PVWX,RHOGX,T_CWX,XGC_CWX,XGWX )
              CALL BO_EQ_OIL( RHONX,PX,POSM_CWX,T_CWX,XNC_CWX )
              SN_CWX = RHOGX*(XGX*(1.D+0-SL_CWX)-XGC_CWX(1)*
     &          (1.D+0-SL_CWX))/(RHONX*(XNC_CWX(1)-XGX) - 
     &          RHOGX*(XGC_CWX(1)-XGX))
              SG_CWX = MAX( 1.D+0-SL_CWX-SN_CWX,0.D+0 )
            ENDIF
          ELSE
            SG_CWX = 1.D+0
            SN_CWX = 0.D+0
          ENDIF
          IZNX = 0
          CALL BO_IC3( POSM_CWX,PVWX,P_CWX,RHOG_CWX,RHON_CWX,
     &      T_CWX,XGC_CWX,XGW_CWX,XNC_CWX,YLSX,IZNX )
          RHOF_CWX = SG_CWX*RHOG_CWX + SN_CWX*RHON_CWX
        ENDIF
!
!---    Convert volumetric flow rate, m^3/s to mass flow rate, kg/s,
!       using fluid density at top of screened interval  ---
!
        IF( IT_CW(NCW).EQ.2 ) THEN
          VAR_CW2X(M) = VAR_CWX(2)*RHOF_CWX
        ELSE
          VAR_CW2X(M) = VAR_CWX(2)
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
          NPHAZX = 0
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            T_CWX = VAR_CWX(4)
          ELSE
            T_CWX = T(2,N)
          ENDIF
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
          P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOF_CWX
!
!---      Adjust the reservoir pressure to the coupled-well node
!         centroid  ---
!
!---      Phase condition #1  ---
!
!         SL = 1.0, SG+SN = 0.0
!
          IF( MOD(NPHAZ(2,N),10).EQ.1 ) THEN
            PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
            PNFX = PN(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(2,N)-MAX( PL(2,N),PN(2,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #2  ---
!
!         SL < 1.0, ST > 0.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.2 ) THEN
            PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
            PNFX = PN(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(2,N)-MAX( PL(2,N),PN(2,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #3  ---
!
!         SL < 1.0, ST > 0.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.3 ) THEN
            PGFX = PG(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(2,N)
            PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
            PNFX = PN(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
!
!---      Phase condition #4  ---
!
!         SL = 0.0, ST = 1.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.4 ) THEN
            PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
            PNFX = PN(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(2,N)-MAX( PL(2,N),PN(2,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #5  ---
!
!         SL = 0.0, ST = 1.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.5 ) THEN
            PGFX = PG(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(2,N)
            PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
            PNFX = PN(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
          ENDIF
!
!---      Well state #1  ---
!
!         SL = 1.0
!         SG+SN = 0.0
!
!         Declared variables:
!
!         VAR_CWX(6) - aqueous salt relative saturation, or
!                           aqueous salt mass fraction
!
          IF( ITS_CW(NTP,NCW)/100.EQ.1 ) THEN
            PVAX = VAR_CWX(5)
            TMSX = VAR_CWX(6)
            ISR5X = MOD(ITS_CW(NTP,NCW)/10,10)
            CALL BO_IC1( P_CWX,RHOL_CWX,TMSX,T_CWX,
     &        XLS_CWX,XLW_CWX,YLSX,ISR5X )
            SG_CWX = 0.D+0
            SN_CWX = 0.D+0
            SL_CWX = 1.D+0
            RHOF_CWX = RHOL_CWX
!
!---      Well state #2  ---
!
!         SL < 1.0
!         SG+SN > 0.0
!
!         Declared variables:
!
!         VAR_CWX(5) - aqueous saturation
!         VAR_CWX(6) - aqueous salt relative saturation, or
!                      aqueous salt mass fraction
!         VAR_CWX(7) - gas mass fraction of oil + gas
!
          ELSEIF( ITS_CW(NTP,NCW)/100.EQ.2 ) THEN
            SL_CWX = VAR_CWX(5)
            TMSX = VAR_CWX(6)
            XGX = VAR_CWX(7)
            RHOGSX = 1.225D+0*GCPP(37,1)
            RHONSX = 1.D+3*GCPP(37,2)
            PX = P_CWX + PATM
            IF( 1.D+0-XGX.GT.EPSL ) THEN
              POSM_CWX = XGX*RHONSX/((1.D+0-XGX)*RHOGSX)
              CALL BO_BP( PBPX,POSM_CWX,T_CWX )
              IF( PX.GT.PBPX ) THEN
                SN_CWX = 1.D+0
                SG_CWX = 0.D+0
              ELSE
                PVWX = 0.D+0
                CALL BO_EQ_GAS( PX,PVWX,RHOGX,T_CWX,XGC_CWX,XGWX )
                CALL BO_EQ_OIL( RHONX,PX,POSM_CWX,T_CWX,XNC_CWX )
                SN_CWX = RHOGX*(XGX*(1.D+0-SL_CWX)-XGC_CWX(1)*
     &            (1.D+0-SL_CWX))/(RHONX*(XNC_CWX(1)-XGX) - 
     &            RHOGX*(XGC_CWX(1)-XGX))
                SG_CWX = MAX( 1.D+0-SL_CWX-SN_CWX,0.D+0 )
              ENDIF
            ELSE
              SG_CWX = 1.D+0
              SN_CWX = 0.D+0
            ENDIF
            ISR5X = MOD(ISRT(NS)/10,10)
            IZNX = 0
            CALL BO_IC2( POSM_CWX,P_CWX,RHOG_CWX,RHOL_CWX,RHON_CWX,
     &        SL_CWX,TMSX,T_CWX,XGC_CWX,XGW_CWX,XLS_CWX,XLW_CWX,XNC_CWX,
     &        YLSX,ISR5X,IZNX )
            RHOF_CWX = SG_CWX*RHOG_CWX + SL_CWX*RHOL_CWX + 
     &        SN_CWX*RHON_CWX
!
!---      Well state #3  ---
!
!         SL = 0.0
!         SG+SN = 1.0
!
!         Declared variables:
!
!         VAR_CWX(5) - water vapor relative saturation
!         VAR_CWX(7) - gas mass fraction of oil + gas
!        
          ELSEIF( ITS_CW(NTP,NCW)/100.EQ.3 ) THEN
            SL_CWX = 0.D+0
            PVWX = VAR_CWX(5)
            YLSX = 0.D+0
            XLS_CWX = 0.D+0
            XGX = VAR_CWX(7)
            RHOGSX = 1.225D+0*GCPP(37,1)
            RHONSX = 1.D+3*GCPP(37,2)
            PX = P_CWX + PATM
            IF( 1.D+0-XGX.GT.EPSL ) THEN
              POSM_CWX = XGX*RHONSX/((1.D+0-XGX)*RHOGSX)
              CALL BO_BP( PBPX,POSM_CWX,T_CWX )
              IF( PX.GT.PBPX ) THEN
                SN_CWX = 1.D+0
                SG_CWX = 0.D+0
              ELSE
                PVWX = 0.D+0
                CALL BO_EQ_GAS( PX,PVWX,RHOGX,T_CWX,XGC_CWX,XGWX )
                CALL BO_EQ_OIL( RHONX,PX,POSM_CWX,T_CWX,XNC_CWX )
                SN_CWX = RHOGX*(XGX*(1.D+0-SL_CWX)-XGC_CWX(1)*
     &            (1.D+0-SL_CWX))/(RHONX*(XNC_CWX(1)-XGX) - 
     &            RHOGX*(XGC_CWX(1)-XGX))
                SG_CWX = MAX( 1.D+0-SL_CWX-SN_CWX,0.D+0 )
              ENDIF
            ELSE
              SG_CWX = 1.D+0
              SN_CWX = 0.D+0
            ENDIF
            IZNX = 0
            CALL BO_IC3( POSM_CWX,PVWX,P_CWX,RHOG_CWX,RHON_CWX,
     &        T_CWX,XGC_CWX,XGW_CWX,XNC_CWX,YLSX,IZNX )
              RHOF_CWX = SG_CWX*RHOG_CWX + SN_CWX*RHON_CWX
          ENDIF
!
!---      Well gas viscosity, Pa s  ---
!
          IF( SG_CWX.GT.0.D+0 ) THEN
            CALL BO_VISG( RHOG_CWX,VISG_CWX )
          ENDIF
!
!---      Well aqueous viscosity, Pa s
!         concentration of CO2 in mol/kg,
!         concentration of salt in mol/kg  ---
!
          IF( SL_CWX.GT.0.D+0 ) THEN
            CMLAX = 0.D+0
            CMLSX = 1.D+3*XLS_CWX*RHOL_CWX/WTMS
            CALL VISC_L( CMLAX,CMLSX,P_CWX,T_CWX,VISL_CWX )
          ENDIF
!
!---      Well nonaqueous-liquid viscosity, Pa s  ---
!
          IF( SN_CWX.GT.0.D+0 ) THEN
            CALL BO_VISN( PX,POSM_CWX,T_CWX,VISN_CWX )
          ENDIF
!
!---      Equivalent field node radius components  ---
!
          PERMX = MAX( PERM(1,IZN),1.D-20 )
          PERMY = MAX( PERM(2,IZN),1.D-20 )
          PERMZ = MAX( PERM(3,IZN),1.D-20 )
          RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---      Cylindrical coordinates with azimuthal symmetry,
!         centrally located wells  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            ROZ = RP(I)
!
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROZ = PAR_CW(3,INVX)
            ENDIF
            RWX = MIN( RWX,9.999D-1*ROZ )
            PERMX = PERMRF(2,N)*PERM(1,IZN)
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
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROX = PAR_CW(3,INVX)
              ROY = PAR_CW(3,INVX)
              ROZ = PAR_CW(3,INVX)
            ENDIF
!
!---        Well index components  ---
!
            PERMX = PERMRF(2,N)*PERM(1,IZN)
            PERMY = PERMRF(2,N)*PERM(2,IZN)
            PERMZ = PERMRF(2,N)*PERM(3,IZN)
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
          DPNX = MAX( P_CWX-PGFX,0.D+0 )
          DPLX = MAX( P_CWX-PGFX,0.D+0 )
!
!---      Zero fluxes from well to reservoir  ---
!
          FXW_CWX(M,NWN) = 0.D+0
          DO 150 IGC = 1,NGC+2
            FXC_CWX(IGC,M,NWN) = 0.D+0
  150     CONTINUE
!
!---      Gas flux from well to reservoir  ---
!
          IF( SG_CWX.GT.0.D+0 ) THEN
            FGX_CWX = SG_CWX*WI_CWX*RHOG_CWX*DPGX/VISG_CWX
!
!---        Water component gas flux from well to reservoir  ---
!
            FXW_CWX(M,NWN) = FXW_CWX(M,NWN) + FGX_CWX*XGW_CWX
!
!---        Gas component gas flux from well to reservoir  ---
!
            FXC_CWX(1,M,NWN) = FXC_CWX(1,M,NWN) +  FGX_CWX*XGC_CWX(1)
          ENDIF
!
!---      Aqueous flux from well to reservoir  ---
!
          IF( SL_CWX.GT.0.D+0 ) THEN
            FLX_CWX = SL_CWX*WI_CWX*RHOL_CWX*DPLX/VISL_CWX
!
!---        Water component aqueous flux from well to reservoir  ---
!
            FXW_CWX(M,NWN) = FXW_CWX(M,NWN) + FLX_CWX*XLW_CWX
          ENDIF
!
!---      Nonaqueous-liquid flux from well to reservoir  ---
!
          IF( SN_CWX.GT.0.D+0 ) THEN
            FNX_CWX = SN_CWX*WI_CWX*RHON_CWX*DPNX/VISN_CWX
!
!---        Gas component Nonaqueous-liquid flux from well to 
!           reservoir  ---
!
            FXC_CWX(1,M,NWN) = FXC_CWX(1,M,NWN) +  FNX_CWX*XNC_CWX(1)
!
!---        Oil component Nonaqueous-liquid flux from well to 
!           reservoir  ---
!
            FXC_CWX(2,M,NWN) = FXC_CWX(2,M,NWN) +  FNX_CWX*XNC_CWX(2)
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
          DO 210 IGC = 1,NGC+2
            QM_CWX(M) = QM_CWX(M) + FXC_CWX(IGC,M,NWN)
  210     CONTINUE
  220   CONTINUE
  300 CONTINUE
!
!---  Consider current well pressure and well pressure limit  ---
!
      IF( ICHK_CWX.EQ.0 ) THEN
        ICHK_CWX = 1
!
!---    Well-limit pressure insufficient to produce specified
!       rate, well is pressure controlled  ---
!
        IF( QM_CWX(2).LT.VAR_CW2X(2) ) THEN
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
!---  End of INJP_BO_WELL group  ---
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
!     STOMP-EOR
!
!     Injection coupled well model
!     
!     Rate controlled or pressure controlled
!
!     Flux of energy, water mass, CO2 mass, CH4 mass, and petroleum
!     component masses from coupled-well nodes to field nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 April 2014.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE HYST
      USE GRID
      USE GLB_PAR
      USE FDVP
      USE FDVGC
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
      REAL*8 VAR_CWX(6+LNGC)
      REAL*8 ZMCX(LNGC)
      REAL*8 XGC_CWX(LNGC),XNC_CWX(LNGC),XMGCX(LNGC),XMNCX(LNGC)
      REAL*8 FXC_CWX(LNGC,3,LWN_CW),FXW_CWX(3,LWN_CW)
      REAL*8 QM_CWX(3),P_CWY(3),VAR_CW2X(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INJP_COUP_WELL'
!
!---  Black-oil option  ---
!
      IF( ISLC(68).EQ.0 ) THEN
        CALL INJP_BO_WELL( NCW )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  K-factors computed from fugacity coefficients  ---
!
      IFK = 0
!
!---  Loop over coupled-well nodes  ---
!
      DO 30 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---    Loop over increment indices  ---
!
        DO 20 M = 1,2
!
!---      Loop over petroleum components  ---
!
          DO 12 IGC = 1,NGC+2
            FXC_CWX(IGC,M,NWN) = 0.D+0
   12     CONTINUE
          FXW_CWX(M,NWN) = 0.D+0
   20   CONTINUE
   30  CONTINUE
!
!---  Injection well time interval ---
!
      TMZ = TM
      IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---  Cyclic time periods  ---
!
      IF( ICC_CW(NCW).EQ.1 ) THEN
!
!---    Loop over the injection well time periods, to find the 
!       final well time  ---
!
        NTX = 0
        DO 40 NTP = 1,IM_CW(NCW)
          NTX = NTX + IMP_CW(NTP,NCW)
   40   CONTINUE
!
!---    Determine time with the cyclic time period  ---
!
        TMZ = MOD( TM,VAR_CW(1,NTX,NCW) )
        IF( TM.GT.VAR_CW(1,NTX,NCW) ) THEN
          IF( TMZ.LT.EPSL ) TMZ = VAR_CW(1,NTX,NCW)
        ENDIF
      ENDIF
!
!---  Injection well is inactive set well pressure to be in 
!     equilibrium with reservoir  ---
!
      IF( TMZ.LE.VAR_CW(1,1,NCW) ) THEN
        CALL EQUIL_COUP_WELL( NCW )
        ID_CW(8,NCW) = 1
        GOTO 500
      ENDIF
!
!---  Loop over the injection well time periods  ---
!
      NS = 1
      DO 120 NTP = 1,IM_CW(NCW)
!
!---    Injection well time period only has one time (start time)  ---
!
        IF( IMP_CW(NTP,NCW).EQ.1 ) THEN
!
!---      Time prior to start time, injection well is inactive, set 
!         well pressure to be in equilibrium with reservoir  ---
!
          IF( TMZ.LE.VAR_CW(1,NS,NCW) ) THEN
            CALL EQUIL_COUP_WELL( NCW )
            ID_CW(8,NCW) = 1
            GOTO 500
!
!---      Time after start time, injection well is active  ---
!
          ELSE
            DO 80 N = 2,NGC+8
              VAR_CWX(N) = VAR_CW(N,1,NCW)
   80       CONTINUE
            GOTO 122
          ENDIF
!
!---    Injection well time period has multiple times  ---
!
        ELSE
          NE = NS + IMP_CW(NTP,NCW) - 1
!
!---      Time outside of injection well time period, go to next 
!         injection well time period  ---
!
          IF( TMZ.LE.VAR_CW(1,NS,NCW) .OR. 
     &      TMZ.GT.VAR_CW(1,NE,NCW) ) GOTO 110 
!
!---      Injection well time period has multiple time points, use  
!         linear interpolation of well parameters between time points  ---
!
          DO 100 M = 2,IMP_CW(NTP,NCW)
            MX = NS + M - 1
            IF( TMZ.LE.VAR_CW(1,MX,NCW) ) THEN
              TD_CW = VAR_CW(1,MX,NCW)-VAR_CW(1,MX-1,NCW)
              DT_CW = MIN( VAR_CW(1,MX,NCW)-TMZ,DT )
              TF_CW = (TMZ-VAR_CW(1,MX-1,NCW))/TD_CW
              DO 90 N = 2,NGC+8
                VAR_CWX(N) = VAR_CW(N,MX-1,NCW) + 
     &            TF_CW*(VAR_CW(N,MX,NCW)-VAR_CW(N,MX-1,NCW))
   90         CONTINUE
              GOTO 122
            ENDIF
  100     CONTINUE
        ENDIF
  110   CONTINUE
        NS = NS + IMP_CW(NTP,NCW)
  120 CONTINUE
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
        NPHAZX = 0
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          T_CWX = VAR_CWX(4)
        ELSE
          T_CWX = T(2,N)
        ENDIF
!
!---    Well state #1  ---
!
!       SL = 1.0
!       SG+SN = 0.0
!
!       Declared variables:
!
!       VAR_CWX(5) - aqueous CO2 relative saturation, or
!                    aqueous CO2 mass fraction
!       VAR_CWX(6) - aqueous salt relative saturation, or
!                         aqueous salt mass fraction
!
        IF( ITS_CW(NTP,NCW)/100.EQ.1 ) THEN
          PVAX = VAR_CWX(5)
          TMSX = VAR_CWX(6)
          ISR4X = MOD(ITS_CW(NTP,NCW),10)
          ISR5X = MOD(ITS_CW(NTP,NCW)/10,10)
          CALL FLH_IC1( P_CWX,PVAX,RHOL_CWX,RHOMLX,TMSX,T_CWX,XLA_CWX,
     &      XLS_CWX,XLW_CWX,XMLAX,XMLSX,XMLWX,YLSX,ZMCX,ISR4X,ISR5X )
          SG_CWX = 0.D+0
          SN_CWX = 0.D+0
          SL_CWX = 1.D+0
          RHOF_CWX = RHOL_CWX
!
!---    Well state #2  ---
!
!       SL < 1.0
!       SG+SN > 0.0
!
!       Declared variables:
!
!       VAR_CWX(5) - aqueous saturation
!       VAR_CWX(6) - aqueous salt relative saturation, or
!                    aqueous salt mass fraction
!       VAR_CWX(6+NGC) - total mole fraction of petroleum comp.
!
        ELSEIF( ITS_CW(NTP,NCW)/100.EQ.2 ) THEN
          SL_CWX = VAR_CWX(5)
          TMSX = VAR_CWX(6)
          ISR5X = MOD(ISRT(NS)/10,10)
          DO 130 IGC = 1,NGC+2
            JGC = 6 + IGC
            ZMCX(IGC) = VAR_CWX(JGC)
  130     CONTINUE
          IZNX = 0
          CALL FLH_IC2( P_CWX,RHOG_CWX,RHOMGX,RHOL_CWX,RHOMLX,
     &      RHON_CWX,RHOMNX,SG_CWX,SL_CWX,SN_CWX,TMSX,T_CWX,
     &      XGC_CWX,
     &      XMGCX,XGW_CWX,XMGWX,XLA_CWX,XLS_CWX,XLW_CWX,XMLAX,XMLSX,
     &      XMLWX,XNC_CWX,XMNCX,YLSX,ZGX,ZMCX,ZNX,ISR5X,IZNX,NPHAZX )
          RHOF_CWX = SG_CWX*RHOG_CWX + SL_CWX*RHOL_CWX + 
     &      SN_CWX*RHON_CWX
!
!---    Well state #3  ---
!
!       SL = 0.0
!       SG+SN = 1.0
!
!       Declared variables:
!
!       VAR_CWX(5) - water vapor relative saturation
!       VAR_CWX(6+) - total mole fraction of petroleum comp.
!      
        ELSEIF( ITS_CW(NTP,NCW)/100.EQ.3 ) THEN
          DO 132 IGC = 1,NGC+2
            JGC = 6 + IGC
            ZMCX(IGC) = VAR_CWX(JGC)
  132     CONTINUE
          SL_CWX = 0.D+0
          PVWX = VAR_CWX(5)
          YLSX = 0.D+0
          IZNX = 0
          CALL FLH_IC3( PVWX,P_CWX,RHOG_CWX,RHOMGX,RHON_CWX,RHOMNX,
     &      SG_CWX,SN_CWX,T_CWX,XGC_CWX,XGW_CWX,XMGCX,XMGWX,
     &      XNC_CWX,
     &      XMNCX,YLSX,ZGX,ZMCX,ZNX,IZNX,NPHAZX )
          RHOF_CWX = SG_CWX*RHOG_CWX + SN_CWX*RHON_CWX
        ENDIF
!
!---    Convert volumetric flow rate, m^3/s to mass flow rate, kg/s,
!       using fluid density at top of screened interval  ---
!
        IF( IT_CW(NCW).EQ.2 ) THEN
          VAR_CW2X(M) = VAR_CWX(2)*RHOF_CWX
        ELSE
          VAR_CW2X(M) = VAR_CWX(2)
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
          NPHAZX = 0
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            T_CWX = VAR_CWX(4)
          ELSE
            T_CWX = T(2,N)
          ENDIF
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
          P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOF_CWX
!
!---      Adjust the reservoir pressure to the coupled-well node
!         centroid  ---
!
!---      Phase condition #1  ---
!
!         SL = 1.0, SG+SN = 0.0
!
          IF( MOD(NPHAZ(2,N),10).EQ.1 ) THEN
            PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
            PNFX = PN(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(2,N)-MAX( PL(2,N),PN(2,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #2  ---
!
!         SL < 1.0, ST > 0.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.2 ) THEN
            PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
            PNFX = PN(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(2,N)-MAX( PL(2,N),PN(2,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #3  ---
!
!         SL < 1.0, ST > 0.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.3 ) THEN
            PGFX = PG(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(2,N)
            PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
            PNFX = PN(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
!
!---      Phase condition #4  ---
!
!         SL = 0.0, ST = 1.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.4 ) THEN
            PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
            PNFX = PN(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(2,N)-MAX( PL(2,N),PN(2,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #5  ---
!
!         SL = 0.0, ST = 1.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.5 ) THEN
            PGFX = PG(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(2,N)
            PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
            PNFX = PN(2,N) - (ZPX(2)-ZPNX)*GRAV*RHON(2,N)
          ENDIF
!
!---      Well state #1  ---
!
!         SL = 1.0
!         SG+SN = 0.0
!
!         Declared variables:
!
!         VAR_CWX(5) - aqueous CO2 relative saturation, or
!                      aqueous CO2 mass fraction
!         VAR_CWX(6) - aqueous salt relative saturation, or
!                           aqueous salt mass fraction
!
          IF( ITS_CW(NTP,NCW)/100.EQ.1 ) THEN
            PVAX = VAR_CWX(5)
            TMSX = VAR_CWX(6)
            ISR4X = MOD(ITS_CW(NTP,NCW),10)
            ISR5X = MOD(ITS_CW(NTP,NCW)/10,10)
            CALL FLH_IC1( P_CWX,PVAX,RHOL_CWX,RHOMLX,TMSX,T_CWX,XLA_CWX,
     &      XLS_CWX,XLW_CWX,XMLAX,XMLSX,XMLWX,YLSX,ZMCX,ISR4X,ISR5X )
            SG_CWX = 0.D+0
            SN_CWX = 0.D+0
            SL_CWX = 1.D+0
            RHOF_CWX = RHOL_CWX
!
!---      Well state #2  ---
!
!         SL < 1.0
!         SG+SN > 0.0
!
!         Declared variables:
!
!         VAR_CWX(5) - aqueous saturation
!         VAR_CWX(6) - aqueous salt relative saturation, or
!                      aqueous salt mass fraction
!         VAR_CWX(6+NGC) - total mole fraction of petroleum comp.
!
          ELSEIF( ITS_CW(NTP,NCW)/100.EQ.2 ) THEN
            SL_CWX = VAR_CWX(5)
            TMSX = VAR_CWX(6)
            ISR5X = MOD(ISRT(NS)/10,10)
            DO 140 IGC = 1,NGC+2
              JGC = 6 + IGC
              ZMCX(IGC) = VAR_CWX(JGC)
  140       CONTINUE
            IZNX = 0
            CALL FLH_IC2( P_CWX,RHOG_CWX,RHOMGX,RHOL_CWX,RHOMLX,
     &        RHON_CWX,RHOMNX,SG_CWX,SL_CWX,SN_CWX,TMSX,T_CWX,
     &        XGC_CWX,
     &        XMGCX,XGW_CWX,XMGWX,XLA_CWX,XLS_CWX,XLW_CWX,XMLAX,XMLSX,
     &        XMLWX,XNC_CWX,XMNCX,YLSX,ZGX,ZMCX,ZNX,ISR5X,IZNX,NPHAZX )
            RHOF_CWX = SG_CWX*RHOG_CWX + SL_CWX*RHOL_CWX + 
     &        SN_CWX*RHON_CWX
!
!---      Well state #3  ---
!
!         SL = 0.0
!         SG+SN = 1.0
!
!         Declared variables:
!
!         VAR_CWX(5) - water vapor relative saturation
!         VAR_CWX(6+) - total mole fraction of petroleum comp.
!        
          ELSEIF( ITS_CW(NTP,NCW)/100.EQ.3 ) THEN
            DO 142 IGC = 1,NGC+2
              JGC = 6 + IGC
              ZMCX(IGC) = VAR_CWX(JGC)
  142       CONTINUE
            SL_CWX = 0.D+0
            PVWX = VAR_CWX(5)
            YLSX = 0.D+0
            IZNX = 0
            CALL FLH_IC3( PVWX,P_CWX,RHOG_CWX,RHOMGX,RHON_CWX,RHOMNX,
     &        SG_CWX,SN_CWX,T_CWX,XGC_CWX,XGW_CWX,XMGCX,XMGWX,
     &        XNC_CWX,
     &        XMNCX,YLSX,ZGX,ZMCX,ZNX,IZNX,NPHAZX )
            RHOF_CWX = SG_CWX*RHOG_CWX + SN_CWX*RHON_CWX
          ENDIF
!
!---      Well gas viscosity, Pa s  ---
!
          IF( SG_CWX.GT.0.D+0 ) THEN
            CALL VISC_G( RHOMGX,T_CWX,VISG_CWX,XMGCX,XMGWX )
          ENDIF
!
!---      Well aqueous viscosity, Pa s
!         concentration of CO2 in mol/kg,
!         concentration of salt in mol/kg  ---
!
          IF( SL_CWX.GT.0.D+0 ) THEN
            CMLAX = 0.D+0
            CMLSX = 1.D+3*XLS_CWX*RHOL_CWX/WTMS
            CALL VISC_L( CMLAX,CMLSX,P_CWX,T_CWX,VISL_CWX )
          ENDIF
!
!---      Well nonaqueous-liquid viscosity, Pa s  ---
!
          IF( SN_CWX.GT.0.D+0 ) THEN
            CALL VISC_N( RHOMNX,T_CWX,VISN_CWX,XMNCX )
          ENDIF
!
!---      Equivalent field node radius components  ---
!
          PERMX = MAX( PERM(1,IZN),1.D-20 )
          PERMY = MAX( PERM(2,IZN),1.D-20 )
          PERMZ = MAX( PERM(3,IZN),1.D-20 )
          RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---      Cylindrical coordinates with azimuthal symmetry,
!         centrally located wells  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            ROZ = RP(I)
!
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROZ = PAR_CW(3,INVX)
            ENDIF
            RWX = MIN( RWX,9.999D-1*ROZ )
            PERMX = PERMRF(2,N)*PERM(1,IZN)
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
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROX = PAR_CW(3,INVX)
              ROY = PAR_CW(3,INVX)
              ROZ = PAR_CW(3,INVX)
            ENDIF
!
!---        Well index components  ---
!
            PERMX = PERMRF(2,N)*PERM(1,IZN)
            PERMY = PERMRF(2,N)*PERM(2,IZN)
            PERMZ = PERMRF(2,N)*PERM(3,IZN)
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
          DPNX = MAX( P_CWX-PGFX,0.D+0 )
          DPLX = MAX( P_CWX-PGFX,0.D+0 )
!
!---      Zero fluxes from well to reservoir  ---
!
          FXW_CWX(M,NWN) = 0.D+0
          DO 150 IGC = 1,NGC+2
            FXC_CWX(IGC,M,NWN) = 0.D+0
  150     CONTINUE
!
!---      Gas flux from well to reservoir  ---
!
          IF( SG_CWX.GT.0.D+0 ) THEN
            FGX_CWX = SG_CWX*WI_CWX*RHOG_CWX*DPGX/VISG_CWX
            FXW_CWX(M,NWN) = FXW_CWX(M,NWN) + FGX_CWX*XGW_CWX
            DO 152 IGC = 1,NGC+2
              FXC_CWX(IGC,M,NWN) = FXC_CWX(IGC,M,NWN) + 
     &          FGX_CWX*XGC_CWX(IGC)
  152       CONTINUE
          ENDIF
!
!---      Aqueous flux from well to reservoir  ---
!
          IF( SL_CWX.GT.0.D+0 ) THEN
            FLX_CWX = SL_CWX*WI_CWX*RHOL_CWX*DPLX/VISL_CWX
            FXW_CWX(M,NWN) = FXW_CWX(M,NWN) + FLX_CWX*XLW_CWX
            FXC_CWX(1,M,NWN) = FXC_CWX(1,M,NWN) + FLX_CWX*XLA_CWX
          ENDIF
!
!---      Nonaqueous-liquid flux from well to reservoir  ---
!
          IF( SN_CWX.GT.0.D+0 ) THEN
            FNX_CWX = SN_CWX*WI_CWX*RHON_CWX*DPNX/VISN_CWX
            DO 154 IGC = 1,NGC+2
              FXC_CWX(IGC,M,NWN) = FXC_CWX(IGC,M,NWN) + 
     &          FNX_CWX*XNC_CWX(IGC)
  154       CONTINUE
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
          DO 210 IGC = 1,NGC+2
            QM_CWX(M) = QM_CWX(M) + FXC_CWX(IGC,M,NWN)
  210     CONTINUE
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
        IF( ID_CW(8,NCW).EQ.1 .AND. NITER.GE.3 ) THEN
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
        ELSEIF( QM_CWX(1).GT.VAR_CW2X(1) ) THEN
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
!     STOMP-EOR
!
!     Modify Jacobian matrix for the coupled-well equations
!     and load Jacobian matrix for the coupled-well equations
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 20 April 2011.
!
#ifdef petsc
!
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
#endif
#ifdef lis
!
!----------------------PETSc Modules-----------------------------------!
!
      USE STOMP_LIS_MODULE
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FDVP
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#endif
#ifdef lis
#include "include/lisf.h"
#endif
!
!----------------------Type Declarations-------------------------------!
!
#ifdef petsc
      REAL*8 BUFFER(LUK+1)
      INTEGER IROW(1), ICOL(LUK+1)
#endif
!#ifdef lis
!      LIS_SCALAR VALUEX
!      LIS_INTEGER IROWX,ICOLX,IERR
!#endif
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
        DO 1500 NCW = 1,N_CW
!
!---      Loop over coupled-well well nodes ---
!
          DO 1200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
            N = IWN_CW(NWN)
            NMD = ABS(IXP(N))
!
!---        Energy equation, unless isothermal simulation  ---
!
            IF( ISLC(30).EQ.0 ) THEN
!
!---          Energy balance equation at field node ---
!
              MP = IM(IEQT,NMD)
!
!---          Change in energy flow into field node with respect
!             to change in field node primary variables  ---
!
              DO 1090 M = 1,ISVC
                MCOL = IM(M,NMD)
                MROW = MP-MCOL+MDC
                ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &            (FXE_CW(1,NWN)-FXE_CW(M+1,NWN))/DNR(M,N)
 1090         CONTINUE
!
!---          Change in energy flow into field node with respect
!             to change in coupled-well pressure  ---
!
              MCOL = JM_CW(NCW)
              MROW = MP-MCOL+MDC
              MX = ISVC+2
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &          (FXE_CW(1,NWN)-FXE_CW(MX,NWN))/DNR_CW(NCW)
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) ALU(MROW,MCOL) = 0.D+0
!
!---          Energy flow into field node, W  ---
!
              BLU(MP) = BLU(MP) + FXE_CW(1,NWN)
              RSDL(IEQT,N) = BLU(MP)
            ENDIF
!
!---        Water mass balance equation at field node ---
!
            MP = IM(IEQW,NMD)
!
!---        Change in water mass flux into field node with respect
!           to change in field node primary variables  ---
!
            DO 1100 M = 1,ISVC
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &          (FXW_CW(1,NWN)-FXW_CW(M+1,NWN))/DNR(M,N)
 1100       CONTINUE
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
!---        Pressure controlled coupled well  ---
!
            IF( ID_CW(8,NCW).EQ.1 ) ALU(MROW,MCOL) = 0.D+0
!
!---        Water mass flux into field node, kg/s  ---
!
            BLU(MP) = BLU(MP) + FXW_CW(1,NWN)
            RSDL(IEQW,N) = BLU(MP)
!
!---        Petroleum component mass balance equation at field node ---
!
            DO 1120 IGC = 1,NGC+2
!
!---          Iso-CO2 option  ---
!
              IF( IGC.EQ.1 .AND. ISLC(45).NE.0 ) CYCLE
!
!---          Iso-CH4 option  ---
!
              IF( IGC.EQ.2 .AND. ISLC(64).NE.0 ) CYCLE
              MP = IM(IEQGC(IGC),NMD)
!
!---          Change in petroleum component mass flux into field node 
!             with respect to change in field node primary variables  ---
!
              DO 1110 M = 1,ISVC
                MCOL = IM(M,NMD)
                MROW = MP-MCOL+MDC
                ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &            (FXC_CW(IGC,1,NWN)-FXC_CW(IGC,M+1,NWN))/DNR(M,N)
 1110         CONTINUE
!
!---          Change in petroleum component mass flux into field node
!             with respect to change in coupled-well pressure  ---
!
              MCOL = JM_CW(NCW)
              MROW = MP-MCOL+MDC
              MX = ISVC+2
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &          (FXC_CW(IGC,1,NWN)-FXC_CW(IGC,MX,NWN))/DNR_CW(NCW)
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) ALU(MROW,MCOL) = 0.D+0
!
!---          Petroleum component mass flux into field node, kg/s  ---
!
              BLU(MP) = BLU(MP) + FXC_CW(IGC,1,NWN)
              RSDL(IEQGC(IGC),N) = BLU(MP)
 1120       CONTINUE
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
          DO 1400 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
            N = IWF_CW(NWF)
            NMD = ABS(IXP(N))
            MX = (NWF-ID_CW(5,NCW))*ISVC + 1
!
!---        Change in coupled-well mass balance with respect to
!           change in field node primary variables  ---
!
            DO 1300 M = 1,ISVC
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &          (RS_CW(MX+M,NCW)-RS_CW(1,NCW))/DNR(M,N)
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) ALU(MROW,MCOL) = 0.D+0
 1300       CONTINUE
 1400     CONTINUE
 1500   CONTINUE
!
!---  SPLib or Lis solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---    Loop over coupled wells ---
!
        DO 3500 NCW = 1,N_CW
!
!---      Loop over coupled-well well nodes ---
!
          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
          DO 3200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
            N = IWN_CW(NWN)
            NMD = ABS(IXP(N))
!
!---        Energy equation, unless isothermal simulation  ---
!
            IF( ISLC(30).EQ.0 ) THEN
!
!---          Energy balance equation at field node ---
!
              MP = IM(IEQT,NMD)
!
!---          Change in energy flow into field node with respect
!             to change in field node primary variables  ---
!
              MA = 0
!
!---          Block refinement scheme  ---
!
              IF( ISLC(33).EQ.1 ) MA = 0
              DO 3090 M = 1,ISVC
                MCOL = KLU(MP,M+MA)
                DLU(MCOL) = DLU(MCOL) + 
     &            (FXE_CW(1,NWN)-FXE_CW(M+1,NWN))/DNR(M,N)
 3090         CONTINUE
!
!---          Change in energy flow with respect
!             to change in coupled-well pressure  ---
!
              MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + 
     &          IEQT + 1
              MCOL = KLU_CW(MA,NCW)
              MX = ISVC+2
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXE_CW(1,NWN)-FXE_CW(MX,NWN))/DNR_CW(NCW)
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 0.D+0
!
!---          Energy flow into field node, W  ---
!
              BLU(MP) = BLU(MP) + FXE_CW(1,NWN)
              RSDL(IEQT,N) = BLU(MP)
            ENDIF
!
!---        Water mass balance equation at field node ---
!
            MP = IM(IEQW,NMD)
!
!---        Change in water mass flux into field node with respect
!           to change in field node primary variables  ---
!
            MA = 0
!
!---        Block refinement scheme  ---
!
            IF( ISLC(33).EQ.1 ) MA = 0
            DO 3100 M = 1,ISVC
              MCOL = KLU(MP,M+MA)
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXW_CW(1,NWN)-FXW_CW(M+1,NWN))/DNR(M,N)
 3100       CONTINUE
!
!---        Change in water fraction of coupled-well flux with respect
!           to change in coupled-well pressure  ---
!
            MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + IEQW + 1
            MCOL = KLU_CW(MA,NCW)
            MX = ISVC+2
            DLU(MCOL) = DLU(MCOL) + 
     &        (FXW_CW(1,NWN)-FXW_CW(MX,NWN))/DNR_CW(NCW)
!
!---        Pressure controlled coupled well  ---
!
            IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 0.D+0
!
!---        Water mass flux into field node, kg/s  ---
!
            BLU(MP) = BLU(MP) + FXW_CW(1,NWN)
            RSDL(IEQW,N) = BLU(MP)
!
!---        Petroleum component mass balance equation at field node ---
!
            DO 3120 IGC = 1,NGC+2
!
!---          Iso-CO2 option  ---
!
              IF( IGC.EQ.1 .AND. ISLC(45).NE.0 ) CYCLE
!
!---          Iso-CH4 option  ---
!
              IF( IGC.EQ.2 .AND. ISLC(64).NE.0 ) CYCLE
              MP = IM(IEQGC(IGC),NMD)
!
!---          Change in CO2 fraction of coupled-well flux with respect
!             to change in field node primary variables  ---
!
              MA = 0
!
!---          Block refinement scheme  ---
!
              IF( ISLC(33).EQ.1 ) MA = 0
              DO 3110 M = 1,ISVC
                MCOL = KLU(MP,M+MA)
                DLU(MCOL) = DLU(MCOL) + 
     &            (FXC_CW(IGC,1,NWN)-FXC_CW(IGC,M+1,NWN))/DNR(M,N)
 3110         CONTINUE
!
!---          Change in petroleum component fraction of coupled-well 
!             flux with respect to change in coupled-well pressure  ---
!
              MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + 
     &          IEQGC(IGC) + 1
              MCOL = KLU_CW(MA,NCW)
              MX = ISVC+2
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXC_CW(IGC,1,NWN)-FXC_CW(IGC,MX,NWN))/DNR_CW(NCW)
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 0.D+0
!
!---          Petroleum component mass flux into field node, kg/s  ---
!
              BLU(MP) = BLU(MP) + FXC_CW(IGC,1,NWN)
              RSDL(IEQGC(IGC),N) = BLU(MP)
 3120       CONTINUE
 3200     CONTINUE
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
          DO 3400 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
            N = IWF_CW(NWF)
            NMD = ABS(IXP(N))
            MX = (NWF-ID_CW(5,NCW))*ISVC + 1
            MA = (NWF-ID_CW(5,NCW))*ISVC + 1
!
!---        Change in coupled-well mass balance with respect to
!           change in field node primary variables  ---
!
            DO 3300 M = 1,ISVC
              MCOL = KLU_CW(M+MA,NCW)
              DLU(MCOL) = DLU(MCOL) + 
     &          (RS_CW(MX+M,NCW)-RS_CW(1,NCW))/DNR(M,N)
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 0.D+0
 3300       CONTINUE
 3400     CONTINUE
 3500   CONTINUE
#ifdef petsc
!
!---  PETSc solver  ---
!
      ELSEIF( ILES.EQ.5 ) THEN
!
!---    Loop over coupled wells ---
!
        DO 5500 NCW = 1,N_CW
!
!---      Loop over coupled-well well nodes ---
!
          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
          DO 5200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
            N = IWN_CW(NWN)
            NMD = ABS(IXP(N))
!
!---        Energy equation, unless nonisothermal simulation  ---
!
            IF( ISLC(30).EQ.0 ) THEN
!
!---          Energy balance equation at field node ---
!
              MP = IM(IEQT,NMD)
              IROW(1) = MP-1
!
!---          Change in energy flow into field node with respect
!             to change in field node primary variables  ---
!
              DO 5090 M = 1,ISVC
                 ICOL(M) = IM(M,NMD)-1
                 BUFFER(M) = (FXE_CW(1,NWN)-FXE_CW(M+1,NWN))/DNR(M,N)
 5090         CONTINUE
!
!---          Change in energy flow with respect
!             to change in coupled-well pressure  ---
!
              MX = ISVC+2
              M = ISVC+1
              ICOL(M) = JM_CW(NCW)-1
              BUFFER(M) = (FXE_CW(1,NWN)-FXE_CW(MX,NWN))/DNR_CW(NCW)
              CALL MatSetValues(F_MAT,1,IROW,ISVC+1,ICOL,BUFFER,
     &                          ADD_VALUES,IERR)
!
!---          Energy flow into field node, W  ---
!
              BLU(MP) = BLU(MP) + FXE_CW(1,NWN)
              RSDL(IEQT,N) = BLU(MP)
            ENDIF
!
!---        Water mass balance equation at field node ---
!
            MP = IM(IEQW,NMD)
            IROW(1) = MP-1
!
!---        Change in water mass flow into field node with respect
!           to change in field node primary variables  ---
!
            DO 5100 M = 1,ISVC
               ICOL(M) = IM(M,NMD)-1
               BUFFER(M) = (FXW_CW(1,NWN)-FXW_CW(M+1,NWN))/DNR(M,N)
 5100       CONTINUE
!
!---        Change in water mass flow with respect
!           to change in coupled-well pressure  ---
!
            MX = ISVC+2
            M = ISVC+1
            ICOL(M) = JM_CW(NCW)-1
            BUFFER(M) = (FXW_CW(1,NWN)-FXW_CW(MX,NWN))/DNR_CW(NCW)
            CALL MatSetValues(F_MAT,1,IROW,ISVC+1,ICOL,BUFFER,
     &                        ADD_VALUES,IERR)
!
!---        Water mass flow into field node, kg/s  ---
!
            BLU(MP) = BLU(MP) + FXW_CW(1,NWN)
            RSDL(IEQW,N) = BLU(MP)
 5102       CONTINUE
!
!---        Petroleum-component mass balance equation at field node ---
!
            DO 5120 IGC = 1,NGC+2
!
!---          Iso-CO2 option  ---
!
              IF( IGC.EQ.1 .AND. ISLC(45).NE.0 ) CYCLE
!
!---          Iso-CH4 option  ---
!
              IF( IGC.EQ.2 .AND. ISLC(64).NE.0 ) CYCLE
              MP = IM(IEQGC(IGC),NMD)
              IROW(1) = MP-1
!
!---          Change in CO2 mass flow with respect
!             to change in field node primary variables  ---
!
              DO 5110 M = 1,ISVC
                ICOL(M) = IM(M,NMD)-1
                BUFFER(M) = (FXC_CW(IGC,1,NWN)-FXC_CW(IGC,M+1,NWN))/
     &            DNR(M,N)
 5110         CONTINUE
!
!---          Change in CO2 mass flow with respect
!             to change in coupled-well pressure  ---
!
              MX = ISVC+2
              M = ISVC+1
              ICOL(M) = JM_CW(NCW)-1
              BUFFER(M) = (FXC_CW(IGC,1,NWN)-FXC_CW(IGC,MX,NWN))/
     &          DNR_CW(NCW)
              CALL MatSetValues(F_MAT,1,IROW,ISVC+1,ICOL,BUFFER,
     &                        ADD_VALUES,IERR)
!
!---          CO2 mass flow into field node, kg/s  ---
!
              BLU(MP) = BLU(MP) + FXC_CW(IGC,1,NWN)
              RSDL(IEQGC(IGC),N) = BLU(MP)
 5120       CONTINUE
 5200     CONTINUE
!
!---      Coupled-well mass balance  ---
!
          MP = JM_CW(NCW)
          IROW(1) = MP-1
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
          M = 1
          ICOL(M) = JM_CW(NCW)-1
          BUFFER(M) = (RS_CW(MX,NCW)-RS_CW(1,NCW))/DNR_CW(NCW)
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) BUFFER(M) = 1.D+0
          CALL MatSetValues(F_MAT,1,IROW,1,ICOL,BUFFER,
     &                      ADD_VALUES,IERR)
!
!---      Loop over field nodes with coupled-well nodes ---
!
          DO 5400 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
            N = IWF_CW(NWF)
            NMD = ABS(IXP(N))
            MX = (NWF-ID_CW(5,NCW))*ISVC + 1
!
!---        Change in coupled-well mass balance with respect to
!           change in field node primary variables  ---
!
            DO 5300 M = 1,ISVC
              ICOL(M) = IM(M,NMD)-1
              BUFFER(M) = (RS_CW(MX+M,NCW)-RS_CW(1,NCW))/DNR(M,N)
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) BUFFER(M) = 0.D+0
 5300       CONTINUE
            CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                        ADD_VALUES,IERR)
 5400     CONTINUE
 5011   CONTINUE
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
!     STOMP-EOR
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
      SUBROUTINE PROD_BO_WELL( NCW )
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
!     Production black-oil well model
!     
!     Bottom-hole pressure controlled
!
!     Flux of energy, water mass, CO2 mass, CH4 mass, and petroleum
!     component masses from field nodes to coupled-well nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 April 2014.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE JACOB
      USE HYST
      USE GRID
      USE GLB_PAR
      USE FDVP
      USE FDVN
      USE FDVGC
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
      REAL*8 VAR_CWX(6+LNGC)
      INTEGER, SAVE :: IALLOC
      DATA IALLOC /0/
      INTEGER, SAVE, DIMENSION(:), ALLOCATABLE :: MCW,MFD
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROD_BO_WELL'
      IF( IALLOC.EQ.0 ) THEN
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
   10   CONTINUE        
        IALLOC = 1
      ENDIF
!
!---  Zero production well fluxes ---
!
      QM_CW(1,NCW) = 0.D+0
      QM_CW(3,NCW) = 0.D+0
      QM_CW(5,NCW) = 0.D+0
      DO 11 IGC = 1,NGC+2
        QM_CW((7+(IGC-1)*2),NCW) = 0.D+0
   11 CONTINUE
!
!---  Loop over production well nodes from bottom to top  ---
!
      DO 30 NWN = ID_CW(4,NCW),ID_CW(3,NCW),-1
!
!---    Zero volumetric production well fluxes
!
!       Q_CW(1,NWN) - total volumetric flux, m^3/s
!       Q_CW(2,NWN) - aqueous volumetric flux, m^3/s
!       Q_CW(3,NWN) - gas volumetric flux, m^3/s
!       Q_CW(4,NWN) - nonaqueous-liquid volumetric flux, m^3/s  ---
!
        DO M = 1,4
          Q_CW(M,NWN) = 0.D+0
        ENDDO
!
!---    Loop over increment indices  ---
!
        DO 20 M = 1,ISVC+2
          DO 12 IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = 0.D+0
   12     CONTINUE
          FXW_CW(M,NWN) = 0.D+0
   20   CONTINUE
   30  CONTINUE
!
!---  Production well time interval ---
!
      TMZ = TM
      IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---  Cyclic time periods  ---
!
      IF( ICC_CW(NCW).EQ.1 ) TMZ = MOD( TM,VAR_CW(1,IM_CW(NCW),NCW) )
!
!---  Time prior to start time, production well is inactive, set 
!     well pressure to be in equilibrium with reservoir  ---
!
      IF( TMZ.LE.VAR_CW(1,1,NCW) ) THEN
        CALL EQUIL_COUP_WELL( NCW )
        GOTO 500
      ENDIF
!
!---  Time after start time, production well is active, with single
!     start time  ---
!
      IF( IM_CW(NCW).EQ.1 ) THEN
        DO 80 N = 2,7
          VAR_CWX(N) = VAR_CW(N,1,NCW)
   80   CONTINUE
         GOTO 122
!
!---  Time after start time, production well is active, with multiple
!     time points, use linear interpolation of well parameters 
!     between time points  ---
!
      ELSE
        DO 100 M = 2,IM_CW(NCW)
          IF( TMZ.LE.VAR_CW(1,M,NCW) ) THEN
            TD_CW = VAR_CW(1,M,NCW)-VAR_CW(1,M-1,NCW)
            DT_CW = MIN( VAR_CW(1,M,NCW)-TMZ,DT )
            TF_CW = (TMZ-VAR_CW(1,M-1,NCW))/TD_CW
            DO 90 N = 2,7
              VAR_CWX(N) = VAR_CW(N,M-1,NCW) + 
     &          TF_CW*(VAR_CW(N,M,NCW)-VAR_CW(N,M-1,NCW))
   90       CONTINUE
            GOTO 122
          ENDIF
  100   CONTINUE
      ENDIF
!
!---  Production well is inactive set well pressure to be in 
!     equilibrium with reservoir  ---
!
      CALL EQUIL_COUP_WELL( NCW )
      GOTO 500
  122 CONTINUE
!
!---  Set bottom-hole pressure  ---
!
      PL_CW(NCW) = VAR_CWX(3) - PATM
      P_CW(2,NCW) = PL_CW(NCW)
      P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
!
!---  Pressure controlled coupled well  ---
!
      ID_CW(8,NCW) = 1
!
!---  Well fluid density based on relative fluid mobility
!     at the previous time step  ---
!
      SRHOFX = 0.D+0
      SVOLX = 0.D+0
!
!---  Store bottom of coupled-well location in previous
!     coupled-well node location  ---
!
      XPX(1) = XTP_CW(2,ID_CW(2,NCW))
      YPX(1) = YTP_CW(2,ID_CW(2,NCW))
      ZPX(1) = ZTP_CW(2,ID_CW(2,NCW))
!
!---  Loop over the nodes in the coupled well from bottom to top ---
!
      DO 130 NWN = ID_CW(4,NCW),ID_CW(3,NCW),-1
        N = IWN_CW(NWN)
        I = ID(N)
        INVX = INV_CW(NWN)
        IZN = IZ(N)
!
!---    Coupled-well node centroids and projections ---
!
        XLX = PLX_CW(NWN)
        YLX = PLY_CW(NWN)
        ZLX = PLZ_CW(NWN)
        XPX(2) = 5.D-1*(XP_CW(2,NWN)+XP_CW(1,NWN))
        YPX(2) = 5.D-1*(YP_CW(2,NWN)+YP_CW(1,NWN))
        ZPX(2) = 5.D-1*(ZP_CW(2,NWN)+ZP_CW(1,NWN))
!
!---    Cylindrical coordinates with azimuthal symmetry,
!       centrally located wells  ---
!
        IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
          XPNX = 0.D+0
          YPNX = 0.D+0
          ZPNX = ZP(N)
!
!---    Cylindrical coordinates  ---
!
        ELSEIF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
          XPNX = XP(N)*COS(YP(N))
          YPNX = XP(N)*SIN(YP(N))
          ZPNX = ZP(N)
!
!---    Cartesian or boundary-fitted orthogonal coordinates  ---
!
        ELSE
          XPNX = XP(N)
          YPNX = YP(N)
          ZPNX = ZP(N)
        ENDIF
!
!---    Equivalent field node radius components  ---
!
        PERMX = MAX( PERM(1,IZN),1.D-20 )
        PERMY = MAX( PERM(2,IZN),1.D-20 )
        PERMZ = MAX( PERM(3,IZN),1.D-20 )
        RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---    Cylindrical coordinates with azimuthal symmetry,
!       centrally located wells  ---
!
        IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
          ROZ = RP(I)
!
!---      Specified pressure-equivalent radius  ---
!
          IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
            ROZ = PAR_CW(3,INVX)
          ENDIF
          RWX = MIN( RWX,9.999D-1*ROZ )
          PERMX = PERMRF(1,N)*PERM(1,IZN)
          PERMY = PERMRF(1,N)*PERM(2,IZN)
          WI_CWX = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &      (LOG(ROZ/RWX)+PAR_CW(1,INVX))
        ELSE
          PERMYZ = SQRT(PERMY/PERMZ)
          PERMZY = SQRT(PERMZ/PERMY)
          DXGFX = DXGF(N)/FF_CW(1,NCW)
          DYGFX = DYGF(N)*RP(I)/FF_CW(2,NCW)
          DZGFX = DZGF(N)/FF_CW(3,NCW)
          ROX = 2.8D-1*SQRT(PERMYZ*(DZGFX**2) + PERMZY*(DYGFX**2))
     &    /(SQRT(PERMYZ)+SQRT(PERMZY))
          PERMZX = SQRT(PERMZ/PERMX)
          PERMXZ = SQRT(PERMX/PERMZ)
          ROY = 2.8D-1*SQRT(PERMZX*(DXGFX**2) + PERMXZ*(DZGFX**2))
     &      /(SQRT(PERMZX)+SQRT(PERMXZ))
          PERMYX = SQRT(PERMY/PERMX)
          PERMXY = SQRT(PERMX/PERMY)
          ROZ = 2.8D-1*SQRT(PERMYX*(DXGFX**2) + PERMXY*(DYGFX**2))
     &      /(SQRT(PERMYX)+SQRT(PERMXY))
!
!---      Specified pressure-equivalent radius  ---
!
          IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
            ROX = PAR_CW(3,INVX)
            ROY = PAR_CW(3,INVX)
            ROZ = PAR_CW(3,INVX)
          ENDIF
!
!---      Well index components  ---
!
          PERMX = PERMRF(1,N)*PERM(1,IZN)
          PERMY = PERMRF(1,N)*PERM(2,IZN)
          PERMZ = PERMRF(1,N)*PERM(3,IZN)
          WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/
     &      (LOG(ROX/RWX)+PAR_CW(1,INVX))
          WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/
     &      (LOG(ROY/RWX)+PAR_CW(1,INVX))
          WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &      (LOG(ROZ/RWX)+PAR_CW(1,INVX))
          WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
        ENDIF
!
!---    Mass fluxes, positive into the well  ---
!
        DPGX = 1.D+0
        DPNX = 1.D+0
        DPLX = 1.D+0
!
!---    Liquid production well, no gas flux from reservoir to
!       well  ---
!
        IF( IT_CW(NCW).EQ.-2 ) THEN
          FGX_CWX = 0.D+0
        ELSE
!
!---      Gas flux from reservoir to well  ---
!
          FGX_CWX = WI_CWX*RKG(1,N)*RHOG(1,N)*DPGX/VISG(1,N)
        ENDIF
!
!---    Aqueous flux from well to reservoir  ---
!
        FLX_CWX = WI_CWX*RKL(1,1,N)*RHOL(1,N)*DPLX/VISL(1,N)
!
!---    Nonaqueous-liquid flux from well to reservoir  ---
!
        FNX_CWX = WI_CWX*RKN(1,N)*RHON(1,N)*DPNX/VISN(1,N)
!
!---    Total mass and volumetric flux into well from reservoir  ---
!
        SRHOFX = SRHOFX + FGX_CWX + FLX_CWX + FNX_CWX
        SVOLX = SVOLX + FGX_CWX/RHOG(1,N) + FLX_CWX/RHOL(1,N) + 
     &    FNX_CWX/RHON(1,N)
!
!---    Store current coupled-well node location in previous
!       coupled-well node location  ---
!
        XPX(1) = XPX(2)
        YPX(1) = YPX(2)
        ZPX(1) = ZPX(2)
  130 CONTINUE
!
!---  Dynamic average well-fluid density  ---
!
      IF( SRHOFX.GT.EPSL ) THEN
        RHOF_CWX = SRHOFX/SVOLX
!
!---  Static average well-fluid density  ---
!
      ELSE
        SRHOFX = 0.D+0
        SVOLX = 0.D+0
!
!---  Loop over the nodes in the coupled well from bottom to top ---
!
        DO 140 NWN = ID_CW(4,NCW),ID_CW(3,NCW),-1
          N = IWN_CW(NWN)
!
!---      Liquid production well, no gas flux from reservoir to
!         well  ---
!
          IF( IT_CW(NCW).EQ.-2 ) THEN
            SRHOFX = SRHOFX + VOL(N)*PORD(1,N)* 
     &        (RHOL(1,N)/VISL(1,N) + RHON(1,N)/VISN(1,N))
            SVOLX = SVOLX + VOL(N)*PORD(1,N)* 
     &        (1.D+0/VISL(1,N) +  1.D+0/VISN(1,N))
          ELSE
            SRHOFX = SRHOFX + VOL(N)*PORD(1,N)*
     &       (RHOG(1,N)/VISG(1,N) + RHOL(1,N)/VISL(1,N) + 
     &        RHON(1,N)/VISN(1,N))
            SVOLX = SVOLX + VOL(N)*PORD(1,N)*
     &       (1.D+0/VISG(1,N) + 1.D+0/VISL(1,N) + 1.D+0/VISN(1,N))
          ENDIF
  140   CONTINUE
        IF( SRHOFX.GT.EPSL ) THEN
          RHOF_CWX = SRHOFX/SVOLX
!
!---    Reference well-fluid density  ---
!
        ELSE
          RHOF_CWX = RHORL
        ENDIF
      ENDIF
!
!---  Record well fluid density for output  ---
!
      RHOF_CW(NCW) = RHOF_CWX
!
!---  Loop over increment indices ---
!
      DO 300 M = 1,ISVC+2
        MW = MCW(M)
        MF = MFD(M)
        N = IWN_CW(ID_CW(3,NCW))
        P_CWX = P_CW(MW,NCW)
!
!---    Store bottom of coupled-well location in previous
!       coupled-well node location  ---
!
        XPX(1) = XTP_CW(2,ID_CW(2,NCW))
        YPX(1) = YTP_CW(2,ID_CW(2,NCW))
        ZPX(1) = ZTP_CW(2,ID_CW(2,NCW))
!
!---    Loop over the nodes in the coupled well from bottom to top ---
!
        DO 200 NWN = ID_CW(4,NCW),ID_CW(3,NCW),-1
          N = IWN_CW(NWN)
          I = ID(N)
          INVX = INV_CW(NWN)
          IZN = IZ(N)
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
          P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOF_CWX
!
!---      Adjust the reservoir pressure to the coupled-well node
!         centroid  ---
!
!---      Phase condition #1  ---
!
!         SL = 1.0, SG+SN = 0.0
!
          IF( MOD(NPHAZ(2,N),10).EQ.1 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #2  ---
!
!         SL < 1.0, ST > 0.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.2 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #3  ---
!
!         SL < 1.0, ST > 0.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.3 ) THEN
            PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(MF,N)
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
!
!---      Phase condition #4  ---
!
!         SL = 0.0, ST = 1.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.4 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #5  ---
!
!         SL = 0.0, ST = 1.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.5 ) THEN
            PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(MF,N)
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
          ENDIF
!
!---      Equivalent field node radius components  ---
!
          PERMX = MAX( PERM(1,IZN),1.D-20 )
          PERMY = MAX( PERM(2,IZN),1.D-20 )
          PERMZ = MAX( PERM(3,IZN),1.D-20 )
          RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---      Cylindrical coordinates with azimuthal symmetry,
!         centrally located wells  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            ROZ = RP(I)
!
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROZ = PAR_CW(3,INVX)
            ENDIF
            RWX = MIN( RWX,9.999D-1*ROZ )
            PERMX = PERMRF(MF,N)*PERM(1,IZN)
            PERMY = PERMRF(MF,N)*PERM(2,IZN)
            WI_CWX = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &        (LOG(ROZ/RWX)+PAR_CW(1,INVX))
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
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROX = PAR_CW(3,INVX)
              ROY = PAR_CW(3,INVX)
              ROZ = PAR_CW(3,INVX)
            ENDIF
!
!---        Well index components  ---
!
            PERMX = PERMRF(MF,N)*PERM(1,IZN)
            PERMY = PERMRF(MF,N)*PERM(2,IZN)
            PERMZ = PERMRF(MF,N)*PERM(3,IZN)
            WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/
     &        (LOG(ROX/RWX)+PAR_CW(1,INVX))
            WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/
     &        (LOG(ROY/RWX)+PAR_CW(1,INVX))
            WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &        (LOG(ROZ/RWX)+PAR_CW(1,INVX))
            WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
          ENDIF
!
!---      Mass fluxes, positive into the well  ---
!
          DPGX = MAX( PGFX-P_CWX,0.D+0 )
          DPNX = MAX( PGFX-P_CWX,0.D+0 )
          DPLX = MAX( PGFX-P_CWX,0.D+0 )
!
!---      Zero fluxes from well to reservoir  ---
!
          FXW_CW(M,NWN) = 0.D+0
          DO 150 IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = 0.D+0
  150     CONTINUE
!
!---      Liquid production well, no gas flux from reservoir to
!         well  ---
!
          IF( IT_CW(NCW).EQ.-2 ) THEN
            FGX_CWX = 0.D+0
          ELSE
!
!---        Gas flux from reservoir to well  ---
!
            FGX_CWX = WI_CWX*RKG(MF,N)*RHOG(MF,N)*DPGX/VISG(MF,N)
          ENDIF
!
!---      Water component gas flux from well to reservoir  ---
!
          FXW_CW(M,NWN) = FXW_CW(M,NWN) - FGX_CWX*XGW(MF,N)
!
!---      Gas component gas flux from well to reservoir  ---
!
          FXC_CW(1,M,NWN) = FXC_CW(1,M,NWN) - FGX_CWX*XGC(1,MF,N)
!
!---      Aqueous flux from well to reservoir  ---
!
          FLX_CWX = WI_CWX*RKL(1,MF,N)*RHOL(MF,N)*DPLX/VISL(MF,N)
!
!---      Water component aqueous fluxes from well to reservoir  ---
!
          FXW_CW(M,NWN) = FXW_CW(M,NWN) - FLX_CWX*XLW(MF,N)
!
!---      Nonaqueous-liquid flux from well to reservoir  ---
!
          FNX_CWX = WI_CWX*RKN(MF,N)*RHON(MF,N)*DPNX/VISN(MF,N)
!
!---      Gas component nonaqueous-liquid fluxes from well to 
!         reservoir  ---
!
          FXC_CW(1,M,NWN) = FXC_CW(1,M,NWN) - FNX_CWX*XNC(1,MF,N)
!
!---      Oil component nonaqueous-liquid fluxes from well to 
!         reservoir  ---
!
          FXC_CW(2,M,NWN) = FXC_CW(2,M,NWN) - FNX_CWX*XNC(2,MF,N)
!
!---      Store current coupled-well node location in previous
!         coupled-well node location  ---
!
          XPX(1) = XPX(2)
          YPX(1) = YPX(2)
          ZPX(1) = ZPX(2)
!
!---      Volumetric production well fluxes  ---
!
          IF( M.EQ.1 ) THEN
            IF( SL(2,N).GT.0.D+0 ) Q_CW(2,NWN) = FLX_CWX/RHOL(2,N)
            IF( SG(2,N).GT.0.D+0 ) Q_CW(3,NWN) = FGX_CWX/RHOG(2,N)
            IF( SN(2,N).GT.0.D+0 ) Q_CW(4,NWN) = FNX_CWX/RHON(2,N)
            Q_CW(1,NWN) = Q_CW(2,NWN) + Q_CW(3,NWN) + Q_CW(4,NWN)
          ENDIF
  200   CONTINUE
  300 CONTINUE
!
!---  Mass balance residuals for production type coupled well  ---
!
      NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
      MX = (NWFX*ISVC)+2
      RS_CW(1,NCW) = 0.D+0
      RS_CW(MX,NCW) = 0.D+0
      QM_CW(1,NCW) = 0.D+0
      QM_CW(3,NCW) = 0.D+0
      QM_CW(5,NCW) = 0.D+0
      DO 310 IGC = 1,NGC+2
        QM_CW((7+(IGC-1)*2),NCW) = 0.D+0
  310 CONTINUE
      QTM_CWX = 0.D+0
!
!---  Loop over coupled-well nodes from bottom to top  ---
!
      DO 400 NWN = ID_CW(4,NCW),ID_CW(3,NCW),-1
        RS_CW(1,NCW) = RS_CW(1,NCW) - FXW_CW(1,NWN)
        RS_CW(MX,NCW) = RS_CW(MX,NCW) - FXW_CW(ISVC+2,NWN)
        QTM_CWX = QTM_CWX - FXW_CW(1,NWN)
        QM_CW(3,NCW) = QM_CW(3,NCW) - FXW_CW(1,NWN)
        DO 320 IGC = 1,NGC+2
          RS_CW(1,NCW) = RS_CW(1,NCW) - FXC_CW(IGC,1,NWN)
          RS_CW(MX,NCW) = RS_CW(MX,NCW) - FXC_CW(IGC,ISVC+2,NWN)
          QTM_CWX = QTM_CWX - FXC_CW(IGC,1,NWN)
          QM_CW(1,NCW) = QM_CW(1,NCW) - FXC_CW(IGC,1,NWN)
          QM_CW((7+(IGC-1)*2),NCW) = QM_CW((7+(IGC-1)*2),NCW) - 
     &      FXC_CW(IGC,1,NWN)
  320   CONTINUE
  400 CONTINUE
!
!---  Loop over field nodes that contain coupled-well nodes  ---
!
      DO 470 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
        M1 = (NWF-ID_CW(5,NCW))*ISVC + 1
        DO 460 M2 = 1,ISVC
          M3 = M1+M2
          RS_CW(M3,NCW) = 0.D+0
!
!---      Loop over coupled-well nodes  ---
!
          DO 450 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---        If coupled-well node is within the current field
!           node, use incremented fluxes  ---
!
            IF( IWF_CW(NWF).EQ.IWN_CW(NWN) ) THEN
              RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXW_CW(M2+1,NWN)                
              DO 440 IGC = 1,NGC+2
                RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXC_CW(IGC,M2+1,NWN)
  440         CONTINUE
!
!---        If coupled-well node is outside the current field
!           node, use un-incremented fluxes  ---
!
            ELSE
              RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXW_CW(1,NWN)
              DO 442 IGC = 1,NGC+2
                RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXC_CW(IGC,1,NWN)
  442         CONTINUE
            ENDIF
  450     CONTINUE
  460   CONTINUE
  470 CONTINUE
  500 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROD_BO_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROD_COUP_WELL( NCW )
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
!     Production coupled well model
!     
!     Bottom-hole pressure controlled
!
!     Flux of energy, water mass, CO2 mass, CH4 mass, and petroleum
!     component masses from field nodes to coupled-well nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 April 2014.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE JACOB
      USE HYST
      USE GRID
      USE GLB_PAR
      USE FDVT
      USE FDVP
      USE FDVN
      USE FDVGC
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
      REAL*8 VAR_CWX(6+LNGC)
      INTEGER, SAVE :: IALLOC
      DATA IALLOC /0/
      INTEGER, SAVE, DIMENSION(:), ALLOCATABLE :: MCW,MFD
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROD_COUP_WELL'
!
!---  Black-oil option  ---
!
      IF( ISLC(68).EQ.0 ) THEN
        CALL PROD_BO_WELL( NCW )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compositional option  ---
!
      IF( IALLOC.EQ.0 ) THEN
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
   10   CONTINUE        
        IALLOC = 1
      ENDIF
!
!---  K-factors computed from fugacity coefficients  ---
!
      IFK = 0
!
!---  Zero production well fluxes ---
!
      QM_CW(1,NCW) = 0.D+0
      QM_CW(3,NCW) = 0.D+0
      QM_CW(5,NCW) = 0.D+0
      DO 11 IGC = 1,NGC+2
        QM_CW((7+(IGC-1)*2),NCW) = 0.D+0
   11 CONTINUE
!
!---  Loop over production well nodes from bottom to top  ---
!
      DO 30 NWN = ID_CW(4,NCW),ID_CW(3,NCW),-1
!
!---    Zero volumetric production well fluxes
!
!       Q_CW(1,NWN) - total volumetric flux, m^3/s
!       Q_CW(2,NWN) - aqueous volumetric flux, m^3/s
!       Q_CW(3,NWN) - gas volumetric flux, m^3/s
!       Q_CW(4,NWN) - nonaqueous-liquid volumetric flux, m^3/s  ---
!
        DO M = 1,4
          Q_CW(M,NWN) = 0.D+0
        ENDDO
!
!---    Loop over increment indices  ---
!
        DO 20 M = 1,ISVC+2
          DO 12 IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = 0.D+0
   12     CONTINUE
          FXE_CW(M,NWN) = 0.D+0
          FXW_CW(M,NWN) = 0.D+0
   20   CONTINUE
   30  CONTINUE
!
!---  Production well time interval ---
!
      TMZ = TM
      IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---  Cyclic time periods  ---
!
      IF( ICC_CW(NCW).EQ.1 ) TMZ = MOD( TM,VAR_CW(1,IM_CW(NCW),NCW) )
!
!---  Time prior to start time, production well is inactive, set 
!     well pressure to be in equilibrium with reservoir  ---
!
      IF( TMZ.LE.VAR_CW(1,1,NCW) ) THEN
        CALL EQUIL_COUP_WELL( NCW )
        GOTO 500
      ENDIF
!
!---  Time after start time, production well is active, with single
!     start time  ---
!
      IF( IM_CW(NCW).EQ.1 ) THEN
        DO 80 N = 2,NGC+8
          VAR_CWX(N) = VAR_CW(N,1,NCW)
   80   CONTINUE
         GOTO 122
!
!---  Time after start time, production well is active, with multiple
!     time points, use linear interpolation of well parameters 
!     between time points  ---
!
      ELSE
        DO 100 M = 2,IM_CW(NCW)
          IF( TMZ.LE.VAR_CW(1,M,NCW) ) THEN
            TD_CW = VAR_CW(1,M,NCW)-VAR_CW(1,M-1,NCW)
            DT_CW = MIN( VAR_CW(1,M,NCW)-TMZ,DT )
            TF_CW = (TMZ-VAR_CW(1,M-1,NCW))/TD_CW
            DO 90 N = 2,NGC+8
              VAR_CWX(N) = VAR_CW(N,M-1,NCW) + 
     &          TF_CW*(VAR_CW(N,M,NCW)-VAR_CW(N,M-1,NCW))
   90       CONTINUE
            GOTO 122
          ENDIF
  100   CONTINUE
      ENDIF
!
!---  Production well is inactive set well pressure to be in 
!     equilibrium with reservoir  ---
!
      CALL EQUIL_COUP_WELL( NCW )
      GOTO 500
  122 CONTINUE
!
!---  Set bottom-hole pressure  ---
!
      PL_CW(NCW) = VAR_CWX(3) - PATM
      P_CW(2,NCW) = PL_CW(NCW)
      P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
!
!---  Pressure controlled coupled well  ---
!
      ID_CW(8,NCW) = 1
!
!---  Well fluid density based on relative fluid mobility
!     at the previous time step  ---
!
      SRHOFX = 0.D+0
      SVOLX = 0.D+0
!
!---  Store bottom of coupled-well location in previous
!     coupled-well node location  ---
!
      XPX(1) = XTP_CW(2,ID_CW(2,NCW))
      YPX(1) = YTP_CW(2,ID_CW(2,NCW))
      ZPX(1) = ZTP_CW(2,ID_CW(2,NCW))
!
!---  Loop over the nodes in the coupled well from bottom to top ---
!
      DO 130 NWN = ID_CW(4,NCW),ID_CW(3,NCW),-1
        N = IWN_CW(NWN)
        I = ID(N)
        INVX = INV_CW(NWN)
        IZN = IZ(N)
!
!---    Coupled-well node centroids and projections ---
!
        XLX = PLX_CW(NWN)
        YLX = PLY_CW(NWN)
        ZLX = PLZ_CW(NWN)
        XPX(2) = 5.D-1*(XP_CW(2,NWN)+XP_CW(1,NWN))
        YPX(2) = 5.D-1*(YP_CW(2,NWN)+YP_CW(1,NWN))
        ZPX(2) = 5.D-1*(ZP_CW(2,NWN)+ZP_CW(1,NWN))
!
!---    Cylindrical coordinates with azimuthal symmetry,
!       centrally located wells  ---
!
        IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
          XPNX = 0.D+0
          YPNX = 0.D+0
          ZPNX = ZP(N)
!
!---    Cylindrical coordinates  ---
!
        ELSEIF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
          XPNX = XP(N)*COS(YP(N))
          YPNX = XP(N)*SIN(YP(N))
          ZPNX = ZP(N)
!
!---    Cartesian or boundary-fitted orthogonal coordinates  ---
!
        ELSE
          XPNX = XP(N)
          YPNX = YP(N)
          ZPNX = ZP(N)
        ENDIF
!
!---    Equivalent field node radius components  ---
!
        PERMX = MAX( PERM(1,IZN),1.D-20 )
        PERMY = MAX( PERM(2,IZN),1.D-20 )
        PERMZ = MAX( PERM(3,IZN),1.D-20 )
        RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---    Cylindrical coordinates with azimuthal symmetry,
!       centrally located wells  ---
!
        IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
          ROZ = RP(I)
!
!---      Specified pressure-equivalent radius  ---
!
          IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
            ROZ = PAR_CW(3,INVX)
          ENDIF
          RWX = MIN( RWX,9.999D-1*ROZ )
          PERMX = PERMRF(1,N)*PERM(1,IZN)
          PERMY = PERMRF(1,N)*PERM(2,IZN)
          WI_CWX = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &      (LOG(ROZ/RWX)+PAR_CW(1,INVX))
        ELSE
          PERMYZ = SQRT(PERMY/PERMZ)
          PERMZY = SQRT(PERMZ/PERMY)
          DXGFX = DXGF(N)/FF_CW(1,NCW)
          DYGFX = DYGF(N)*RP(I)/FF_CW(2,NCW)
          DZGFX = DZGF(N)/FF_CW(3,NCW)
          ROX = 2.8D-1*SQRT(PERMYZ*(DZGFX**2) + PERMZY*(DYGFX**2))
     &    /(SQRT(PERMYZ)+SQRT(PERMZY))
          PERMZX = SQRT(PERMZ/PERMX)
          PERMXZ = SQRT(PERMX/PERMZ)
          ROY = 2.8D-1*SQRT(PERMZX*(DXGFX**2) + PERMXZ*(DZGFX**2))
     &      /(SQRT(PERMZX)+SQRT(PERMXZ))
          PERMYX = SQRT(PERMY/PERMX)
          PERMXY = SQRT(PERMX/PERMY)
          ROZ = 2.8D-1*SQRT(PERMYX*(DXGFX**2) + PERMXY*(DYGFX**2))
     &      /(SQRT(PERMYX)+SQRT(PERMXY))
!
!---      Specified pressure-equivalent radius  ---
!
          IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
            ROX = PAR_CW(3,INVX)
            ROY = PAR_CW(3,INVX)
            ROZ = PAR_CW(3,INVX)
          ENDIF
!
!---      Well index components  ---
!
          PERMX = PERMRF(1,N)*PERM(1,IZN)
          PERMY = PERMRF(1,N)*PERM(2,IZN)
          PERMZ = PERMRF(1,N)*PERM(3,IZN)
          WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/
     &      (LOG(ROX/RWX)+PAR_CW(1,INVX))
          WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/
     &      (LOG(ROY/RWX)+PAR_CW(1,INVX))
          WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &      (LOG(ROZ/RWX)+PAR_CW(1,INVX))
          WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
        ENDIF
!
!---    Mass fluxes, positive into the well  ---
!
        DPGX = 1.D+0
        DPNX = 1.D+0
        DPLX = 1.D+0
!
!---    Liquid production well, no gas flux from reservoir to
!       well  ---
!
        IF( IT_CW(NCW).EQ.-2 ) THEN
          FGX_CWX = 0.D+0
        ELSE
!
!---      Gas flux from reservoir to well  ---
!
          FGX_CWX = WI_CWX*RKG(1,N)*RHOG(1,N)*DPGX/VISG(1,N)
        ENDIF
!
!---    Aqueous flux from well to reservoir  ---
!
        FLX_CWX = WI_CWX*RKL(1,1,N)*RHOL(1,N)*DPLX/VISL(1,N)
!
!---    Nonaqueous-liquid flux from well to reservoir  ---
!
        FNX_CWX = WI_CWX*RKN(1,N)*RHON(1,N)*DPNX/VISN(1,N)
!
!---    Total mass and volumetric flux into well from reservoir  ---
!
        SRHOFX = SRHOFX + FGX_CWX + FLX_CWX + FNX_CWX
        SVOLX = SVOLX + FGX_CWX/RHOG(1,N) + FLX_CWX/RHOL(1,N) + 
     &    FNX_CWX/RHON(1,N)
!
!---    Store current coupled-well node location in previous
!       coupled-well node location  ---
!
        XPX(1) = XPX(2)
        YPX(1) = YPX(2)
        ZPX(1) = ZPX(2)
  130 CONTINUE
!
!---  Dynamic average well-fluid density  ---
!
      IF( SRHOFX.GT.EPSL ) THEN
        RHOF_CWX = SRHOFX/SVOLX
!
!---  Static average well-fluid density  ---
!
      ELSE
        SRHOFX = 0.D+0
        SVOLX = 0.D+0
!
!---  Loop over the nodes in the coupled well from bottom to top ---
!
        DO 140 NWN = ID_CW(4,NCW),ID_CW(3,NCW),-1
          N = IWN_CW(NWN)
!
!---      Liquid production well, no gas flux from reservoir to
!         well  ---
!
          IF( IT_CW(NCW).EQ.-2 ) THEN
            SRHOFX = SRHOFX + VOL(N)*PORD(1,N)* 
     &        (RHOL(1,N)/VISL(1,N) + RHON(1,N)/VISN(1,N))
            SVOLX = SVOLX + VOL(N)*PORD(1,N)* 
     &        (1.D+0/VISL(1,N) +  1.D+0/VISN(1,N))
          ELSE
            SRHOFX = SRHOFX + VOL(N)*PORD(1,N)*
     &       (RHOG(1,N)/VISG(1,N) + RHOL(1,N)/VISL(1,N) + 
     &        RHON(1,N)/VISN(1,N))
            SVOLX = SVOLX + VOL(N)*PORD(1,N)*
     &       (1.D+0/VISG(1,N) + 1.D+0/VISL(1,N) + 1.D+0/VISN(1,N))
          ENDIF
  140   CONTINUE
        IF( SRHOFX.GT.EPSL ) THEN
          RHOF_CWX = SRHOFX/SVOLX
!
!---    Reference well-fluid density  ---
!
        ELSE
          RHOF_CWX = RHORL
        ENDIF
      ENDIF
!
!---  Record well fluid density for output  ---
!
      RHOF_CW(NCW) = RHOF_CWX
!
!---  Loop over increment indices ---
!
      DO 300 M = 1,ISVC+2
        MW = MCW(M)
        MF = MFD(M)
        N = IWN_CW(ID_CW(3,NCW))
        P_CWX = P_CW(MW,NCW)
!
!---    Store bottom of coupled-well location in previous
!       coupled-well node location  ---
!
        XPX(1) = XTP_CW(2,ID_CW(2,NCW))
        YPX(1) = YTP_CW(2,ID_CW(2,NCW))
        ZPX(1) = ZTP_CW(2,ID_CW(2,NCW))
!
!---    Loop over the nodes in the coupled well from bottom to top ---
!
        SRHOFX = 0.D+0
        SVOLX = 0.D+0
        DO 200 NWN = ID_CW(4,NCW),ID_CW(3,NCW),-1
          N = IWN_CW(NWN)
          I = ID(N)
          INVX = INV_CW(NWN)
          IZN = IZ(N)
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
          P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOF_CWX
!
!---      Adjust the reservoir pressure to the coupled-well node
!         centroid  ---
!
!---      Phase condition #1  ---
!
!         SL = 1.0, SG+SN = 0.0
!
          IF( MOD(NPHAZ(2,N),10).EQ.1 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #2  ---
!
!         SL < 1.0, ST > 0.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.2 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #3  ---
!
!         SL < 1.0, ST > 0.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.3 ) THEN
            PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(MF,N)
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
!
!---      Phase condition #4  ---
!
!         SL = 0.0, ST = 1.0, SG = 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.4 ) THEN
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PTX = MAX( PLFX,PNFX )
            CPGTX = MAX( (PG(MF,N)-MAX( PL(MF,N),PN(MF,N) )),0.D+0 )
            PGFX = PTX + CPGTX
!
!---      Phase condition #5  ---
!
!         SL = 0.0, ST = 1.0, SG > 0.0
!
          ELSEIF( MOD(NPHAZ(2,N),10).EQ.5 ) THEN
            PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(MF,N)
            PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
            PNFX = PN(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHON(MF,N)
          ENDIF
!
!---      Equivalent field node radius components  ---
!
          PERMX = MAX( PERM(1,IZN),1.D-20 )
          PERMY = MAX( PERM(2,IZN),1.D-20 )
          PERMZ = MAX( PERM(3,IZN),1.D-20 )
          RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---      Cylindrical coordinates with azimuthal symmetry,
!         centrally located wells  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            ROZ = RP(I)
!
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROZ = PAR_CW(3,INVX)
            ENDIF
            RWX = MIN( RWX,9.999D-1*ROZ )
            PERMX = PERMRF(MF,N)*PERM(1,IZN)
            PERMY = PERMRF(MF,N)*PERM(2,IZN)
            WI_CWX = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &        (LOG(ROZ/RWX)+PAR_CW(1,INVX))
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
!---        Specified pressure-equivalent radius  ---
!
            IF( PAR_CW(3,INVX)/EPSL.GT.EPSL ) THEN
              ROX = PAR_CW(3,INVX)
              ROY = PAR_CW(3,INVX)
              ROZ = PAR_CW(3,INVX)
            ENDIF
!
!---        Well index components  ---
!
            PERMX = PERMRF(MF,N)*PERM(1,IZN)
            PERMY = PERMRF(MF,N)*PERM(2,IZN)
            PERMZ = PERMRF(MF,N)*PERM(3,IZN)
            WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/
     &        (LOG(ROX/RWX)+PAR_CW(1,INVX))
            WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/
     &        (LOG(ROY/RWX)+PAR_CW(1,INVX))
            WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &        (LOG(ROZ/RWX)+PAR_CW(1,INVX))
            WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
          ENDIF
!
!---      Mass fluxes, positive into the well  ---
!
          DPGX = MAX( PGFX-P_CWX,0.D+0 )
          DPNX = MAX( PGFX-P_CWX,0.D+0 )
          DPLX = MAX( PGFX-P_CWX,0.D+0 )
!
!---      Zero fluxes of water, energy, and nonaqueous components
!         from well to reservoir  ---
!
          FXE_CW(M,NWN) = 0.D+0
          FXW_CW(M,NWN) = 0.D+0
          DO 150 IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = 0.D+0
  150     CONTINUE
!
!---      Liquid production well, no gas flux from reservoir to
!         well  ---
!
          IF( IT_CW(NCW).EQ.-2 ) THEN
            FGX_CWX = 0.D+0
          ELSE
!
!---        Gas flux from reservoir to well  ---
!
            FGX_CWX = WI_CWX*RKG(MF,N)*RHOG(MF,N)*DPGX/VISG(MF,N)
          ENDIF
!
!---      Energy, water and oil component fluxes
!         from well to reservoir via gas phase  ---
!
          IF( ISLC(30).EQ.0 )
     &      FXE_CW(M,NWN) = FXE_CW(M,NWN) - FGX_CWX*HG(MF,N)
          FXW_CW(M,NWN) = FXW_CW(M,NWN) - FGX_CWX*XGW(MF,N)
          DO 152 IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = FXC_CW(IGC,M,NWN) - 
     &        FGX_CWX*XGC(IGC,MF,N)
  152     CONTINUE
!
!---      Aqueous flux from reservoir to well  ---
!
          FLX_CWX = WI_CWX*RKL(1,MF,N)*RHOL(MF,N)*DPLX/VISL(MF,N)
!
!---      Energy, water and oil component fluxes
!         from well to reservoir via aqueous phase  ---
!
          IF( ISLC(30).EQ.0 )
     &      FXE_CW(M,NWN) = FXE_CW(M,NWN) - FLX_CWX*HL(MF,N)
          FXW_CW(M,NWN) = FXW_CW(M,NWN) - FLX_CWX*XLW(MF,N)
          FXC_CW(1,M,NWN) = FXC_CW(1,M,NWN) - FLX_CWX*XLA(MF,N)
!
!---      Nonaqueous-liquid flux from reservoir to well  ---
!
          FNX_CWX = WI_CWX*RKN(MF,N)*RHON(MF,N)*DPNX/VISN(MF,N)
!
!---      Energy, water and oil component fluxes
!         from well to reservoir via nonaqueous-liquid phase  ---
!
          IF( ISLC(30).EQ.0 )
     &      FXE_CW(M,NWN) = FXE_CW(M,NWN) - FNX_CWX*HN(MF,N)
          DO 154 IGC = 1,NGC+2
            FXC_CW(IGC,M,NWN) = FXC_CW(IGC,M,NWN) - 
     &        FNX_CWX*XNC(IGC,MF,N)
  154     CONTINUE
!
!---      Total mass and volumetric flux into well from reservoir  ---
!
          SRHOFX = SRHOFX + FGX_CWX + FLX_CWX + FNX_CWX
          SVOLX = SVOLX + FGX_CWX/RHOG(MF,N) + FLX_CWX/RHOL(MF,N) + 
     &      FNX_CWX/RHON(MF,N)
!
!---      Store current coupled-well node location in previous
!         coupled-well node location  ---
!
          XPX(1) = XPX(2)
          YPX(1) = YPX(2)
          ZPX(1) = ZPX(2)
!
!---      Volumetric production well fluxes  ---
!
          IF( M.EQ.1 ) THEN
            IF( SL(2,N).GT.0.D+0 ) Q_CW(2,NWN) = FLX_CWX/RHOL(2,N)
            IF( SG(2,N).GT.0.D+0 ) Q_CW(3,NWN) = FGX_CWX/RHOG(2,N)
            IF( SN(2,N).GT.0.D+0 ) Q_CW(4,NWN) = FNX_CWX/RHON(2,N)
            Q_CW(1,NWN) = Q_CW(2,NWN) + Q_CW(3,NWN) + Q_CW(4,NWN)
          ENDIF
  200   CONTINUE
  300 CONTINUE
!
!---  Mass balance residuals for production type coupled well  ---
!
      NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
      MX = (NWFX*ISVC)+2
      RS_CW(1,NCW) = 0.D+0
      RS_CW(MX,NCW) = 0.D+0
      QM_CW(1,NCW) = 0.D+0
      QM_CW(3,NCW) = 0.D+0
      QM_CW(5,NCW) = 0.D+0
      DO 310 IGC = 1,NGC+2
        QM_CW((7+(IGC-1)*2),NCW) = 0.D+0
  310 CONTINUE
      QTM_CWX = 0.D+0
!
!---  Loop over coupled-well nodes from bottom to top  ---
!
      DO 400 NWN = ID_CW(4,NCW),ID_CW(3,NCW),-1
        RS_CW(1,NCW) = RS_CW(1,NCW) - FXW_CW(1,NWN)
        RS_CW(MX,NCW) = RS_CW(MX,NCW) - FXW_CW(ISVC+2,NWN)
        QTM_CWX = QTM_CWX - FXW_CW(1,NWN)
        QM_CW(3,NCW) = QM_CW(3,NCW) - FXW_CW(1,NWN)
        IF( ISLC(30).EQ.0 )
     &    QM_CW(5,NCW) = QM_CW(5,NCW) - FXE_CW(1,NWN)
        DO 320 IGC = 1,NGC+2
          RS_CW(1,NCW) = RS_CW(1,NCW) - FXC_CW(IGC,1,NWN)
          RS_CW(MX,NCW) = RS_CW(MX,NCW) - FXC_CW(IGC,ISVC+2,NWN)
          QTM_CWX = QTM_CWX - FXC_CW(IGC,1,NWN)
          QM_CW(1,NCW) = QM_CW(1,NCW) - FXC_CW(IGC,1,NWN)
          QM_CW((7+(IGC-1)*2),NCW) = QM_CW((7+(IGC-1)*2),NCW) - 
     &      FXC_CW(IGC,1,NWN)
  320   CONTINUE
  400 CONTINUE
!
!---  Loop over field nodes that contain coupled-well nodes  ---
!
      DO 470 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
        M1 = (NWF-ID_CW(5,NCW))*ISVC + 1
        DO 460 M2 = 1,ISVC
          M3 = M1+M2
          RS_CW(M3,NCW) = 0.D+0
!
!---      Loop over coupled-well nodes  ---
!
          DO 450 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---        If coupled-well node is within the current field
!           node, use incremented fluxes  ---
!
            IF( IWF_CW(NWF).EQ.IWN_CW(NWN) ) THEN
              RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXW_CW(M2+1,NWN)                
              DO 440 IGC = 1,NGC+2
                RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXC_CW(IGC,M2+1,NWN)
  440         CONTINUE
!
!---        If coupled-well node is outside the current field
!           node, use un-incremented fluxes  ---
!
            ELSE
              RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXW_CW(1,NWN)
              DO 442 IGC = 1,NGC+2
                RS_CW(M3,NCW) = RS_CW(M3,NCW) - FXC_CW(IGC,1,NWN)
  442         CONTINUE
            ENDIF
  450     CONTINUE
  460   CONTINUE
  470 CONTINUE
  500 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROD_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROJ_COUP_WELL( N,NWN )
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
      USE GLB_PAR
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
      REAL*8 ROT_MAT(3,3),PLGCX(3),PLLCX(3)
      REAL*8 XVX(4),YVX(4),ZVX(4)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROJ_COUP_WELL'
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
      PLGCX(1) = ABS(XP_CW(2,NWN)-XP_CW(1,NWN))
      PLGCX(2) = ABS(YP_CW(2,NWN)-YP_CW(1,NWN))
      PLGCX(3) = ABS(ZP_CW(2,NWN)-ZP_CW(1,NWN))
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
      PLX_CW(NWN) = PLLCX(1)*WLGCX/WLLCX
      PLY_CW(NWN) = PLLCX(2)*WLGCX/WLLCX
      PLZ_CW(NWN) = PLLCX(3)*WLGCX/WLLCX
!
!---  End of PROJ_COUP_WELL group  ---
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
!     STOMP-EOR
!
!     Reads the coupled well card.
!     
!     FF_CW(1,LN_CW) - x-direction well fraction factor for well
!     FF_CW(2,LN_CW) - y-direction well fraction factor for well
!     FF_CW(3,LN_CW) - z-direction well fraction factor for well
!     ICC_CW(LN_CW) - cyclic well time index for well
!     IM_CW(LN_CW) - number of time points for well
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
!     ISOLU_CW(LSOLU_CW,LN_CW) - solute number
!     VARC_CW(LSOLU_CW,LWT_CW,LN_CW) - solute conc. in well fluid, 1/m^3
!     ISPC_CW(LSPC_CW,LN_CW) - species number
!     VARSP_CW(LSPC_CW,LWT_CW,LN_CW) - species conc. in well fluid, 
!       mol/m^3
!
!----------------------Compositional Option----------------------------!
!
!     ITS_CW(LN_CW) - type index for injection well state
!     Injection type wells
!     122 - State #1, salt rel. sat., CO2 rel. sat.
!     123 - State #1, salt rel. sat., CO2 mass frac.
!     132 - State #1, salt mass frac., CO2 rel. sat.
!     133 - State #1, salt mass frac., CO2 mass frac.
!     220 - State #2, salt rel. sat.
!     230 - State #2, salt mass frac.
!     300 - State #3
!
!     JM_CW(LN_CW) - location of the well equation in the Jacobian
!     N_CW - number of coupled wells
!     PAR_CW(1,LWI_CW) - skin factor for well interval
!     VAR_CW(1,LWT_CW,LN_CW) - well time, s
!     VAR_CW(2,LWT_CW,LN_CW) - mass rate, kg/s
!     VAR_CW(3,LWT_CW,LN_CW) - press. limit, Pa
!     VAR_CW(4,LWT_CW,LN_CW) - temperature, C
!     VAR_CW(5,LWT_CW,LN_CW) - aqueous CO2 relative saturation
!     VAR_CW(5,LWT_CW,LN_CW) - aqueous CO2 mass fraction
!     VAR_CW(5,LWT_CW,LN_CW) - aqueous saturation
!     VAR_CW(5,LWT_CW,LN_CW) - water vapor relative saturation
!     VAR_CW(6,LWT_CW,LN_CW) - aqueous salt relative saturation
!     VAR_CW(6,LWT_CW,LN_CW) - aqueous salt mass fraction
!     VAR_CW(7,LWT_CW,LN_CW) - CO2 nonaqueous mole fraction
!     VAR_CW(8,LWT_CW,LN_CW) - CH4 nonaqueous mole fraction
!     VAR_CW(8+I,LWT_CW,LN_CW) - petroleum comp. I nonaqueous mole frac.
!     PAR_CW(2,LWI_CW) - well bore radius for well interval, m
!     XTP_CW(2,LWI_CW) - x-transition points for well interval, m
!     YTP_CW(2,LWI_CW) - y-transition points for well interval, m
!     ZTP_CW(2,LWI_CW) - z-transition points for well interval, m
!     PAR_CW(3,LWI_CW) - pressure equivalent radius, m
!
!----------------------Black-Oil Option----------------------------!
!
!     ITS_CW(LN_CW) - type index for injection well state
!     Injection type wells
!     122 - State #1, salt rel. sat.
!     123 - State #1, salt rel. sat.
!     220 - State #2, salt rel. sat.
!     230 - State #2, salt mass frac.
!     300 - State #3
!
!     JM_CW(LN_CW) - location of the well equation in the Jacobian
!     N_CW - number of coupled wells
!     PAR_CW(1,LWI_CW) - skin factor for well interval
!     VAR_CW(1,LWT_CW,LN_CW) - well time, s
!     VAR_CW(2,LWT_CW,LN_CW) - mass rate, kg/s
!     VAR_CW(3,LWT_CW,LN_CW) - press. limit, Pa
!     VAR_CW(4,LWT_CW,LN_CW) - temperature, C
!     VAR_CW(5,LWT_CW,LN_CW) - aqueous saturation
!     VAR_CW(5,LWT_CW,LN_CW) - water vapor relative saturation
!     VAR_CW(6,LWT_CW,LN_CW) - aqueous salt relative saturation
!     VAR_CW(6,LWT_CW,LN_CW) - aqueous salt mass fraction
!     VAR_CW(7,LWT_CW,LN_CW) - gas mass fraction of gas + oil
!     PAR_CW(2,LWI_CW) - well bore radius for well interval, m
!     XTP_CW(2,LWI_CW) - x-transition points for well interval, m
!     YTP_CW(2,LWI_CW) - y-transition points for well interval, m
!     ZTP_CW(2,LWI_CW) - z-transition points for well interval, m
!     PAR_CW(3,LWI_CW) - pressure equivalent radius, m
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 March 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FD
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE GRID
      USE GLB_PAR
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
      SAVE FORM1
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
        IF( INDEX(ADUM(1:),'injection').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Coupled Mass-Rate Injection Well'
          IT_CW(NCW) = 1
          WNM_CW(NCW) = 'IW'
          ICX = ICOUNT(NCW)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(WNM_CW(NCW)(3:3+ICX-1),FORM1) NCW
        ELSEIF( INDEX(ADUM(1:),'injection').NE.0 .AND.
     &    INDEX(ADUM(1:),'volu').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Coupled Volumetric-Rate Injection Well'
          IT_CW(NCW) = 2
          WNM_CW(NCW) = 'IW'
          ICX = ICOUNT(NCW)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(WNM_CW(NCW)(3:3+ICX-1),FORM1) NCW
        ELSEIF( INDEX(ADUM(1:),'production').NE.0 .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Coupled Liquid Production Well'
          IT_CW(NCW) = -2
          WNM_CW(NCW) = 'PW'
          ICX = ICOUNT(NCW)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(WNM_CW(NCW)(3:3+ICX-1),FORM1) NCW
        ELSEIF( INDEX(ADUM(1:),'production').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Coupled Fluid Production Well'
          IT_CW(NCW) = -1
          WNM_CW(NCW) = 'PW'
          ICX = ICOUNT(NCW)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(WNM_CW(NCW)(3:3+ICX-1),FORM1) NCW
!
!---    Unrecognized well type  ---
!
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
!---    Check for solutes or species, if solute transport
!       or reactive transport is active  ---
!
        IF( IEQC.NE.0 .OR. ISLC(40).EQ.1 ) THEN
          NC = 0
          DO
            CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
            ISX = ISTART
            ICX = ICOMMA
            IF( INDX.EQ.0 ) EXIT
!
!---        Check for solute or species names  ---
!
            ADUM(1:) = ' '
            VARB = 'Solute or Species Name'
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
#ifdef ecke
!
!---        Aqueous species  ---
!
            DO M = 1,NSPL
              IF( SPNML(M).EQ.ADUM ) THEN
                NC = NC + 1
                ISPC_CW(NC,NCW) = M
                IFIND = 1
                EXIT
              ENDIF
            ENDDO
!
!---        Gas species  ---
!
            DO M = 1,NSPG
              IF( SPNMG(M).EQ.ADUM ) THEN
                NC = NC + 1
                ISPC_CW(NC,NCW) = M
                IFIND = 1
                EXIT
              ENDIF
            ENDDO
#endif
            IF( IFIND.EQ.0 ) EXIT
          ENDDO
          NSP_CW(NCW) = NC
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
!
!---      Check for pressure equivalent radius  ---
!
          VARB = 'Pressure-Equivalent Radius and Effective ' // 
     &      'Permeability-Thickness Product'
          CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
          IF( INDX.EQ.1 ) THEN
!
!---        Read pressure-equivalent radius  ---
!
            VARB = 'Pressure-Equivalent Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(3,NICW))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_CW(3,NICW)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_CW(3,NICW),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_CW(3,NICW),', m)'
          ENDIF
  100   CONTINUE
!
!---    Read number of well time periods  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Number of Well Time Periods'
        CALL RDINT(ISTART,ICOMMA,CHDUM,IM_CW(NCW))
!
!---    Check number of well time points parameter  ---
!
        IF( IM_CW(NCW).GT.LWT_CW ) THEN
          INDX = 5
          CHMSG = 'Number of Well Time Periods > LWT_CW'
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
!---    Loop over the number of well time periods  ---
!
        NC = 0
        DO 300 M = 1,IM_CW(NCW)
!
!---      Injection well, read number of time points for the time
!         period and state of the injection well for the time
!         period  ---
!
          IF( IT_CW(NCW).GT.0 ) THEN
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
            ISTART = 1
            VARB = 'Number of Time Points in Time Period'
            CALL RDINT(ISTART,ICOMMA,CHDUM,IMP_CW(M,NCW))
            VARB = 'Coupled Injection Well State'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM(1))
!
!---        Aqueous saturated injection well state  ---
!
            IF( INDEX(BDUM(1)(1:),'cw1').NE.0 ) THEN
              ITS_CW(M,NCW) = 100
              WRITE(IWR,'(2X,A)') 'Coupled Well State #1'
!
!---          Compositional option  ---
!
              IF( ISLC(68).NE.0 ) THEN
!
!---          Read aqueous CO2 option  ---
!
                VARB = 'Aqueous CO2 Option'
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(2))
                IF( INDEX(BDUM(2)(1:),'co2').NE.0 .AND.
     &            INDEX(BDUM(2)(1:),'rel').NE.0 .AND.
     &            INDEX(BDUM(2)(1:),'sat').NE.0 ) THEN
                  ITS_CW(M,NCW) = ITS_CW(M,NCW) + 2
                  WRITE(IWR,'(2X,2A)') 'Aqueous CO2: Relative ' // 
     &              'Saturation'
                ELSEIF( INDEX(BDUM(2)(1:),'co2').NE.0 .AND.
     &            INDEX(BDUM(2)(1:),'mass').NE.0 .AND.
     &            INDEX(BDUM(2)(1:),'frac').NE.0 ) THEN
                  ITS_CW(M,NCW) = ITS_CW(M,NCW) + 3
                  WRITE(IWR,'(2X,2A)') 'Aqueous CO2: Mass Fraction'
                ELSE
                  INDX = 4
                  CHMSG = 'Unrecognized Aqueous CO2 Option: '
     &              // BDUM(2)
                  CALL WRMSGS( INDX )
                ENDIF
              ENDIF
!
!---          Read aqueous salt option  ---
!
              VARB = 'Aqueous Salt Option'
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(3))
              IF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &          INDEX(BDUM(3)(1:),'rel').NE.0 .AND.
     &          INDEX(BDUM(3)(1:),'sat').NE.0 ) THEN
                ITS_CW(M,NCW) = ITS_CW(M,NCW) + 20
                WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Relative Saturation'
              ELSEIF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &          INDEX(BDUM(3)(1:),'mass').NE.0 .AND.
     &          INDEX(BDUM(3)(1:),'frac').NE.0 ) THEN
                ITS_CW(M,NCW) = ITS_CW(M,NCW) + 30
                WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Mass Fraction'
              ELSE
                INDX = 4
                CHMSG = 'Unrecognized Aqueous Salt Option: '
     &            // BDUM(3)
                CALL WRMSGS( INDX )
              ENDIF
!
!---        Three-phase injection well state   ---
!
            ELSEIF( INDEX(BDUM(1)(1:),'cw2').NE.0 ) THEN
              ITS_CW(M,NCW) = 200
              WRITE(IWR,'(2X,A)') 'Coupled Well State #2'
!
!---          Read aqueous salt option  ---
!
              VARB = 'Aqueous Salt Option'
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(3))
              IF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &          INDEX(BDUM(3)(1:),'rel').NE.0 .AND.
     &          INDEX(BDUM(3)(1:),'sat').NE.0 ) THEN
                ITS_CW(M,NCW) = ITS_CW(M,NCW) + 20
                WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Relative Saturation'
              ELSEIF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &          INDEX(BDUM(3)(1:),'mass').NE.0 .AND.
     &          INDEX(BDUM(3)(1:),'frac').NE.0 ) THEN
                ITS_CW(M,NCW) = ITS_CW(M,NCW) + 30
                WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Mass Fraction'
              ELSE
                INDX = 4
                CHMSG = 'Unrecognized Aqueous Salt Option: '
     &            // BDUM(3)
                CALL WRMSGS( INDX )
              ENDIF
!
!---        Gas and nonaqueous saturated injection well state  ---
!
            ELSEIF( INDEX(BDUM(1)(1:),'cw3').NE.0 ) THEN
              ITS_CW(M,NCW) = 300
              WRITE(IWR,'(2X,A)') 'Coupled Well State #3'
!
!---        Unrecognized injection well state  ---
!
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Coupled Injection Well State: ' // 
     &          BDUM(1)(1:NCHB)
              CALL WRMSGS( INDX )
            ENDIF
!
!---      Production well, set the number of time points for the
!         time period to 1  ---
!
          ELSE
            IMP_CW(M,NCW) = 1       
          ENDIF
!
!---      Loop over the number of time points in time period  ---
!
          DO 200 MX = 1,IMP_CW(M,NCW)
            NC = NC + 1
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
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(1,NC,NCW))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR_CW(1,NC,NCW)
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
!---        Injection type well  ---
!
            IF( IT_CW(NCW).GT.0 ) THEN
!
!---          Read injection temperature, unless isothermal option  ---
!
              IF( ISLC(30).EQ.0 ) THEN
                VARB = 'Temperature'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(4,NC,NCW))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',VAR_CW(4,NC,NCW)
                IUNK = 1
                CALL RDUNIT(UNTS,VAR_CW(4,NC,NCW),INDX)
                WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(4,NC,NCW),', C)'
              ENDIF
!
!---          Read injection mass rate  ---
!
              IF( IT_CW(NCW).EQ.1 ) THEN
                VARB = 'Injection Mass Rate'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(2,NC,NCW))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',VAR_CW(2,NC,NCW)
                IUNKG = 1
                IUNS = -1
                CALL RDUNIT(UNTS,VAR_CW(2,NC,NCW),INDX)
                WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(2,NC,NCW),
     &            ', kg/s)'
!
!---          Read injection volumetric rate  ---
!
              ELSEIF( IT_CW(NCW).EQ.2 ) THEN
                VARB = 'Injection Volumetric Rate'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(2,NC,NCW))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',VAR_CW(2,NC,NCW)
                IUNM = 3
                IUNS = -1
                CALL RDUNIT(UNTS,VAR_CW(2,NC,NCW),INDX)
                WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(2,NC,NCW),
     &            ', m^3/s)'
              ENDIF
!
!---          Read maximum well-top pressure  ---
!
              VARB = 'Maximum Well-Top Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(3,NC,NCW))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_CW(3,NC,NCW)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_CW(3,NC,NCW),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(3,NC,NCW),', Pa)'
!
!---          Well state #1  ---
!
!             SL = 1.0
!             SG+SN = 0.0
!
!             Declared variables:
!
!             VAR_CW(5,NC,NCW) - aqueous CO2 relative saturation, or
!                                aqueous CO2 mass fraction
!             VAR_CW(6,NC,NCW) - aqueous salt relative saturation, or
!                                aqueous salt mass fraction
!
              IF( ITS_CW(M,NCW)/100.EQ.1 ) THEN
!
!---            Compositional option  ---
!
                IF( ISLC(68).NE.0 ) THEN
!
!---              Read aqueous CO2  ---
!
                  IF( MOD(ITS_CW(M,NCW),10).EQ.2 ) THEN
                    VARB = 'Aqueous CO2 Relative Saturation'
                    CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(5,NC,NCW))
                    WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &                ': ',VAR_CW(5,NC,NCW)
                  ELSEIF( MOD(ITS_CW(M,NCW),10).EQ.3 ) THEN
                    VARB = 'Aqueous CO2 Mass Fraction'
                    CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(5,NC,NCW))
                    WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &                ': ',VAR_CW(5,NC,NCW)
                  ENDIF
                ENDIF
!
!---            Read aqueous salt  ---
!
                IF( MOD((ITS_CW(M,NCW)/10),10).EQ.2 ) THEN
                  VARB = 'Aqueous Salt Relative Saturation'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(6,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(6,NC,NCW)
                ELSEIF( MOD((ITS_CW(M,NCW)/10),10).EQ.3 ) THEN
                  VARB = 'Aqueous Salt Mass Fraction'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(6,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(6,NC,NCW)
                ENDIF
!
!---          Well state #2  ---
!
!             SL < 1.0
!             SG+SN > 0.0
!
!             Declared variables:
!
!             VAR_CW(5,NC,NCW) - aqueous saturation
!             VAR_CW(6,NC,NCW) - aqueous salt relative saturation, or
!                                 aqueous salt mass fraction
!             VAR_CW(6+NGC,NC,NCW) - total mole fraction of petro. comp.
!
              ELSEIF( ITS_CW(M,NCW)/100.EQ.2 ) THEN
!
!---            Read aqueous saturation  ---
!
                VARB = 'Aqueous Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(5,NC,NCW))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_CW(5,NC,NCW)
!
!---            Read aqueous salt  ---
!
                IF( MOD((ITS_CW(M,NCW)/10),10).EQ.2 ) THEN
                  VARB = 'Aqueous Salt Relative Saturation'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(6,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(6,NC,NCW)
                ELSEIF( MOD((ITS_CW(M,NCW)/10),10).EQ.3 ) THEN
                  VARB = 'Aqueous Salt Mass Fraction'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(6,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(6,NC,NCW)
                ENDIF
!
!---            Black-oil option  ---
!
                IF( ISLC(68).EQ.0 ) THEN
                  VARB = 'Gas Mass Fraction of Gas + Oil'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(7,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(7,NC,NCW)
!
!---            Compositional option  ---
!
                ELSE
!
!---              Read total mole fraction of petroleum components  ---
!
                  SZMCX = 0.D+0
                  DO 130 IGC = 1,NGC+2
                    JGC = 6 + IGC
                    VARB = 'Total Mole Fraction of ' // GCNM(IGC)
                    CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(JGC,NC,NCW))
                    WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &                ': ',VAR_CW(JGC,NC,NCW)
                    SZMCX = SZMCX + VAR_CW(JGC,NC,NCW)
  130             CONTINUE
                  IF( ABS(1.D+0-SZMCX).GT.EPSL ) THEN
                    INDX = 9
                    RLMSG = SZMCX
                    CHMSG = 'Summation of Total Mole Fraction ' // 
     &                'of Petroleum Components'
                    CALL WRMSGS( INDX )
                  ELSE
                    DO 132 IGC = 1,NGC+2
                      VAR_CW(JGC,NC,NCW) = VAR_CW(JGC,NC,NCW)/SZMCX
  132               CONTINUE
                  ENDIF
                ENDIF
!
!---          Well state #3  ---
!
!             SL = 0.0
!             SG+SN = 1.0
!
!             Declared variables:
!
!             VAR_CW(5,NC,NCW) - water vapor relative saturation
!             VAR_CW(6+,NC,NCW) - total mole fraction of petroleum comp.
!            
              ELSEIF( ITS_CW(M,NCW)/100.EQ.3 ) THEN
!
!---            Read water vapor relative saturation  ---
!
                VARB = 'Water Vapor Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(5,NC,NCW))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_CW(5,NC,NCW)
!
!---            Black-oil option  ---
!
                IF( ISLC(68).EQ.0 ) THEN
                  VARB = 'Gas Mass Fraction of Gas + Oil'
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(7,NC,NCW))
                  WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &              ': ',VAR_CW(7,NC,NCW)
!
!---            Compositional option  ---
!
                ELSE
!
!---              Read total mole fraction of petroleum components  ---
!
                  SZMCX = 0.D+0
                  DO 140 IGC = 1,NGC+2
                    JGC = 6 + IGC
                    VARB = 'Total Mole Fraction of ' // GCNM(IGC)
                    CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(JGC,NC,NCW))
                    WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &                ': ',VAR_CW(JGC,NC,NCW)
                    SZMCX = SZMCX + VAR_CW(JGC,NC,NCW)
  140             CONTINUE
                  IF( ABS(1.D+0-SZMCX).GT.EPSL ) THEN
                    INDX = 9
                    RLMSG = SZMCX
                    CHMSG = 'Summation of Total Mole Fraction ' // 
     &                'of Petroleum Components'
                    CALL WRMSGS( INDX )
                  ELSE
                    DO 142 IGC = 1,NGC+2
                      VAR_CW(JGC,NC,NCW) = VAR_CW(JGC,NC,NCW)/SZMCX
  142               CONTINUE
                  ENDIF
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
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VARC_CW(NSL_CW,NC,NCW))
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &              UNTS(1:NCH),': ',VARC_CW(NSL_CW,NC,NCW)
                  INDX = 0
                  IUNM = -3
                  CALL RDUNIT(UNTS,VARC_CW(NSL_CW,NC,NCW),INDX)
                  WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &              VARC_CW(NSL_CW,NC,NCW),', 1/m^3)'
                ENDDO
              ENDIF
#ifdef ecke
!
!---          Simulation with reactive species transport  ---
!
              IF( ISLC(40).EQ.1 ) THEN
                ISKIP = 1
                DO NSPCW = 1,LSPC_CW
!
!---              Species active in well  ---
!
                  IF( ISPC_CW(NSPCW,NCW).EQ.0 ) EXIT
                  NSP = ISPC_CW(NSPCW,NCW)
!
!---              Well state #1  ---
!
!                 SL = 1.0
!                 SG+SN = 0.0
!
!                 Aqueous species only
!
                  IF( ITS_CW(M,NCW)/100.EQ.1 ) THEN
                    IF( NSP.LE.NSPL ) ISKIP = 0
!
!---              Well state #2  ---
!
!                 SL < 1.0
!                 SG+SN > 0.0
!
!                 Aqueous, gas, and nonaqueous-liquid species
!
                  ELSEIF( ITS_CW(M,NCW)/100.EQ.2 ) THEN
                    IF( NSP.LE.NSPL .OR. (NSP.GT.NSPL+NSPS .AND. 
     &                NSP.LE.NSPL+NSPS+NSPG) ) ISKIP = 0
!
!---              Well state #3  ---
!
!                 SL = 0.0
!                 SG+SN = 1.0
!
!                 Gas and nonaqueous-liquid species only
!            
                  ELSEIF( ITS_CW(M,NCW)/100.EQ.3 ) THEN
                    IF( NSP.GT.NSPL+NSPS .AND. 
     &                NSP.LE.NSPL+NSPS+NSPG ) ISKIP = 0
                  ENDIF
!
!---              Read species concentration in well fluid, depending
!                 on species type and active fluids  ---
!
                  IF( ISKIP.EQ.0 ) THEN
                    VARB = 'Species Concentration in Well Fluid'
                    CALL RDDPR(ISTART,ICOMMA,CHDUM,
     &                VARSP_CW(NSPCW,NC,NCW))
                    CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                    WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &                UNTS(1:NCH),': ',VARSP_CW(NSPCW,NC,NCW)
                    INDX = 0
                    IUNMOL = 1
                    IUNM = -3
                    CALL RDUNIT(UNTS,VARSP_CW(NSPCW,NC,NCW),INDX)
                    WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &               VARSP_CW(NSPCW,NC,NCW),', kmol/m^3 (i.e., mol/L))'
!
!---                Convert species concentration to mol/m^3 fluid  ---
!
                    VARSP_CW(NSPCW,NC,NCW) = 
     &                1.D+3*VARSP_CW(NSPCW,NC,NCW)
!
!---                Identify associated conservation components in 
!                   coupled well  ---
!
                    DO NEQ = 1,NEQC
                      DO IX = 1,IEQ_C(1,NEQ)
                        NSPC = IEQ_C(IX+1,NEQ)
                        IF( NSP.EQ.NSPC ) ISOLC_CW(NEQ,NCW) = 1
                      ENDDO
                    ENDDO
!
!---                Identify associated kinetic components in 
!                   coupled well  ---
!
                    DO NEQ = 1,NEQK
                      DO IX = 1,IEQ_K(1,NEQ)
                        NSPK = IEQ_K(IX+1,NEQ)
                        IF( NSP.EQ.NSPK ) ISOLK_CW(NEQ,NCW) = 1
                      ENDDO
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF
#endif
!
!---        Read production well-bottom pressure  ---
!
            ELSEIF( IT_CW(NCW).EQ.-1 .OR. IT_CW(NCW).EQ.-2 ) THEN
!
!---          Read production well-bottom pressure  ---
!
              VARB = 'Production Well-Bottom Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(3,NC,NCW))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_CW(3,NC,NCW)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_CW(3,NC,NCW),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(3,NC,NCW),', Pa)'
#ifdef ecke
!
!---          Simulation with reactive species transport  ---
!
              IF( ISLC(40).EQ.1 ) THEN
!
!---            Identify all conservation components in 
!               production type coupled well  ---
!
                DO NEQ = 1,NEQC
                  ISOLC_CW(NEQ,NCW) = 1
                ENDDO
!
!---            Identify all kinetic components in 
!               production type coupled well  ---
!
                DO NEQ = 1,NEQK
                  ISOLK_CW(NEQ,NCW) = 1
                ENDDO
              ENDIF
#endif
            ENDIF
  200     CONTINUE
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
!     STOMP-EOR
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
      USE GLB_PAR
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
            RSDFX = 0.D+0
            RSDPX = ABS(DP_CWX)/(ABS(P_CW(2,NCW))+PATM)
            RSDX = MAX( RSDPX,RSDFX )
            RSDX = ABS( RS_CW(1,NCW)/(QM_CW(1,NCW)+QM_CW(3,NCW)) )
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
            RSDFX = 0.D+0
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
!     Written by M.D. White, PNNL, 16 May 2016 (First National
!     Challenge Cup in 1914, Brooklyn Field Club defeats Brooklyn
!     Celtic 2-1).
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
!---    Solute not found in well cycle to next well ---
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            IROW = petsc_offset(N)
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
!
!---        Solute produced via nonaqueous-liquid phase production  ---
!
            CNX = PORD(2,N)*SN(2,N)
            IF( CNX.GT.EPSL ) SORTX = SORTX + Q_CW(4,NWN)/CNX
          ENDIF
!
!---      Load Jacobian  ---
!
          IF( ILES.EQ.1 ) THEN
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + SORTX
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,SORTX,
     &                        ADD_VALUES,IERR)
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
            ITS_CWX = ITS_CW(NTP,NCW)
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
!---        Well state #1  ---
!
!           SL = 1.0
!           SG+SN = 0.0
!
!           Aqueous species only
!
            IF( ITS_CWX/100.EQ.1 ) THEN
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
!---        Well state #2  ---
!
!           SL < 1.0
!           SG+SN > 0.0
!
!           Aqueous, gas, and nonaqueous-liquid species
!
            ELSEIF( ITS_CWX/100.EQ.2 ) THEN
              DO NSPCW = 1,NSP_CW(NCW)
                NSP = ISPC_CW(NSPCW,NCW)
                IF( NSP.EQ.0 ) EXIT
                IF( NSP.LE.NSPL .OR. (NSP.GT.NSPL+NSPS .AND. 
     &            NSP.LE.NSPL+NSPS+NSPG) ) THEN
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
!---        Well state #3  ---
!
!           SL = 0.0
!           SG+SN = 1.0
!
!           Gas and nonaqueous-liquid species only
!      
            ELSEIF( ITS_CWX/100.EQ.3 ) THEN
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
            ITS_CWX = ITS_CW(NTP,NCW)
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
!---            Well state #1  ---
!
!               SL = 1.0
!               SG+SN = 0.0
!
!               Aqueous species only
!
                IF( ITS_CWX/100.EQ.1 ) THEN
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
!---            Well state #2  ---
!
!               SL < 1.0
!               SG+SN > 0.0
!
!               Aqueous, gas, and nonaqueous-liquid species
!
                ELSEIF( ITS_CWX/100.EQ.2 ) THEN
                  DO NSPCW = 1,NSP_CW(NCW)
                    NSP = ISPC_CW(NSPCW,NCW)
                    IF( NSP.EQ.0 ) EXIT
                    IF( NSP.LE.NSPL .OR. (NSP.GT.NSPL+NSPS .AND. 
     &                NSP.LE.NSPL+NSPS+NSPG) ) THEN
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
!---            Well state #3  ---
!
!               SL = 0.0
!               SG+SN = 1.0
!
!               Gas and nonaqueous-liquid species only
!          
                ELSEIF( ITS_CWX/100.EQ.3 ) THEN
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            IROW = petsc_offset(N)
#endif
          ENDIF
          SORTX = 0.D+0
!
!---      Injection well (volumetric fluxes are positive from well) ---
!
!         Q_CW(1,NWN) - total volumetric flux, m^3/s
!         Q_CW(2,NWN) - aqueous volumetric flux, m^3/s
!         Q_CW(3,NWN) - gas volumetric flux, m^3/s
!         Q_CW(4,NWN) - nonaqueous-liquid volumetric flux, m^3/s
!
          IF( IT_CW(NCW).GT.0 ) THEN
!
!---        Well state #1  ---
!
!           SL = 1.0
!           SG+SN = 0.0
!
!           Aqueous species only,
!           mol/s = m^3 fluid/s * mol/m^3 fluid
!
            IF( ITS_CWX/100.EQ.1 ) THEN
              BLU(MP) = BLU(MP) + Q_CW(2,NWN)*VARC_CWX
!
!---        Well state #2  ---
!
!           SL < 1.0
!           SG+SN > 0.0
!
!           Aqueous, gas, and nonaqueous-liquid species,
!           mol/s = m^3 fluid/s * mol/m^3 fluid
!
            ELSEIF( ITS_CWX/100.EQ.2 ) THEN
              BLU(MP) = BLU(MP) + Q_CW(1,NWN)*VARC_CWX
!
!---        Well state #3  ---
!
!           SL = 0.0
!           SG+SN = 1.0
!
!           Gas and nonaqueous-liquid species only,
!           mol/s = m^3 fluid/s * mol/m^3 fluid
!            
            ELSEIF( ITS_CWX/100.EQ.3 ) THEN
              BLU(MP) = BLU(MP) + (Q_CW(3,NWN)+Q_CW(4,NWN))*VARC_CWX
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
!
!---        Solute produced via nonaqueous-liquid phase production  ---
!
            CNX = PORD(2,N)*SN(2,N)
            IF( CNX.GT.EPSL ) SORTX = SORTX + Q_CW(4,NWN)/CNX
          ENDIF
!
!---      Load Jacobian  ---
!
          IF( ILES.EQ.1 ) THEN
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + SORTX
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,SORTX,
     &                        ADD_VALUES,IERR)
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
            ITS_CWX = ITS_CW(NTP,NCW)
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
!---        Well state #1  ---
!
!           SL = 1.0
!           SG+SN = 0.0
!
!           Aqueous species only
!
            IF( ITS_CWX/100.EQ.1 ) THEN
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
!---        Well state #2  ---
!
!           SL < 1.0
!           SG+SN > 0.0
!
!           Aqueous, gas, and nonaqueous-liquid species
!
            ELSEIF( ITS_CWX/100.EQ.2 ) THEN
              DO NSPCW = 1,NSP_CW(NCW)
                NSP = ISPC_CW(NSPCW,NCW)
                IF( NSP.EQ.0 ) EXIT
                IF( NSP.LE.NSPL .OR. (NSP.GT.NSPL+NSPS .AND. 
     &            NSP.LE.NSPL+NSPS+NSPG) ) THEN
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
!---        Well state #3  ---
!
!           SL = 0.0
!           SG+SN = 1.0
!
!           Gas and nonaqueous-liquid species only
!      
            ELSEIF( ITS_CWX/100.EQ.3 ) THEN
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
            ITS_CWX = ITS_CW(NTP,NCW)
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
!---            Well state #1  ---
!
!               SL = 1.0
!               SG+SN = 0.0
!
!               Aqueous species only
!
                IF( ITS_CWX/100.EQ.1 ) THEN
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
!---            Well state #2  ---
!
!               SL < 1.0
!               SG+SN > 0.0
!
!               Aqueous, gas, and nonaqueous-liquid species
!
                ELSEIF( ITS_CWX/100.EQ.2 ) THEN
                  DO NSPCW = 1,NSP_CW(NCW)
                    NSP = ISPC_CW(NSPCW,NCW)
                    IF( NSP.EQ.0 ) EXIT
                    IF( NSP.LE.NSPL .OR. (NSP.GT.NSPL+NSPS .AND. 
     &                NSP.LE.NSPL+NSPS+NSPG) ) THEN
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
!---            Well state #3  ---
!
!               SL = 0.0
!               SG+SN = 1.0
!
!               Gas and nonaqueous-liquid species only
!          
                ELSEIF( ITS_CWX/100.EQ.3 ) THEN
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            IROW = petsc_offset(N)
#endif
          ENDIF
          SORTX = 0.D+0
!
!---      Injection well (volumetric fluxes are positive from well) ---
!
!         Q_CW(1,NWN) - total volumetric flux, m^3/s
!         Q_CW(2,NWN) - aqueous volumetric flux, m^3/s
!         Q_CW(3,NWN) - gas volumetric flux, m^3/s
!         Q_CW(4,NWN) - nonaqueous-liquid volumetric flux, m^3/s
!
          IF( IT_CW(NCW).GT.0 ) THEN
!
!---        Well state #1  ---
!
!           SL = 1.0
!           SG+SN = 0.0
!
!           Aqueous species only,
!           mol/s = m^3 fluid/s * mol/m^3 fluid
!
            IF( ITS_CWX/100.EQ.1 ) THEN
              BLU(MP) = BLU(MP) + Q_CW(2,NWN)*VARC_CWX
!
!---        Well state #2  ---
!
!           SL < 1.0
!           SG+SN > 0.0
!
!           Aqueous, gas, and nonaqueous-liquid species,
!           mol/s = m^3 fluid/s * mol/m^3 fluid
!
            ELSEIF( ITS_CWX/100.EQ.2 ) THEN
              BLU(MP) = BLU(MP) + Q_CW(1,NWN)*VARC_CWX
!
!---        Well state #3  ---
!
!           SL = 0.0
!           SG+SN = 1.0
!
!           Gas and nonaqueous-liquid species only,
!           mol/s = m^3 fluid/s * mol/m^3 fluid
!            
            ELSEIF( ITS_CWX/100.EQ.3 ) THEN
              BLU(MP) = BLU(MP) + (Q_CW(3,NWN)+Q_CW(4,NWN))*VARC_CWX
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
!
!---        Solute produced via nonaqueous-liquid phase production  ---
!
            CNX = PORD(2,N)*SN(2,N)
            IF( CNX.GT.EPSL ) SORTX = SORTX + Q_CW(4,NWN)/CNX
          ENDIF
!
!---      Load Jacobian  ---
!
          IF( ILES.EQ.1 ) THEN
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + SORTX
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,SORTX,
     &                        ADD_VALUES,IERR)
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
!     STOMP-EOR
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
      USE GLB_PAR
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
        DPX = MAX( 1.D-1*(P_CW(2,NCW)+PATM)/SQRT(REAL(NITER)),5.D+5 )
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
      SUBROUTINE WRCOUP_WELL
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
!     Write well.dat file.
!     
!     FF_CW(1,LN_CW) - x-direction well fraction factor for well
!     FF_CW(2,LN_CW) - y-direction well fraction factor for well
!     FF_CW(3,LN_CW) - z-direction well fraction factor for well
!     ICC_CW(LN_CW) - cyclic well time index for well
!     IM_CW(LN_CW) - number of time points for well
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
      USE GLB_PAR
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
      SAVE FORM1
      DATA FORM1 /'(I1)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRCOUP_WELL'
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
!---  Close 'well.dat' file  ---
!
      CLOSE(UNIT=26)
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRCOUP_WELL group  ---
!
      RETURN
      END


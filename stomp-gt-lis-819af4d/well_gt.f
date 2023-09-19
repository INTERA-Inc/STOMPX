!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_GRLP_WELL
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
!     STOMP-GT
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
!      REAL*8 PAX(3),PBX(3),PCX(3),PBCX(3)
      REAL*8 AJ(3,3),BJ(3)
      INTEGER MSX(4,6),IJ(3)
      INTEGER N1X(4),N2X(4)
      CHARACTER*9 FORM1
      CHARACTER*8 FORM2
      CHARACTER*23 FORM5
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
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
      SUB_LOG(ISUB_LOG) = '/CHK_GRLP_WELL'
      EPSLX = 1.D-12
!
!---  Loop over ground-loop wells ---
!
      DO 600 NCW = 1,N_GLW
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
          DO 480 N = NFLD,1,-1
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
!!
!!---          Loop over node surfaces,
!!             (bottom,south,west,east,north,top)  ---
!!
!              ICWX = 0
!              DO 180 NS = 1,6
!!
!!---            Define the five surface points, four corners
!!               and one centroid---
!!
!                DO 110 NP = 1,4
!                  MX = MSX(NP,NS)
!!
!!---              Cylindrical coordinates  ---
!!
!                  IF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
!                    XPX(NP) = XE(MX,N)*COS(YE(MX,N))
!                    YPX(NP) = XE(MX,N)*SIN(YE(MX,N))
!                    ZPX(NP) = ZE(MX,N)
!                  ELSE
!                    XPX(NP) = XE(MX,N)
!                    YPX(NP) = YE(MX,N)
!                    ZPX(NP) = ZE(MX,N)
!                  ENDIF
!  110           CONTINUE
!                NP = 4
!                CALL PGCNTRD( NP,XPX(1),YPX(1),ZPX(1),
!     &            XPX(5),YPX(5),ZPX(5) )
!!
!!---            Loop over the four triangular planes on the 
!!               surface face  ---
!!
!                DO 130 NT = 1,4
!!
!!---              Built vectors between transition point
!!                 and triangular plane points  ---
!!
!                  PAX(1) = XPX(5)-XTP_CW(NPT,NICW)
!                  PAX(2) = YPX(5)-YTP_CW(NPT,NICW)
!                  PAX(3) = ZPX(5)-ZTP_CW(NPT,NICW)
!                  PBX(1) = XPX(N1X(NT))-XTP_CW(NPT,NICW)
!                  PBX(2) = YPX(N1X(NT))-YTP_CW(NPT,NICW)
!                  PBX(3) = ZPX(N1X(NT))-ZTP_CW(NPT,NICW)
!                  PCX(1) = XPX(N2X(NT))-XTP_CW(NPT,NICW)
!                  PCX(2) = YPX(N2X(NT))-YTP_CW(NPT,NICW)
!                  PCX(3) = ZPX(N2X(NT))-ZTP_CW(NPT,NICW)
!                  CALL VCROSSP( PBX,PCX,PBCX )
!                  SX = VDOTP( PAX,PBCX )
!!
!!---              Clockwise rotation  ---
!!
!                  IF( SX.GT.EPSL ) THEN
!!
!!---                Opposing rotations found, point outside hexahedron  ---
!!
!                    IF( ICWX.EQ.-1 ) THEN
!                      GOTO 190
!!
!!---                Similar rotations found, continue searching  ---
!!
!                    ELSE
!                      ICWX = 1
!                    ENDIF
!!
!!---              Counterclockwise rotation  ---
!!
!                  ELSEIF( SX.LT.-EPSL ) THEN
!!
!!---                Opposing rotations found, point outside hexahedron  ---
!!
!                    IF( ICWX.EQ.1 ) THEN
!                      GOTO 190
!!
!!---                Similar rotations found, continue searching  ---
!!
!                    ELSE
!                      ICWX = -1
!                    ENDIF
!                  ENDIF
!  130           CONTINUE
!  180         CONTINUE
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
                    CHMSGX = 'Three Distinct Coupled Well Points'
     &                //  ' at Node'
                    IMSGSX = N
                    NMSGX = 0
                    SUBLOGX = 'CHK_GRLP_WELL'
                    RLMSGX = 0.D+0
                    CALL WRMSGX(RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX)
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
                  CHMSGX = 'Both Well Trajectory Points within '
     &              //  'Surface Plane at Node'
                  IMSGX = N
                  NMSGX = 0
                  SUBLOGX = 'CHK_GRLP_WELL'
                  RLMSGX = 0.D+0
                  CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
  358           CONTINUE
  360         CONTINUE
  362         CONTINUE
!
!---          Loop over well nodes in current well  ---
!
              DO 370 NWN = ID_CW(3,NCW),NWN_CW
!!
!!---          Loop over all well nodes  ---
!!
!              DO 370 NWN = 1,NWN_CW
!
!---            Well node previously counted  ---
!
                IF( IWN_CW(NWN).EQ.N ) THEN
                  GOTO 372
                ENDIF
  370         CONTINUE
!
!---          Increment the number of field nodes with
!             ground-loop-well nodes  ---
!
              NWF_CW = NWF_CW + 1
              IF( NWF_CW.GT.LWF_CW ) THEN
                INDX = 5
                CHMSGX = 'Number of Field Nodes with Coupled Wells ' //
     &            '> Parameter LWF_CW'
                IMSGX = 0
                NMSGX = 0
                SUBLOGX = 'CHK_GRLP_WELL'
                RLMSGX = 0.D+0
                CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
              ENDIF
  372         CONTINUE
!
!---          Define a new well node and set the trajectory
!             points  ---
!
              NWN_CW = NWN_CW + 1
              IF( NWN_CW.GT.LWN_CW ) THEN
               INDX = 5
               CHMSGX = 'Number of Coupled-Well Nodes ' //
     &           '> Parameter LWN_CW'
                IMSGX = 0
                NMSGX = 0
                SUBLOGX = 'CHK_GRLP_WELL'
                RLMSGX = 0.D+0
                CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
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
              CALL PROJ_GRLP_WELL( N,NWN )
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
!---  Load array for field nodes with ground-loop-well nodes  ---
!
      NC = 0
      DO 640 NCW = 1,N_GLW
        NUK_CW = (ISVC*(ABS(ID_CW(6,NCW)-ID_CW(5,NCW))+1))
        IF( NUK_CW.GT.LUK_CW ) THEN
          INDX = 5
          CHMSGX = 'Number of Well Equation Unknowns ' //
     &      '> Parameter LUK_CW'
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'CHK_GRLP_WELL'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
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
!
!---    Set well head temperature equal to node temperature  ---
!
        N = IWF_CW(ID_CW(3,NCW))
        P_CW(2,NCW) = T(2,N)
  640 CONTINUE
!
!---  Record coupled-well nodes to output file ---
!
      WRITE(IWR,'(//,A,/)') ' --- Ground-Loop-Well Node Record  ---'
      DO 720 NCW = 1,N_GLW
        NCH = INDEX( WNM_CW(NCW),'  ' )-1
        WRITE(IWR,'(2X,A,I4,2A)') 'Nodes With Ground-Loop-Well ' //
     &    'Number: ',NCW,' Ground-Loop-Well Name: ',WNM_CW(NCW)(1:NCH)
        DO 710 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
          NX = IWF_CW(NWF)
          IX = ID(NX)
          JX = JD(NX)
          KX = KD(NX)
          WRITE(IWR,FORM5) ' (',IX,',',JX,',',KX,':',NX,')'
  710   CONTINUE
  720 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_GRLP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE GRLP_WELL( NCW )
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
!     STOMP-GT
!
!     Ground-loop well model
!
!     Flux of energy from well nodes to field nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 January 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE JACOB
      USE GRID
      USE FDVT
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
      REAL*8 VAR_CWX(6+LNGC)
      REAL*8 C1FX(LWN_CW)
      REAL*8 AJ(2*LWN_CW,2*LWN_CW),BJ(2*LWN_CW)
      INTEGER IJ(2*LWN_CW)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/GRLP_WELL'
!
!---  Loop over ground-loop-well nodes  ---
!
      DO 30 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---    Loop over increment indices  ---
!
        DO 20 M = 1,ISVC+2
          FXE_CW(M,NWN) = 0.D+0
   20   CONTINUE
   30  CONTINUE
!
!---  Well time interval ---
!
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        IF( ICC_CW(NCW).EQ.1 ) TMZ = MOD( TM,VAR_CW(1,IM_CW(NCW),NCW) )
!
!---    Coupled well is inactive  ---
!
        IF( TMZ.LE.VAR_CW(1,1,NCW) ) THEN
          GOTO 500
        ENDIF
        IF( IM_CW(NCW).EQ.1 ) THEN
          DO 80 N = 2,7
            VAR_CWX(N) = VAR_CW(N,1,NCW)
   80     CONTINUE
        ELSE
          DO 100 M = 2,IM_CW(NCW)
            IF( TMZ.LE.VAR_CW(1,M,NCW) ) THEN
              TD_CW = VAR_CW(1,M,NCW)-VAR_CW(1,M-1,NCW)
              DT_CW = MIN( VAR_CW(1,M,NCW)-TMZ,DT )
              TF_CW = (TMZ-VAR_CW(1,M-1,NCW))/TD_CW
              DO 90 N = 2,7
                VAR_CWX(N) = VAR_CW(N,M-1,NCW) + 
     &            TF_CW*(VAR_CW(N,M,NCW)-VAR_CW(N,M-1,NCW))
   90         CONTINUE
              GOTO 110
            ENDIF
  100     CONTINUE
!
!---      Coupled well is inactive  ---
!
          GOTO 500
        ENDIF
  110   CONTINUE
!
!---  Field node index ---
!
      N = IWN_CW(ID_CW(3,NCW))
!
!---  Store top of ground-loop-well location in previous
!     ground-loop-well node location  ---
!
      XPX(1) = XTP_CW(1,ID_CW(1,NCW))
      YPX(1) = YTP_CW(1,ID_CW(1,NCW))
      ZPX(1) = ZTP_CW(1,ID_CW(1,NCW))
!
!---  Index offset between annulus and inner temperatures  ---
!
      I12X = ID_CW(4,NCW) - ID_CW(3,NCW) + 1
!
!---  Zero coefficients and problem vectors  ---
!
      DO 152 J = 1,2*I12X
        DO 150 I = 1,2*I12X
          AJ(I,J) = 0.D+0
  150   CONTINUE
        BJ(J) = 0.D+0
  152 CONTINUE
!
!---  Loop over the well nodes in the ground-loop well ---
!
      DO 200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
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
        INDX = 1
        THKEX = THKE_L( IZN,SL(1,N),THKL(1,N),PORD(1,N),PORT(1,N),INDX )
        INDX = 2
        THKEY = THKE_L( IZN,SL(1,N),THKL(1,N),PORD(1,N),PORT(1,N),INDX )
        INDX = 3
        THKEZ = THKE_L( IZN,SL(1,N),THKL(1,N),PORD(1,N),PORT(1,N),INDX )
        R1X = MAX( PAR_CW(1,INVX),1.D-20 )
!
!---    Cylindrical coordinates with azimuthal symmetry,
!       centrally located wells  ---
!
        IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
          ROZ = RP(I)
          R1X = MIN( R1X,9.999D-1*ROZ )
          WI_CWX = 2.D+0*GPI*THKEX*ZLX/(LOG(ROZ/R1X)+PAR_CW(1,INVX))
        ELSE
          THKEYZ = SQRT(THKEY/THKEZ)
          THKEZY = SQRT(THKEZ/THKEY)
          DXGFX = DXGF(N)/FF_CW(1,NCW)
          DYGFX = DYGF(N)*RP(I)/FF_CW(2,NCW)
          DZGFX = DZGF(N)/FF_CW(3,NCW)
          ROX = 2.8D-1*SQRT(THKEYZ*(DZGFX**2) + THKEZY*(DYGFX**2))
     &    /(SQRT(THKEYZ)+SQRT(THKEZY))
          THKEZX = SQRT(THKEZ/THKEX)
          THKEXZ = SQRT(THKEX/THKEZ)
          ROY = 2.8D-1*SQRT(THKEZX*(DXGFX**2) + THKEXZ*(DZGFX**2))
     &      /(SQRT(THKEZX)+SQRT(THKEXZ))
          THKEYX = SQRT(THKEY/THKEX)
          THKEXY = SQRT(THKEX/THKEY)
          ROZ = 2.8D-1*SQRT(THKEYX*(DXGFX**2) + THKEXY*(DYGFX**2))
     &      /(SQRT(THKEYX)+SQRT(THKEXY))
!
!---      Well index components  ---
!
          WIX = 2.D+0*GPI*SQRT(THKEY*THKEZ)*XLX/LOG(ROX/R1X)
          WIY = 2.D+0*GPI*SQRT(THKEX*THKEZ)*YLX/LOG(ROY/R1X)
          WIZ = 2.D+0*GPI*SQRT(THKEX*THKEY)*ZLX/LOG(ROZ/R1X)
          WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
        ENDIF
!
!---    Store current ground-loop-well node location in previous
!       ground-loop-well node location  ---
!
        XPX(1) = XPX(2)
        YPX(1) = YPX(2)
        ZPX(1) = ZPX(2)     
!
!---    Single-pipe model  ---
!
        IF( IT_CW(NCW).EQ.1 ) THEN
!
!---      Load well coefficient matrix and problem vector  ---
!
          DLX = SQRT((XLX**2) + (YLX**2) + (ZLX**2))
          DLX = DLX*PAR_CW(3,INVX)
          R2X = PAR_CW(2,INVX)
          THKPX = PAR_CW(4,INVX)
          H1FX = PAR_CW(5,INVX)
!
!---      Load well coefficient matrix and problem vector  ---
!
          DLX = SQRT((XLX**2) + (YLX**2) + (ZLX**2))
          CPMX = VAR_CWX(3)*VAR_CWX(4)*VAR_CWX(5)
          A1FX = 2.D+0*R1X*GPI*DLX
          RPX = LOG( R1X/R2X )/(2.D+0*GPI*DLX*THKPX)
          C1FX(NWN) = 1.D+0/(1.D+0/(A1FX*H1FX) + RPX + 1.D+0/WI_CWX)
          CMX = ABS(CPMX)
!
!---      First well node  ---
!
          IF( NWN.EQ.ID_CW(3,NCW) ) THEN
!
!---        Row index for pipe temperature at well node  ---
!
            NRX = NWN - ID_CW(3,NCW) + 1
!
!---        Coefficient for pipe temperature at well node  ---
!
            NCX = NRX
            AJ(NRX,NCX) = C1FX(NWN) + CMX
!
!---        Coefficient for pipe temperature at next well node,
!           only if next node exists and flow is negative  ---
!
            IF( I12X.GT.1 .AND. CPMX.LT.0.D+0 ) THEN
              NCX = NRX + 1
              AJ(NRX,NCX) = -CMX
            ENDIF
!
!---        Problem vector for pipe temperature at well node if
!           flow is positive  ---
!
            IF( CPMX.GE.0.D+0 ) THEN
              BJ(NRX) = C1FX(NWN)*T(1,N) + CMX*VAR_CWX(2)
            ELSE
              BJ(NRX) = C1FX(NWN)*T(1,N)
            ENDIF
!
!---      Remaining well nodes  ---
!
          ELSE
!
!---        Row index for pipe temperature at well node  ---
!
            NRX = NWN - ID_CW(3,NCW) + 1
!
!---        Coefficient for pipe temperature at well node  ---
!
            NCX = NRX
            AJ(NRX,NCX) = C1FX(NWN) + CMX
!
!---        Coefficient for pipe temperature at prior well node if
!           flow is positive  ---
!
            IF( CPMX.GE.0.D+0 ) THEN
              NCX = NRX - 1
              AJ(NRX,NCX) = -CMX
!
!---        Coefficient for pipe temperature at next well node if
!           flow is negative  ---
!
            ELSE
              NCX = NRX + 1
              AJ(NRX,NCX) = -CMX
            ENDIF
!
!---        Problem vector for pipe temperature at well node  ---
!
            BJ(NRX) = C1FX(NWN)*T(1,N)
          ENDIF
!
!---    Double-pipe model  ---
!
        ELSEIF( IT_CW(NCW).EQ.2 ) THEN
!
!---      Load well coefficient matrix and problem vector  ---
!
          DLX = SQRT((XLX**2) + (YLX**2) + (ZLX**2))
          H1FX = PAR_CW(2,INVX)
          H12X = PAR_CW(4,INVX)
          CPMX = VAR_CWX(3)*VAR_CWX(4)*VAR_CWX(5)
          A1FX = 2.D+0*R1X*GPI*DLX
          R2X = MAX( PAR_CW(3,INVX),1.D-20 )
          A12X = 2.D+0*R2X*GPI*DLX
          C1FX(NWN) = 1.D+0/(1.D+0/(A1FX*H1FX) + 1.D+0/WI_CWX)
          CMX = ABS(CPMX)
!
!---      First well node  ---
!
          IF( NWN.EQ.ID_CW(3,NCW) ) THEN
!
!---        Row index for annulus temperature at well node  ---
!
            NRX = NWN - ID_CW(3,NCW) + 1
!
!---        Coefficient for annulus temperature at well node  ---
!
            NCX = NRX
            AJ(NRX,NCX) = C1FX(NWN) + A12X*H12X + CMX
!
!---        Coefficient for annulus temperature at next well node,
!           only if next node exists and flow is negative  ---
!
            IF( I12X.GT.1 .AND. CPMX.LT.0.D+0 ) THEN
              NCX = NRX + 1
              AJ(NRX,NCX) = -CMX
            ENDIF
!
!---        Coefficient for inner temperature at well node  ---
!
            NCX = NRX + I12X
            AJ(NRX,NCX) = -A12X*H12X
!
!---        Problem vector for annulus temperature at well node if
!           flow is positive  ---
!
            IF( CPMX.GE.0.D+0 ) THEN
              BJ(NRX) = C1FX(NWN)*T(1,N) + CMX*VAR_CWX(2)
            ELSE
              BJ(NRX) = C1FX(NWN)*T(1,N)
            ENDIF
!
!---        Row index for inner temperature at well node  ---
!
            NRX = NRX + I12X
!
!---        Coefficient for inner temperature at well node  ---
!
            NCX = NRX
            AJ(NRX,NCX) = A12X*H12X + CMX
!
!---        Coefficient for inner temperature at next well node,
!           only if next node exists and flow is positive  ---
!
            IF( I12X.GT.1 .AND. CPMX.GE.0.D+0 ) THEN
              NCX = NRX + 1
              AJ(NRX,NCX) = -CMX
            ENDIF
!
!---        Coefficient for annulus temperature at well node  ---
!
            NCX = NRX - I12X
            AJ(NRX,NCX) = -A12X*H12X
!
!---        Problem vector for inner temperature at well node if
!           flow is negative  ---
!
            IF( CPMX.LT.0.D+0 ) THEN
              BJ(NRX) = CMX*VAR_CWX(2)
            ENDIF     
!
!---      Last well node  ---
!
          ELSEIF( NWN.EQ.ID_CW(4,NCW) ) THEN
!
!---        Row index for annulus temperature at well node  ---
!
            NRX = NWN - ID_CW(3,NCW) + 1
!
!---        Coefficient for annulus temperature at well node  ---
!
            NCX = NRX
            AJ(NRX,NCX) = C1FX(NWN) + A12X*H12X + CMX
!
!---        Coefficient for annulus temperature at prior well node if
!           flow is positive  ---
!
            IF( CPMX.GE.0.D+0 ) THEN
              NCX = NRX - 1
              AJ(NRX,NCX) = -CMX
            ENDIF
!
!---        Coefficient for inner temperature at well node if flow
!           is positive  ---
!
            IF( CPMX.GE.0.D+0 ) THEN
              NCX = NRX + I12X
              AJ(NRX,NCX) = -A12X*H12X
!
!---        Coefficient for inner temperature at well node if flow
!           is negative, inner feeds annulus at last node  ---
!
            ELSE
              NCX = NRX + I12X
              AJ(NRX,NCX) = -(A12X*H12X + CMX)
            ENDIF
!
!---        Problem vector for annulus temperature at well node  ---
!
            BJ(NRX) = C1FX(NWN)*T(1,N)
!
!---        Row index for inner temperature at well node  ---
!
            NRX = NRX + I12X
!
!---        Coefficient for inner temperature at well node  ---
!
            NCX = NRX
            AJ(NRX,NCX) = A12X*H12X + CMX
!
!---        Coefficient for inner temperature at prior well node if
!           flow is negative  ---
!
            IF( CPMX.LT.0.D+0 ) THEN
              NCX = NRX - 1
              AJ(NRX,NCX) = -CMX
            ENDIF
!
!---        Coefficient for annulus temperature at well node if flow
!           is positive, annulus feeds inner at last node  ---
!
            IF( CPMX.GE.0.D+0 ) THEN
              NCX = NRX - I12X
              AJ(NRX,NCX) = -(A12X*H12X + CMX)
!
!---        Coefficient for annulus temperature at well node if flow
!           is negative  ---
!
            ELSE
              NCX = NRX - I12X
              AJ(NRX,NCX) = -A12X*H12X
            ENDIF
!
!---        Problem vector for inner temperature at well node  ---
!
            BJ(NRX) = 0.D+0            
!
!---      Intermediate well nodes  ---
!
          ELSE
!
!---        Row index for annulus temperature at well node  ---
!
            NRX = NWN - ID_CW(3,NCW) + 1
!
!---        Coefficient for annulus temperature at well node  ---
!
            NCX = NRX
            AJ(NRX,NCX) = C1FX(NWN) + A12X*H12X + CMX
!
!---        Coefficient for annulus temperature at prior well node if
!           flow is positive  ---
!
            IF( CPMX.GE.0.D+0 ) THEN
              NCX = NRX - 1
              AJ(NRX,NCX) = -CMX
!
!---        Coefficient for annulus temperature at next well node if
!           flow is negative  ---
!
            ELSE
              NCX = NRX + 1
              AJ(NRX,NCX) = -CMX
            ENDIF
!
!---        Coefficient for inner temperature at well node  ---
!
            NCX = NRX + I12X
            AJ(NRX,NCX) = -A12X*H12X
!
!---        Problem vector for annulus temperature at well node  ---
!
            BJ(NRX) = C1FX(NWN)*T(1,N)
!
!---        Row index for inner temperature at well node  ---
!
            NRX = NRX + I12X
!
!---        Coefficient for inner temperature at well node  ---
!
            NCX = NRX
            AJ(NRX,NCX) = A12X*H12X + CMX
!
!---        Coefficient for inner temperature at next well node, if
!           flow is positive  ---
!
            IF( CPMX.GE.0.D+0 ) THEN
              NCX = NRX + 1
              AJ(NRX,NCX) = -CMX
!
!---        Coefficient for inner temperature at prior well node, if
!           flow is negative  ---
!
            ELSE
              NCX = NRX - 1
              AJ(NRX,NCX) = -CMX
            ENDIF
!
!---        Coefficient for annulus temperature at well node  ---
!
            NCX = NRX - I12X
            AJ(NRX,NCX) = -A12X*H12X
!
!---        Problem vector for inner temperature at well node  ---
!
            BJ(NRX) = 0.D+0            
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  Single-pipe model  ---
!
      IF( IT_CW(NCW).EQ.1 ) THEN
        JP = I12X
        KP = 2*LWN_CW
!
!---    Solve linear system  ---
!
        CALL LUDCMP( AJ,JP,KP,IJ,DJ )
        CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---    Loop over well nodes  ---
!
        DO 220 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          NRX = NWN - ID_CW(3,NCW) + 1
          N = IWN_CW(NWN)
!
!---      Load pipe temperatures  ---
!
          T_CW(1,NWN) = BJ(NRX)
!
!---      Energy fluxes from well node annulus to field nodes  ---
!
          DO 210 M = 2,ISVC+2
            FXE_CW(M,NWN) = C1FX(NWN)*(T_CW(1,NWN)-T(1,N))
            SRCT(M,N) = SRCT(M,N) + C1FX(NWN)*(T_CW(1,NWN)-T(1,N))
  210     CONTINUE
  220   CONTINUE
!
!---  Double-pipe model  ---
!
      ELSEIF( IT_CW(NCW).EQ.2 ) THEN
        JP = 2*I12X
        KP = 2*LWN_CW
!
!---    Solve linear system  ---
!
        CALL LUDCMP( AJ,JP,KP,IJ,DJ )
        CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---    Loop over well nodes  ---
!
        DO 240 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          NRX = NWN - ID_CW(3,NCW) + 1
          N = IWN_CW(NWN)
!
!---      Load annulus and inner well temperatures  ---
!
          T_CW(1,NWN) = BJ(NRX)
          T_CW(2,NWN) = BJ(NRX+I12X)
!
!---      Energy fluxes from well node annulus to field nodes  ---
!
          DO 230 M = 2,ISVC+2
            FXE_CW(M,NWN) = C1FX(NWN)*(T_CW(1,NWN)-T(1,N))
            SRCT(M,N) = SRCT(M,N) + C1FX(NWN)*(T_CW(1,NWN)-T(1,N))
  230     CONTINUE
  240   CONTINUE
      ENDIF
!
!---  For single-pass well load pipe outlet temperature for output  ---
!
      IF( IT_CW(NCW).EQ.1 ) THEN
       P_CW(2,NCW) = T_CW(1,ID_CW(4,NCW))
!
!---  For double-pass well load top of well temperature for output  ---
!
      ELSEIF( IT_CW(NCW).EQ.2 ) THEN
        IF( CPMX.GE.0.D+0 ) THEN
          P_CW(2,NCW) = T_CW(2,ID_CW(3,NCW))
        ELSE
          P_CW(2,NCW) = T_CW(1,ID_CW(3,NCW))
        ENDIF
      ENDIF
  500 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of GRLP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROJ_GRLP_WELL( N,NWN )
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
!     STOMP-GT
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
      SUB_LOG(ISUB_LOG) = '/PROJ_GRLP_WELL'
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
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROJ_GRLP_WELL group  ---
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
!---  End of PROJ_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RD_GRLP_WELL
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
!     STOMP-GT
!
!     Reads the ground-loop well card.
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
!     1 - single pipe, single pass
!     2 - double pipe, douple pass, annular and inner flow channels
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
!     N_GLW - number of ground-loop wells
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
!     N_GLW - number of ground-loop wells
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
      CHARACTER*64 ADUM,UNTS
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
      SUB_LOG(ISUB_LOG) = '/RD_GRLP_WELL'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Ground-Loop Well Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of wells  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Ground-Loop Wells'
      CALL RDINT(ISTART,ICOMMA,CHDUM,N_GLW)
!
!---  Check number of wells parameter  ---
!
      IF( N_GLW.GT.LN_CW ) THEN
        INDX = 5
        CHMSG = 'Number of Ground-Loop Wells > LN_CW'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Loop over number of ground-loop wells  ---
!
      NIT_CW = 0
      DO 400 NCW = 1,N_GLW
!
!---    Read well type  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Ground-Loop Well Type'
        CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
!
!---    Single-pass, single-pipe well  ---
!
        IF( INDEX(ADUM(1:),'single').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Single-Pass Ground-Loop Well'
          IT_CW(NCW) = 1
          WNM_CW(NCW) = 'SP'
          ICX = ICOUNT(NCW)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(WNM_CW(NCW)(3:3+ICX-1),FORM1) NCW
!
!---    Double-pass, double-pipe well  ---
!
        ELSEIF( INDEX(ADUM(1:),'double').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Double-Pass Double-Pipe Well'
          WRITE(IWR,'(2X,A)') 'Annulus and Inner Pipe Design'
          IT_CW(NCW) = 2
          WNM_CW(NCW) = 'DP'
          ICX = ICOUNT(NCW)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(WNM_CW(NCW)(3:3+ICX-1),FORM1) NCW
!
!---    Unrecognized well type  ---
!
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Well Type: ' // ADUM(1:NCH)
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
!---      Single-pipe, single-pass system  ---
!
          IF( IT_CW(NCW).EQ.1 ) THEN
!
!---        Read pipe outside radius  ---
!
            VARB = 'Interval Outside Pipe Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(1,NICW))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_CW(1,NICW)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_CW(1,NICW),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_CW(1,NICW),', m)'
!
!---        Read pipe inside radius  ---
!
            VARB = 'Interval Inside Pipe Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(2,NICW))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_CW(2,NICW)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_CW(2,NICW),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_CW(2,NICW),', m)'
!
!---        Read pipe length per trajectory length  ---
!
            VARB = 'Interval Pipe Length per Trajectory Length'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(3,NICW))
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),': ',
     &        PAR_CW(3,NICW)
!
!---        Read pipe thermal conductivity ---
!
            VARB = 'Interval Pipe Thermal Conductivity'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(4,NICW))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_CW(4,NICW)
            INDX = 0
            IUNKG = 1
            IUNM = 1
            IUNS = -3
            IUNK = -1
            CALL RDUNIT(UNTS,PAR_CW(4,NICW),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_CW(4,NICW),', W/m K)'
!
!---        Read fluid-pipe heat transfer coefficient  ---
!
            VARB = 'Interval Fluid-Pipe Heat Transfer Coefficient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(5,NICW))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_CW(5,NICW)
            INDX = 0
            IUNKG = 1
            IUNS = -3
            IUNK = -1
            CALL RDUNIT(UNTS,PAR_CW(5,NICW),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_CW(5,NICW),', W/m^2 K)'
!
!---      Double-pipe well  ---
!
          ELSEIF( IT_CW(NCW).EQ.2 ) THEN
!
!---        Read well-bore radius  ---
!
            VARB = 'Interval Well-Bore Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(1,NICW))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_CW(1,NICW)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_CW(1,NICW),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_CW(1,NICW),', m)'
!
!---        Read fluid-casing-field overall heat transfer 
!           coefficient  ---
!
            VARB = 'Interval Outer-Casing Overall Heat Transfer Coeff.'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(2,NICW))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_CW(2,NICW)
            INDX = 0
            IUNKG = 1
            IUNS = -3
            IUNK = -1
            CALL RDUNIT(UNTS,PAR_CW(2,NICW),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_CW(2,NICW),', W/m^2 K)'
!
!---        Read inner-pipe radius  ---
!
            VARB = 'Interval Inner-Pipe Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(3,NICW))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_CW(3,NICW)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_CW(3,NICW),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_CW(3,NICW),', m)'
!
!---        Read fluid-inner pipe-fluid overall heat transfer 
!           coefficient  ---
!
            VARB = 'Interval Overall Heat Transfer Coeff.'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_CW(4,NICW))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_CW(4,NICW)
            INDX = 0
            IUNKG = 1
            IUNS = -3
            IUNK = -1
            CALL RDUNIT(UNTS,PAR_CW(4,NICW),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_CW(4,NICW),', W/m^2 K)'
          ENDIF
  100   CONTINUE
!
!---    Read number of well time points  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Number of Well Time Points'
        CALL RDINT(ISTART,ICOMMA,CHDUM,IM_CW(NCW))
!
!---    Check number of well time points parameter  ---
!
        IF( IM_CW(NCW).GT.LWT_CW ) THEN
          INDX = 5
          CHMSG = 'Number of Well Time Points > LWT_CW'
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
!---    Loop over the number of well times  ---
!
        DO 200 M = 1,IM_CW(NCW)
!
!---      Read new line  ---
!
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
!
!---      Read well time  ---
!
          VARB = 'Well Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(1,M,NCW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',VAR_CW(1,M,NCW)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,VAR_CW(1,M,NCW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(1,M,NCW),', s)'
          IF( M.GT.1 ) THEN
            IF( VAR_CW(1,M,NCW).LT.VAR_CW(1,M-1,NCW) ) THEN
              INDX = 7
              CHMSG = 'Nonascending Well Times: Well'
              IMSG = NCW
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
!
!---      Well inlet temperature  ---
!
          VARB = 'Well Inlet Temperature'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(2,M,NCW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',VAR_CW(2,M,NCW)
          IUNK = 1
          CALL RDUNIT(UNTS,VAR_CW(2,M,NCW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(2,M,NCW),
     &      ', C)'
!
!---      Well volumetric rate  ---
!
          VARB = 'Well Volumetric Rate'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(3,M,NCW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',VAR_CW(3,M,NCW)
          IUNM = 3
          IUNS = -1
          CALL RDUNIT(UNTS,VAR_CW(3,M,NCW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(3,M,NCW),
     &      ', m^3/s)'
!
!---      Well fluid density  ---
!
          VARB = 'Well Fluid Density'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(4,M,NCW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',VAR_CW(4,M,NCW)
          IUNKG = 1
          IUNM = -3
          CALL RDUNIT(UNTS,VAR_CW(4,M,NCW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(4,M,NCW),
     &      ', kg/m^3)'
!
!---      Well fluid specific heat  ---
!
          VARB = 'Well Fluid Specific Heat'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_CW(5,M,NCW))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',VAR_CW(5,M,NCW)
          IUNM = 2
          IUNS = -2
          IUNK = -1
          CALL RDUNIT(UNTS,VAR_CW(5,M,NCW),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_CW(5,M,NCW),
     &      ', J/kg K)'
  200   CONTINUE
  400 CONTINUE
      WRITE(IWR,'(/)')
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RD_GRLP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WR_GRLP_WELL
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
!     STOMP-GT
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
!     N_GLW - number of ground-loop wells
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
      SUB_LOG(ISUB_LOG) = '/WR_GRLP_WELL'
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
!---  Loop over number of ground-loop wells  ---
!
      DO 400 NCW = 1,N_GLW
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
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WR_GRLP_WELL group  ---
!
      RETURN
      END



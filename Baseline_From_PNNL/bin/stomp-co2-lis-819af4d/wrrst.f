!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRRST
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
!     Write restart files
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 15, 1997.




!     wrrst.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PARM_BH
      USE LEAK_WELL
      USE HYST
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE GEO_MECH
      USE FILES
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP_BH
      USE FDVP
      USE FDVH
      USE FDVGC_FRC
      USE FDVGC
      USE FDVD
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*32 FN
      CHARACTER*4 FORM1
      CHARACTER*19 FORM2
      CHARACTER*37 FORM3
      CHARACTER*38 FORM4
      CHARACTER*20 FORM5
      CHARACTER*38 FORM6
      CHARACTER*39 FORM7
      CHARACTER*22 FORM8
      CHARACTER*40 FORM9
      CHARACTER*41 FORM10
      CHARACTER*16 FORM11
      CHARACTER*17 FORM12
      CHARACTER*39 FORM13
      CHARACTER*40 FORM14
      CHARACTER*42 FORM15




!      REAL*8 SP_CX(LSPR)

!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM1,FORM2,FORM3,FORM4,FORM5,FORM6,FORM7
      SAVE FORM8,FORM9,FORM10,FORM11,FORM12,FORM13
      SAVE FORM14,FORM15
      DATA FORM1 /'(I6)'/
      DATA FORM2 /'(1(1PE22.15,1X),I3)'/
      DATA FORM3 /'(1(1PE22.15,1X),I3,1X,1(1PE22.15,1X))'/
      DATA FORM4 /'(1(1PE22.15,1X),I3,1X,10(1PE22.15,1X))'/
      DATA FORM5 /'(10(1PE22.15,1X),I3)'/
      DATA FORM6 /'(10(1PE22.15,1X),I3,1X,1(1PE22.15,1X))'/
      DATA FORM7 /'(10(1PE22.15,1X),I3,1X,10(1PE22.15,1X))'/
      DATA FORM8 /'(1(1PE22.15,1X),2(I3))'/
      DATA FORM9 /'(1(1PE22.15,1X),2(I3),1X,1(1PE22.15,1X))'/
      DATA FORM10 /'(1(1PE22.15,1X),2(I3),1X,10(1PE22.15,1X))'/
      DATA FORM11 /'(1(1PE22.15,1X))'/
      DATA FORM12 /'(10(1PE22.15,1X))'/
      DATA FORM13 /'(1(1PE22.15,1X),I3,1X,100(1PE22.15,1X))'/
      DATA FORM14 /'(10(1PE22.15,1X),I3,1X,100(1PE22.15,1X))'/
      DATA FORM15 /'(1(1PE22.15,1X),2(I3),1X,100(1PE22.15,1X))'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRRST'
!
!---  Create a new restart file with number of time steps
!     as the file name extension  ---
!
      FN(1:8) = 'restart.'
      N1X = MIN( ICOUNT( MXSTEP+NRST ),9 )
      DO 10 N2X = 1,N1X
        N3X = N2X + 8
        FN(N3X:N3X) = '0'
   10 CONTINUE
      N4X = ICOUNT( NSTEP )
      WRITE(FORM1(3:3),'(I1)') N4X
      N5X = 9 + N1X - N4X
      WRITE( FN(N5X:),FORM1) NSTEP
      OPEN(UNIT=IRS, FILE=FN, FORM='FORMATTED')
      CLOSE(UNIT=IRS, STATUS='DELETE')
      OPEN(UNIT=IRS, FILE=FN, STATUS='NEW', FORM='FORMATTED')
!
!---  Create a new restart_bf file with number of time steps
!     as the file name extension  ---
!
      IF( ISLC(74).NE.0 ) THEN
        FN(1:11) = 'restart_bf.'
        N1X = MIN( ICOUNT( MXSTEP+NRST ),9 )
        DO N2X = 1,N1X
          N3X = N2X + 11
          FN(N3X:N3X) = '0'
        ENDDO
        N4X = ICOUNT( NSTEP )
        WRITE(FORM1(3:3),'(I1)') N4X
        N5X = 12 + N1X - N4X
        WRITE( FN(N5X:),FORM1) NSTEP
        OPEN(UNIT=IRS_BF, FILE=FN, FORM='FORMATTED')
        CLOSE(UNIT=IRS_BF, STATUS='DELETE')
        OPEN(UNIT=IRS_BF, FILE=FN, STATUS='NEW', FORM='FORMATTED')
      ENDIF

!
!---  Write header  ---
!
      WRITE(IRS,'(A,//)')' Welcome to ...'
      WRITE(IRS,'(A)')   '                           STOMP'
      WRITE(IRS,'(A,//)')'        Subsurface Transport Over Multiple Pha
     &ses'
      WRITE(IRS,'(A)')   ' This file was produced by STOMP, a numerical
     &simulator'
      WRITE(IRS,'(A)')   ' developed by the Pacific Northwest Laboratory
     &, with'
      WRITE(IRS,'(A)')   ' support from the VOC-Arid Integrated Demonstr
     &ation Project,'
      WRITE(IRS,'(A)')   ' Office of Technology Development, U.S. Depart
     &ment of Energy.'
      WRITE(IRS,'(A)')   ' Results from this version of STOMP should not
     & be used for'
      WRITE(IRS,'(A,/)') ' license related applications.'
      WRITE(IRS,'(A,/)') ' For inquiries or assistance:  Call (509) 372-
     &6070'
      WRITE(IRS,'(A,//)')'                       ---  RESTART  ---'

      WRITE(IRS,'(2A)') 'Version: ',CH_VRSN




!      WRITE(IRS,'(2A)') 'Date: ',CHDATE
!      WRITE(IRS,'(2A)') 'Time: ',CHTIME
      WRITE(IRS,'(A)') 'Date: '
      WRITE(IRS,'(A)') 'Time: '

!
!---  Write timing data, field data by node numbers  ---
!

      NSPR = NSPL + NSPS + NSPE + NSPG + NSPN
      WRITE(IRS,'(/,7(1PE22.15),10(I9))') TM,DT,DTMX,DTMN,DTAF,
     &  RSDMX,DTCF,NRIMX,NSTEP,NFLD,NSOLU,IOM,ISLC(5),NSPR,NSPS,
     &  NFBN,NWN_LW






!
!---  Fracture flow data  ---
!
      IF( ISLC(74).NE.0 ) THEN
!
!---    Write header  ---
!
        WRITE(IRS_BF,'(A,//)')' Welcome to ...'
        WRITE(IRS_BF,'(A)')   '                           STOMP'
        WRITE(IRS_BF,'(A,//)')'        Subsurface Transport Over ' // 
     &    'Multiple Phases'
        WRITE(IRS_BF,'(A)')   ' This file was produced by STOMP, a ' // 
     &    'numerical simulator'
        WRITE(IRS_BF,'(A)')   ' developed by the Pacific Northwest ' // 
     &    'Laboratory  , with'
        WRITE(IRS_BF,'(A)')   ' support from the VOC-Arid ' // 
     &    'Integrated Demonstration Project,'
        WRITE(IRS_BF,'(A)')   ' Office of Technology Development, ' // 
     &    'U.S. Department of Energy.'
        WRITE(IRS_BF,'(A)')   ' Results from this version of STOMP ' // 
     &    'should not be used for'
        WRITE(IRS_BF,'(A,/)') ' license related applications.'
        WRITE(IRS_BF,'(A,/)') ' For inquiries or assistance:  ' // 
     &    'Call (509) 372-6070'
        WRITE(IRS_BF,'(A,//)')'                      ---  ' // 
     &    'RESTART_BF  ---'

        WRITE(IRS_BF,'(2A)') 'Version: ',CH_VRSN
        WRITE(IRS_BF,'(A)') 'Date: '
        WRITE(IRS_BF,'(A)') 'Time: '
!
!---    Write fracture geometry data  ---
!
        WRITE(IRS_BF,'(/,I9)') NF_FRC
        DO NFX = 1,NF_FRC
          WRITE(IRS_BF,'(2(I9))') IP_FRC(1,NFX),IP_FRC(2,NFX)
        ENDDO
!
!---    Write borehole geometry data  ---
!
        WRITE(IRS_BF,'(I9)') N_BH
        DO NBH = 1,N_BH
          WRITE(IRS_BF,'(2(I9))') ID_BH(3,NBH),ID_BH(4,NBH)
        ENDDO
      ENDIF
!
!---  Water Operational Mode ---
!
      IF( IOM.EQ.1 ) THEN
        NRSV = 8
        IF( (NSOLU+NSPR+NSPS).LE.0 ) THEN
          WRITE( FORM8(2:2),'(I1)' ) NRSV
          DO 110 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM8) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),IPH(2,N)
  110     CONTINUE
        ELSEIF( (NSOLU+NSPR+NSPS).LT.10 ) THEN
          WRITE( FORM9(2:2),'(I1)' ) NRSV
          WRITE( FORM9(26:26),'(I1)' ) NSOLU+NSPR+NSPS
          DO 114 N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM9) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),IPH(2,N),
     &        (C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)





  114     CONTINUE
        ELSEIF( (NSOLU+NSPR+NSPS).LT.100 ) THEN
          WRITE( FORM10(2:2),'(I1)' ) NRSV
          WRITE( FORM10(26:27),'(I2)' ) NSOLU+NSPR+NSPS
          DO 117 N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM10) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),IPH(2,N),
     &        (C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)





  117     CONTINUE
        ELSE
          WRITE( FORM15(2:2),'(I1)' ) NRSV
          WRITE( FORM15(26:28),'(I3)' ) NSOLU+NSPR+NSPS
          DO 119 N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM15) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),IPH(2,N),
     &        (C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)





  119     CONTINUE
        ENDIF
!
!---  STOMP-GT Operational Mode ---
!
      ELSEIF( IOM.EQ.3 ) THEN
        NRSV = 14
        IF( NSOLU.LE.0 ) THEN
          WRITE( FORM5(2:3),'(I2)' ) NRSV
          DO 130 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM5) T(2,N),PL(2,N),PG(2,N),PVW(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SG(2,N),PVA(2,N),XLA(2,N),YLS(2,N),
     &        TMS(2,N),PCMP(N),TCMP(N),NPHAZ(2,N)
  130     CONTINUE
!
!---      Fracture flow data  ---
!
          IF( ISLC(74).NE.0 ) THEN
            NRSV_BF = 10
            WRITE( FORM5(2:3),'(I2)' ) NRSV_BF
!
!---       Loop over fractures  ---
!
            DO NFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
                WRITE(IRS_BF,FORM5) T_FRC(2,NTX),PL_FRC(2,NTX),
     &            PG_FRC(2,NTX),PVW_FRC(2,NTX),SL_FRC(2,NTX),
     &            SG_FRC(2,NTX),PVA_FRC(2,NTX),XLA_FRC(2,NTX),
     &            YLS_FRC(2,NTX),TMS_FRC(2,NTX),NPHAZ_FRC(2,NTX)
              ENDDO
            ENDDO
!
!---        Loop over boreholes  ---
!
            DO NBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
                WRITE(IRS_BF,FORM5) T_BH(2,NBN),PL_BH(2,NBN),
     &            PG_BH(2,NBN),PVW_BH(2,NBN),SL_BH(2,NBN),
     &            SG_BH(2,NBN),PVA_BH(2,NBN),XLA_BH(2,NBN),
     &            YLS_BH(2,NBN),TMS_BH(2,NBN),NPHAZ_BH(2,NBN)
              ENDDO
            ENDDO
          ENDIF
        ELSEIF( NSOLU.LT.10 ) THEN
          WRITE( FORM6(2:3),'(I2)' ) NRSV
          WRITE( FORM6(24:24),'(I1)' ) NSOLU
          DO 132 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM6) T(2,N),PL(2,N),PG(2,N),PVW(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SG(2,N),PVA(2,N),XLA(2,N),YLS(2,N),
     &        TMS(2,N),PCMP(N),TCMP(N),
     &        NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU)
  132     CONTINUE
!
!---      Fracture flow data  ---
!
          IF( ISLC(74).NE.0 ) THEN
            NRSV_BF = 10
            WRITE( FORM6(2:3),'(I2)' ) NRSV_BF
            WRITE( FORM6(24:24),'(I1)' ) NSOLU
!
!---       Loop over fractures  ---
!
            DO NFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
                WRITE(IRS_BF,FORM6) T_FRC(2,NTX),PL_FRC(2,NTX),
     &            PG_FRC(2,NTX),PVW_FRC(2,NTX),SL_FRC(2,NTX),
     &            SG_FRC(2,NTX),PVA_FRC(2,NTX),XLA_FRC(2,NTX),
     &            YLS_FRC(2,NTX),TMS_FRC(2,NTX),NPHAZ_FRC(2,NTX),
     &            (C_FRC(NTX,NSL),NSL=1,NSOLU)
              ENDDO
            ENDDO
!
!---        Loop over boreholes  ---
!
            DO NBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
                WRITE(IRS_BF,FORM6) T_BH(2,NBN),PL_BH(2,NBN),
     &            PG_BH(2,NBN),PVW_BH(2,NBN),SL_BH(2,NBN),
     &            SG_BH(2,NBN),PVA_BH(2,NBN),XLA_BH(2,NBN),
     &            YLS_BH(2,NBN),TMS_BH(2,NBN),NPHAZ_BH(2,NBN),
     &            (C_BH(NBN,NSL),NSL=1,NSOLU)
              ENDDO
            ENDDO
          ENDIF
        ELSE
          WRITE( FORM7(2:3),'(I2)' ) NRSV
          WRITE( FORM7(24:25),'(I2)' ) NSOLU
          DO 134 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM7) T(2,N),PL(2,N),PG(2,N),PVW(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SG(2,N),PVA(2,N),XLA(2,N),YLS(2,N),
     &        TMS(2,N),PCMP(N),TCMP(N),
     &        NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU)
  134     CONTINUE
!
!---      Fracture flow data  ---
!
          IF( ISLC(74).NE.0 ) THEN
            NRSV_BF = 10
            WRITE( FORM7(2:3),'(I2)' ) NRSV_BF
            WRITE( FORM7(24:25),'(I2)' ) NSOLU
!
!---       Loop over fractures  ---
!
            DO NFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
                WRITE(IRS_BF,FORM7) T_FRC(2,NTX),PL_FRC(2,NTX),
     &            PG_FRC(2,NTX),PVW_FRC(2,NTX),SL_FRC(2,NTX),
     &            SG_FRC(2,NTX),PVA_FRC(2,NTX),XLA_FRC(2,NTX),
     &            YLS_FRC(2,NTX),TMS_FRC(2,NTX),NPHAZ_FRC(2,NTX),
     &            (C_FRC(NTX,NSL),NSL=1,NSOLU)
              ENDDO
            ENDDO
!
!---       Loop over boreholes  ---
!
            DO NBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
                WRITE(IRS_BF,FORM7) T_BH(2,NBN),PL_BH(2,NBN),
     &            PG_BH(2,NBN),PVW_BH(2,NBN),SL_BH(2,NBN),
     &            SG_BH(2,NBN),PVA_BH(2,NBN),XLA_BH(2,NBN),
     &            YLS_BH(2,NBN),TMS_BH(2,NBN),NPHAZ_BH(2,NBN),
     &            (C_BH(NBN,NSL),NSL=1,NSOLU)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
!
!---  Water-Oil Operational Mode ---
!
      ELSEIF( IOM.EQ.4 ) THEN
        NRSV = 17
        IF( NSOLU.LE.0 ) THEN
          WRITE( FORM5(2:3),'(I2)' ) NRSV
          DO 140 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM5) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SNT(2,N),SN(2,N),XMLO(2,N),
     &        PN(2,N),SGTL(N),SGTN(N),ASTMIN(2,N),TRPNL(2,N),
     &        ASTMAX(2,N),PCMP(N),TCMP(N),NPHAZ(2,N)
  140     CONTINUE
        ELSEIF( NSOLU.LT.10 ) THEN
          WRITE( FORM6(2:3),'(I2)' ) NRSV
          WRITE( FORM6(24:24),'(I1)' ) NSOLU
          DO 142 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM6) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SNT(2,N),SN(2,N),XMLO(2,N),
     &        PN(2,N),SGTL(N),SGTN(N),ASTMIN(2,N),TRPNL(2,N),
     &        ASTMAX(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),
     &        (C(N,NSL),NSL=1,NSOLU)
  142     CONTINUE
        ELSE
          WRITE( FORM7(2:3),'(I2)' ) NRSV
          WRITE( FORM7(24:25),'(I2)' ) NSOLU
          DO 144 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM7) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SNT(2,N),SN(2,N),XMLO(2,N),
     &        PN(2,N),SGTL(N),SGTN(N),ASTMIN(2,N),TRPNL(2,N),
     &        ASTMAX(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),
     &        (C(N,NSL),NSL=1,NSOLU)
  144     CONTINUE
        ENDIF
!
!---  Water-Oil-Air Operational Mode ---
!
      ELSEIF( IOM.EQ.5 ) THEN
        NRSV = 18
        IF( NSOLU.LE.0 ) THEN
          WRITE( FORM5(2:3),'(I2)' ) NRSV
          DO 150 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM5) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SNT(2,N),SN(2,N),XMLO(2,N),
     &        XMLA(2,N),PN(2,N),SGTL(N),SGTN(N),ASTMIN(2,N),
     &        TRPNL(2,N),ASTMAX(2,N),PCMP(N),TCMP(N),NPHAZ(2,N)
  150     CONTINUE
        ELSEIF( NSOLU.LT.10 ) THEN
          WRITE( FORM6(2:3),'(I2)' ) NRSV
          WRITE( FORM6(24:24),'(I1)' ) NSOLU
          DO 152 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM6) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SNT(2,N),SN(2,N),XMLO(2,N),
     &        XMLA(2,N),PN(2,N),SGTL(N),SGTN(N),ASTMIN(2,N),
     &        TRPNL(2,N),ASTMAX(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),
     &        (C(N,NSL),NSL=1,NSOLU)
  152     CONTINUE
        ELSE
          WRITE( FORM7(2:3),'(I2)' ) NRSV
          WRITE( FORM7(24:25),'(I2)' ) NSOLU
          DO 154 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM7) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SNT(2,N),SN(2,N),XMLO(2,N),
     &        XMLA(2,N), PN(2,N),SGTL(N),SGTN(N),ASTMIN(2,N),
     &        TRPNL(2,N),ASTMAX(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),
     &        (C(N,NSL),NSL=1,NSOLU)
  154     CONTINUE
        ENDIF
!
!---  STOMP-CO2 and STOMP-CO2E Operational Modes ---
!
      ELSEIF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
        NRSV = 13
        IF( (NSOLU+NSPR+NSPS).LE.0 ) THEN
          WRITE( FORM5(2:3),'(I2)' ) NRSV
          DO 320 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM5) T(2,N),PL(2,N),PG(2,N),PVW(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SG(2,N),XLA(2,N),YLS(2,N),TMS(2,N),
     &        PCMP(N),TCMP(N),NPHAZ(2,N)
  320     CONTINUE
        ELSEIF( (NSOLU+NSPR+NSPS).LT.10 ) THEN
          WRITE( FORM6(2:3),'(I2)' ) NRSV
          WRITE( FORM6(24:24),'(I1)' ) NSOLU+NSPR+NSPS
          DO 324 N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM6) T(2,N),PL(2,N),PG(2,N),PVW(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SG(2,N),XLA(2,N),YLS(2,N),TMS(2,N),
     &        PCMP(N),TCMP(N),NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)





  324     CONTINUE
        ELSEIF( (NSOLU+NSPR+NSPS).LT.100 ) THEN
          WRITE( FORM7(2:3),'(I2)' ) NRSV
          WRITE( FORM7(24:25),'(I2)' ) NSOLU+NSPR+NSPS
          DO 327 N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM7) T(2,N),PL(2,N),PG(2,N),PVW(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SG(2,N),XLA(2,N),YLS(2,N),TMS(2,N),
     &        PCMP(N),TCMP(N),NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)





  327     CONTINUE
        ELSE
          WRITE( FORM14(2:3),'(I2)' ) NRSV
          WRITE( FORM14(23:25),'(I3)' ) NSOLU+NSPR+NSPS
          DO 329 N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM14) T(2,N),PL(2,N),PG(2,N),PVW(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SG(2,N),XLA(2,N),YLS(2,N),TMS(2,N),
     &        PCMP(N),TCMP(N),NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)





  329     CONTINUE
        ENDIF
!
!---  STOMP-SEQ Operational Mode ---
!
      ELSEIF( IOM.EQ.34 ) THEN
        NRSV = 15
        IF( (NSOLU+NSPR+NSPS).LE.0 ) THEN
          WRITE( FORM5(2:3),'(I2)' ) NRSV
          DO N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM5) T(2,N),PL(2,N),PG(2,N),PSO(2,N),PVW(2,N),
     &        SL(2,N),SGT(2,N),ASLMIN(2,N),SG(2,N),XLA(2,N),XMLO(2,N),
     &        YLS(2,N),TMS(2,N),PCMP(N),TCMP(N),NPHAZ(2,N)
          ENDDO
        ELSEIF( (NSOLU+NSPR+NSPS).LT.10 ) THEN
          WRITE( FORM6(2:3),'(I2)' ) NRSV
          WRITE( FORM6(24:24),'(I1)' ) NSOLU+NSPR+NSPS
          DO N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM6) T(2,N),PL(2,N),PG(2,N),PSO(2,N),PVW(2,N),
     &        SL(2,N),SGT(2,N),ASLMIN(2,N),SG(2,N),XLA(2,N),XMLO(2,N),
     &        YLS(2,N),TMS(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),
     &        (C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)






          ENDDO
        ELSEIF( (NSOLU+NSPR+NSPS).LT.100 ) THEN
          WRITE( FORM7(2:3),'(I2)' ) NRSV
          WRITE( FORM7(24:25),'(I2)' ) NSOLU+NSPR+NSPS
          DO N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM7) T(2,N),PL(2,N),PG(2,N),PSO(2,N),PVW(2,N),
     &        SL(2,N),SGT(2,N),ASLMIN(2,N),SG(2,N),XLA(2,N),XMLO(2,N),
     &        YLS(2,N),TMS(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),
     &        (C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)






          ENDDO
        ELSE
          WRITE( FORM14(2:3),'(I2)' ) NRSV
          WRITE( FORM14(23:25),'(I3)' ) NSOLU+NSPR+NSPS
          DO N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM14) T(2,N),PL(2,N),PG(2,N),PSO(2,N),PVW(2,N),
     &        SL(2,N),SGT(2,N),ASLMIN(2,N),SG(2,N),XLA(2,N),XMLO(2,N),
     &        YLS(2,N),TMS(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),
     &        (C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)     






          ENDDO
        ENDIF
!
!---  HYD Operational Modes ---
!
      ELSEIF( IOM.EQ.37 ) THEN
        NRSV = 14
        IF( NSOLU.LE.0 ) THEN
          WRITE( FORM5(2:3),'(I2)' ) NRSV
          DO 360 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM5) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SG(2,N),SH(2,N),YLS(2,N),PVA(2,N),PVO(2,N),
     &        SN(2,N),PN(2,N),YMGO(2,N),
     &        PCMP(N),TCMP(N),NPHAZ(2,N)
  360     CONTINUE
        ELSEIF( NSOLU.LT.10 ) THEN
          WRITE( FORM6(2:3),'(I2)' ) NRSV
          WRITE( FORM6(24:24),'(I1)' ) NSOLU
          DO 362 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM6) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SG(2,N),SH(2,N),YLS(2,N),PVA(2,N),PVO(2,N),
     &        SN(2,N),PN(2,N),YMGO(2,N),
     &        PCMP(N),TCMP(N),NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU)
  362     CONTINUE
        ELSE
          WRITE( FORM7(2:3),'(I2)' ) NRSV
          WRITE( FORM7(24:25),'(I2)' ) NSOLU
          DO 364 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM7) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SG(2,N),SH(2,N),YLS(2,N),PVA(2,N),PVO(2,N),
     &        SN(2,N),PN(2,N),YMGO(2,N),
     &        PCMP(N),TCMP(N),NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU)
  364     CONTINUE
        ENDIF
!
!---  HYD-KE Operational Mode ---
!
      ELSEIF( IOM.EQ.38 ) THEN
        NRSV = 19
        IF( NSOLU.LE.0 ) THEN
          WRITE( FORM5(2:3),'(I2)' ) NRSV
          DO 380 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM5) T(2,N),PL(2,N),PG(2,N),PN(2,N),
     &        PVA(2,N),PVO(2,N),SG(2,N),SL(2,N),SH(2,N),SN(2,N),
     &        YLS(2,N),PVHA(2,N),PVHO(2,N),YMGO(2,N),YMHGO(2,N),
     &        TMHA(2,N),TMHO(2,N),PCMP(N),TCMP(N),NPHAZ(2,N)
  380     CONTINUE
        ELSEIF( NSOLU.LT.10 ) THEN
          WRITE( FORM6(2:3),'(I2)' ) NRSV
          WRITE( FORM6(24:24),'(I1)' ) NSOLU
          DO 382 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM6) T(2,N),PL(2,N),PG(2,N),PN(2,N),
     &        PVA(2,N),PVO(2,N),SG(2,N),SL(2,N),SH(2,N),SN(2,N),
     &        YLS(2,N),PVHA(2,N),PVHO(2,N),YMGO(2,N),YMHGO(2,N),
     &        TMHA(2,N),TMHO(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),
     &        (C(N,NSL),NSL=1,NSOLU)
  382     CONTINUE
        ELSE
          WRITE( FORM7(2:3),'(I2)' ) NRSV
          WRITE( FORM7(24:25),'(I2)' ) NSOLU
          DO 384 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM7) T(2,N),PL(2,N),PG(2,N),PN(2,N),
     &        PVA(2,N),PVO(2,N),SG(2,N),SL(2,N),SH(2,N),SN(2,N),
     &        YLS(2,N),PVHA(2,N),PVHO(2,N),YMGO(2,N),YMHGO(2,N),
     &        TMHA(2,N),TMHO(2,N),PCMP(N),TCMP(N),NPHAZ(2,N),
     &        (C(N,NSL),NSL=1,NSOLU)
  384     CONTINUE
        ENDIF
!
!---  STOMP-HYDT-KE Operational Mode ---
!
      ELSEIF( IOM.EQ.39 ) THEN
        NRSV = 25
        IF( NSOLU.LE.0 ) THEN
          WRITE( FORM5(2:3),'(I2)' ) NRSV
          DO 390 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM5) T(2,N),PL(2,N),PG(2,N),PN(2,N),PSO(2,N),
     &        PVA(2,N),PVO(2,N),PVN(2,N),
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),SG(2,N),SL(2,N),SN(2,N),
     &        SH(2,N),SI(2,N),YLS(2,N),TMHA(2,N),TMHO(2,N),TMHN(2,N),
     &        YMHGA(2,N),YMHGO(2,N),YMHGN(2,N),PCMP(N),TCMP(N),
     &        NPHAZ(2,N)
  390     CONTINUE
        ELSEIF( NSOLU.LT.10 ) THEN
          WRITE( FORM6(2:3),'(I2)' ) NRSV
          WRITE( FORM6(24:24),'(I1)' ) NSOLU
          DO 392 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM6) T(2,N),PL(2,N),PG(2,N),PN(2,N),PSO(2,N),
     &        PVA(2,N),PVO(2,N),PVN(2,N),
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),SG(2,N),SL(2,N),SN(2,N),
     &        SH(2,N),SI(2,N),YLS(2,N),TMHA(2,N),TMHO(2,N),TMHN(2,N),
     &        YMHGA(2,N),YMHGO(2,N),YMHGN(2,N),PCMP(N),TCMP(N),
     &        NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU)
  392     CONTINUE
        ELSE
          WRITE( FORM7(2:3),'(I2)' ) NRSV
          WRITE( FORM7(24:25),'(I2)' ) NSOLU
          DO 394 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM7) T(2,N),PL(2,N),PG(2,N),PN(2,N),PSO(2,N),
     &        PVA(2,N),PVO(2,N),PVN(2,N),
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),SG(2,N),SL(2,N),SN(2,N),
     &        SH(2,N),SI(2,N),YLS(2,N),TMHA(2,N),TMHO(2,N),TMHN(2,N),
     &        YMHGA(2,N),YMHGO(2,N),YMHGN(2,N),PCMP(N),TCMP(N),
     &        NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU)
  394     CONTINUE
        ENDIF
!
!---  H2O-NaCl-NComponent-Energy Operational Mode ---
!
      ELSEIF( IOM.EQ.40 ) THEN
        NRSV = 10+NGC*2
        IF( (NSOLU+NSPR+NSPS).LE.0 ) THEN
          WRITE( FORM5(2:3),'(I2)' ) NRSV
          DO 400 N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM5) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SG(2,N),(PVC(IGC,2,N),IGC=1,NGC),
     &        (XLC(IGC,2,N),IGC=1,NGC),YLS(2,N),PCMP(N),TCMP(N),
     &        NPHAZ(2,N)
  400     CONTINUE
        ELSEIF( (NSOLU+NSPR+NSPS).LT.10 ) THEN
          WRITE( FORM6(2:3),'(I2)' ) NRSV
          WRITE( FORM6(24:24),'(I1)' ) NSOLU+NSPR+NSPS
          DO 404 N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM6) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SG(2,N),(PVC(IGC,2,N),IGC=1,NGC),
     &        (XLC(IGC,2,N),IGC=1,NGC),YLS(2,N),
     &        PCMP(N),TCMP(N),NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)






  404     CONTINUE
        ELSEIF( (NSOLU+NSPR+NSPS).LT.100 ) THEN
          WRITE( FORM7(2:3),'(I2)' ) NRSV
          WRITE( FORM7(24:25),'(I2)' ) NSOLU+NSPR+NSPS
          DO 407 N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM7) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SG(2,N),(PVC(IGC,2,N),IGC=1,NGC),
     &        (XLC(IGC,2,N),IGC=1,NGC),YLS(2,N),
     &        PCMP(N),TCMP(N),NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)






  407     CONTINUE
        ELSE
          WRITE( FORM14(2:3),'(I2)' ) NRSV
          WRITE( FORM14(24:26),'(I3)' ) NSOLU+NSPR+NSPS
          DO 409 N = 1,NFBN+NWN_LW

!
!---        Species concentrations  ---
!
            WRITE(IRS,FORM14) T(2,N),PL(2,N),PG(2,N),SL(2,N),
     &        SGT(2,N),ASLMIN(2,N),SG(2,N),(PVC(IGC,2,N),IGC=1,NGC),
     &        (XLC(IGC,2,N),IGC=1,NGC),YLS(2,N),
     &        PCMP(N),TCMP(N),NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)






  409     CONTINUE
        ENDIF
!
!---  EOR Operational Mode ---
!
      ELSEIF( IOM.EQ.43 ) THEN
        NRSV = 14 + 2*(NGC+2)
        IF( (NSOLU+NSPR+NSPS).LE.0 ) THEN
          WRITE( FORM5(2:3),'(I2)' ) NRSV
          DO N = 1,NFBN+NWN_LW
            WRITE(IRS,FORM5) T(2,N),PG(2,N),PL(2,N),PN(2,N),POSM(2,N),
     &        PSO(2,N),PVA(2,N),SG(2,N),SL(2,N),SN(2,N),YLS(2,N),
     &        TMS(2,N),(TMC(IGC,2,N),IGC=1,NGC+2),
     &        (ZMC(IGC,2,N),IGC=1,NGC+2),PCMP(N),TCMP(N),NPHAZ(2,N)
          ENDDO
!
!---      Fault flow data  ---
!
          IF( ISLC(74).NE.0 ) THEN
            NRSV_BF = 14 + 2*(NGC+2)
            WRITE( FORM5(2:3),'(I2)' ) NRSV_BF
!
!---       Loop over fault  ---
!
            DO NFX = 1,NF_FRC
!
!---          Loop over fault triangles  ---
!
              DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
                WRITE(IRS_BF,FORM5) T_FRC(2,NTX),PG_FRC(2,NTX),
     &            PL_FRC(2,NTX),PN_FRC(2,NTX),POSM_FRC(2,NTX),
     &            PSO_FRC(2,NTX),PVA_FRC(2,NTX),SG_FRC(2,NTX),
     &            SL_FRC(2,NTX),SN_FRC(2,NTX),YLS_FRC(2,NTX),
     &            TMS_FRC(2,NTX),(TMC(IGC,2,N),IGC=1,NGC+2),
     &            (ZMC_FRC(IGC,2,NTX),IGC=1,NGC+2),PCMP_FRC(NTX),
     &            TCMP_FRC(NTX),NPHAZ_FRC(2,NTX)
              ENDDO
            ENDDO
          ENDIF
        ELSEIF( (NSOLU+NSPR+NSPS).LT.10 ) THEN
          WRITE( FORM6(2:3),'(I2)' ) NRSV
          WRITE( FORM6(24:24),'(I1)' ) NSOLU+NSPR+NSPS
          DO N = 1,NFBN+NWN_LW

!
!---        Solute and species concentrations  ---
!
            WRITE(IRS,FORM6) T(2,N),PG(2,N),PL(2,N),PN(2,N),POSM(2,N),
     &        PSO(2,N),
     &        PVA(2,N),SG(2,N),SL(2,N),SN(2,N),YLS(2,N),TMS(2,N),
     &        (TMC(IGC,2,N),IGC=1,NGC+2),(ZMC(IGC,2,N),IGC=1,NGC+2),
     &        PCMP(N),TCMP(N),NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)

          ENDDO
!
!---      Fault flow data  ---
!
          IF( ISLC(74).NE.0 ) THEN
            NRSV_BF = 14 + 2*(NGC+2)
            WRITE( FORM6(2:3),'(I2)' ) NRSV_BF
            WRITE( FORM6(24:24),'(I1)' ) NSOLU+NSPR
!
!---        Loop over faults  ---
!
            DO NFX = 1,NF_FRC
!
!---          Loop over fault triangles  ---
!
              DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)

!
!---            Solute and species concentrations  ---
!
                WRITE(IRS_BF,FORM6) T_FRC(2,NTX),PG_FRC(2,NTX),
     &            PL_FRC(2,NTX),PN_FRC(2,NTX),POSM_FRC(2,NTX),
     &            PSO_FRC(2,NTX),PVA_FRC(2,NTX),SG_FRC(2,NTX),
     &            SL_FRC(2,NTX),SN_FRC(2,NTX),YLS_FRC(2,NTX),
     &            TMS_FRC(2,NTX),(TMC(IGC,2,N),IGC=1,NGC+2),
     &            (ZMC_FRC(IGC,2,NTX),IGC=1,NGC+2),PCMP_FRC(NTX),
     &            TCMP_FRC(NTX),NPHAZ_FRC(2,NTX),
     &            (C_FRC(NTX,NSL),NSL=1,NSOLU),
     &            (SP_C_FRC(NTX,NSP),NSP=1,NSPR)

              ENDDO
            ENDDO
          ENDIF
        ELSEIF( (NSOLU+NSPR+NSPS).LT.100 ) THEN
          WRITE( FORM7(2:3),'(I2)' ) NRSV
          WRITE( FORM7(24:25),'(I2)' ) NSOLU+NSPR+NSPS
          DO N = 1,NFBN+NWN_LW

!
!---        Solute and species concentrations  ---
!
            WRITE(IRS,FORM7) T(2,N),PG(2,N),PL(2,N),PN(2,N),POSM(2,N),
     &        PSO(2,N),
     &        PVA(2,N),SG(2,N),SL(2,N),SN(2,N),YLS(2,N),TMS(2,N),
     &        (TMC(IGC,2,N),IGC=1,NGC+2),(ZMC(IGC,2,N),IGC=1,NGC+2),
     &        PCMP(N),TCMP(N),NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)

          ENDDO
!
!---      Fault flow data  ---
!
          IF( ISLC(74).NE.0 ) THEN
            NRSV_BF = 14 + 2*(NGC+2)
            WRITE( FORM7(2:3),'(I2)' ) NRSV_BF
            WRITE( FORM7(24:25),'(I2)' ) NSOLU+NSPR
!
!---        Loop over faults  ---
!
            DO NFX = 1,NF_FRC
!
!---          Loop over fault triangles  ---
!
              DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)

!
!---            Solute and species concentrations  ---
!
                WRITE(IRS_BF,FORM7) T_FRC(2,NTX),PG_FRC(2,NTX),
     &            PL_FRC(2,NTX),PN_FRC(2,NTX),POSM_FRC(2,NTX),
     &            PSO_FRC(2,NTX),PVA_FRC(2,NTX),SG_FRC(2,NTX),
     &            SL_FRC(2,NTX),SN_FRC(2,NTX),YLS_FRC(2,NTX),
     &            TMS_FRC(2,NTX),(TMC(IGC,2,N),IGC=1,NGC+2),
     &            (ZMC_FRC(IGC,2,NTX),IGC=1,NGC+2),PCMP_FRC(NTX),
     &            TCMP_FRC(NTX),NPHAZ_FRC(2,NTX),
     &            (C_FRC(NTX,NSL),NSL=1,NSOLU),
     &            (SP_C_FRC(NTX,NSP),NSP=1,NSPR)

              ENDDO
            ENDDO
          ENDIF
        ELSE
          WRITE( FORM14(2:3),'(I2)' ) NRSV
          WRITE( FORM14(23:25),'(I3)' ) NSOLU+NSPR+NSPS
          DO N = 1,NFBN+NWN_LW

!
!---        Solute and species concentrations  ---
!
            WRITE(IRS,FORM14) T(2,N),PG(2,N),PL(2,N),PN(2,N),POSM(2,N),
     &        PSO(2,N),
     &        PVA(2,N),SG(2,N),SL(2,N),SN(2,N),YLS(2,N),TMS(2,N),
     &        (TMC(IGC,2,N),IGC=1,NGC+2),(ZMC(IGC,2,N),IGC=1,NGC+2),
     &        PCMP(N),TCMP(N),NPHAZ(2,N),(C(N,NSL),NSL=1,NSOLU),
     &        (SP_C(N,NSP),NSP=1,NSPR),(SP_CMN(N,NSP),NSP=1,NSPS)

          ENDDO
!
!---      Fault flow data  ---
!
          IF( ISLC(74).NE.0 ) THEN
            NRSV_BF = 14 + 2*(NGC+2)
            WRITE( FORM14(2:3),'(I2)' ) NRSV_BF
            WRITE( FORM14(23:25),'(I3)' ) NSOLU+NSPR
!
!---        Loop over faults  ---
!
            DO NFX = 1,NF_FRC
!
!---          Loop over fault triangles  ---
!
              DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)

!
!---            Solute and species concentrations  ---
!
                WRITE(IRS_BF,FORM14) T_FRC(2,NTX),PG_FRC(2,NTX),
     &            PL_FRC(2,NTX),PN_FRC(2,NTX),POSM_FRC(2,NTX),
     &            PSO_FRC(2,NTX),PVA_FRC(2,NTX),SG_FRC(2,NTX),
     &            SL_FRC(2,NTX),SN_FRC(2,NTX),YLS_FRC(2,NTX),
     &            TMS_FRC(2,NTX),(TMC(IGC,2,N),IGC=1,NGC+2),
     &            (ZMC_FRC(IGC,2,NTX),IGC=1,NGC+2),PCMP_FRC(NTX),
     &            TCMP_FRC(NTX),NPHAZ_FRC(2,NTX),
     &            (C_FRC(NTX,NSL),NSL=1,NSOLU),
     &            (SP_C_FRC(NTX,NSP),NSP=1,NSPR)

              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF
!
!---  Geomechanics data  ---
!
      IF( ISLC(50).NE.0 ) THEN
        WRITE(IRS,'(A)') 'Geomechanics Model Data'
        WRITE(FORM1(3:3),'(I1)') ICOUNT(NFEN_GM)
        WRITE(IRS,FORM1) NFEN_GM
        DO NFEN = 1,NFEN_GM
          WRITE(IRS,'(3(1PE22.15,1X))') U_GM(1,NFEN),V_GM(1,NFEN),
     &      W_GM(1,NFEN)
        ENDDO
      ENDIF
!
!---  Coupled-well data  ---
!
      IF( N_CW.GT.0 ) THEN
        WRITE(IRS,'(A)') 'Coupled-Well Model Data'
        WRITE(FORM1(3:3),'(I1)') ICOUNT(N_CW)
        WRITE(IRS,FORM1) N_CW
        DO NCW = 1,N_CW
          WRITE(IRS,'(1PE22.15)') P_CW(2,NCW)
        ENDDO
      ENDIF
!
!---  Close the restart files  ---
!
      CLOSE( UNIT=IRS )
      IF( ISLC(74).NE.0 ) CLOSE( UNIT=IRS_BF )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRRST group
!
      RETURN
      END

